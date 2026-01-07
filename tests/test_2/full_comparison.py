import os
import re
import shutil
import subprocess
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Any

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from tqdm import tqdm

# --- Constants ---
C_CGS = 2.99792458e10
G_CGS = 6.6730831e-8
SQRT_KAPPA = np.sqrt(1.0e-15 * C_CGS * C_CGS / G_CGS)
CU_TO_HZ = C_CGS / SQRT_KAPPA
CU_TO_KHZ = CU_TO_HZ / 1000.0

# --- Configuration Template ---
CONFIG_TEMPLATE = """# Configuration file for RNS
# WARNING: Parameters are auto-generated.

# SETTINGS
EoS_BM_type         0
EoS_BM_file         {eos_bm}
EoS_DM_type         0
EoS_DM_file         {eos_dm}
accuracy            1e-10

# PARAMETERS
BM_central_energy   {e_center}      # Central energy density 
DM_fraction         0.0             # Target dark matter fraction 

BM_rotation_type    2               # 2 -> Uryu 8
BM_A                0.7
BM_B                0.5
BM_lambda1          {lambda1}       # Omega_max / Omega_central
BM_lambda2          {lambda2}       # Omega_equator / Omega_central

r_ratio_BM          {r_ratio}       # Polar/Equatorial radius ratio
r_step              0.025

1Doutput            0
2Doutput            1
output_name         {out_name}
verbose             1
"""

def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='[%(levelname)s] %(message)s')

class ReferenceDataLoader:
    """Handles parsing of the reference .d data file."""
    
    COLUMNS = [
        "rp/re", "εc", "εmax", "M", "M0", "J", "T/|W|", 
        "Ωc", "Ωmax", "Ωe", "ΩK", "Re", "re", "GRV3"
    ]

    @staticmethod
    def load(filename: str) -> Dict[Tuple[float, float], Dict[str, List[float]]]:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Reference data file not found: {filename}")

        data = {}
        current_lambda = None
        
        logging.info(f"Loading reference data from {filename}...")
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Detect lambda pair definition: {l1, l2}
                lambda_match = re.match(r"\{([\d\.]+),\s*([\d\.]+)\}", line)
                if lambda_match:
                    current_lambda = tuple(map(float, lambda_match.groups()))
                    data[current_lambda] = {col: [] for col in ReferenceDataLoader.COLUMNS}
                else:
                    if current_lambda is None:
                        continue # specific logic to skip headers if format varies
                    
                    try:
                        values = list(map(float, line.split()))
                        if len(values) != len(ReferenceDataLoader.COLUMNS):
                            continue # Skip malformed lines
                            
                        for col, val in zip(ReferenceDataLoader.COLUMNS, values):
                            data[current_lambda][col].append(val)
                    except ValueError:
                        continue

        logging.info(f"Loaded {len(data)} configuration sets.")
        return data

class SimulationRunner:
    """Manages the generation of configs and execution of the C binary."""

    def __init__(self, executable: str, config_dir: str, output_dir: str, eos_paths: Dict[str, str]):
        self.exec_path = executable
        self.config_dir = config_dir
        self.output_dir = output_dir
        self.eos_paths = eos_paths
        
        os.makedirs(self.config_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_filename(self, l1, l2, ec, r_ratio):
        """Strict naming convention."""
        return f"1fluid_l1{l1}_l2{l2}_e{ec:.4f}_r{r_ratio:.4f}"

    def prepare_config(self, l1, l2, ec, r_ratio) -> Tuple[str, str]:
        fname = self.generate_filename(l1, l2, ec, r_ratio)
        config_content = CONFIG_TEMPLATE.format(
            e_center=ec,
            r_ratio=r_ratio,
            lambda1=l1,
            lambda2=l2,
            out_name=fname,
            eos_bm=self.eos_paths['bm'],
            eos_dm=self.eos_paths['dm']
        )
        
        config_path = os.path.join(self.config_dir, fname + ".d")
        with open(config_path, "w") as f:
            f.write(config_content)
            
        return fname, config_path

    @staticmethod
    def _execute_single_run(args):
        """Static method for ProcessPool serialization."""
        exec_path, fname, config_path, out_dir = args
        full_out_path = os.path.join(out_dir, fname + ".out")
        
        # Skip if output exists
        if os.path.isfile(full_out_path):
            return f"SKIP: {fname}"

        cmd = [exec_path, "-c", config_path]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            
            # Write stdout/stderr to log file
            with open(full_out_path, "w") as f:
                f.write("=== STDOUT ===\n")
                f.write(result.stdout)
                f.write("\n\n=== STDERR ===\n")
                f.write(result.stderr)

            if result.returncode != 0:
                return f"FAIL: {fname} (Exit Code {result.returncode})"

            # Post-process: Move files
            # 1. Move config
            shutil.move(config_path, os.path.join(out_dir, fname + ".d"))
            
            # 2. Move .dat file (RNS default output)
            dat_file = os.path.join(".", fname + ".dat")
            if os.path.exists(dat_file):
                shutil.move(dat_file, os.path.join(out_dir, fname + ".dat"))

            # 3. Find and move HDF5
            h5_path = None
            for line in result.stdout.splitlines():
                if "Saving H5 file in:" in line:
                    pass 
                if ".h5" in line and "/" in line:
                    possible = line.strip().split()[-1] 
                    if possible.endswith(".h5"):
                        h5_path = possible

            lines = result.stdout.splitlines()
            for i, line in enumerate(lines):
                if "Saving H5 file in:" in line and i + 1 < len(lines):
                     h5_path = lines[i + 1].strip()
                     break
            
            if h5_path and os.path.exists(h5_path):
                shutil.move(h5_path, os.path.join(out_dir, os.path.basename(h5_path)))
            else:
                return f"WARN: {fname} (HDF5 missing)"

            return f"DONE: {fname}"

        except Exception as e:
            return f"ERR: {fname} ({str(e)})"

    def run_all(self, tasks, max_workers=8):
        logging.info(f"Starting simulation of {len(tasks)} cases with {max_workers} workers.")
        
        # Prepare arguments for static method
        # task structure: (l1, l2, ec, r_ratio)
        exec_args = []
        for l1, l2, ec, r_ratio in tasks:
            fname, cfg_path = self.prepare_config(l1, l2, ec, r_ratio)
            exec_args.append((self.exec_path, fname, cfg_path, self.output_dir))

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._execute_single_run, arg) for arg in exec_args]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing", unit="run"):
                res = future.result()
                if "FAIL" in res or "ERR" in res:
                    logging.error(res)
                elif "WARN" in res:
                    logging.warning(res)
                else:
                    logging.debug(res)

class ResultAnalyzer:
    """Parses output files and compares with reference."""
    
    def __init__(self, out_dir: str):
        self.out_dir = out_dir

    def parse_output(self, filename: str) -> Dict[str, Any] | None:
        """Parses .out and .h5 files for a single run."""
        base_path = os.path.join(self.out_dir, filename)
        out_file = base_path + ".out"
        h5_file = base_path + ".h5"
        
        result = {}

        # 1. Parse Text Output (.out)
        try:
            with open(out_file, "r") as f:
                lines = f.readlines()
            
            # Search backwards for baryonic matter data
            found_baryonic = False
            for i in range(len(lines) - 1, -1, -1):
                line = lines[i].strip()
                
                if "Baryonic Matter:" in line:
                    if "Target axis ratio" in line:
                        continue

                    if i + 2 < len(lines):
                        target_line = lines[i + 2].strip()
                        parts = target_line.split()
                        
                        # Safety check: ensure we actually have numbers
                        try:
                            # Try converting the first element to verify it's a number line
                            float(parts[0]) 
                            
                            if len(parts) >= 4:
                                result["M"] = float(parts[0])
                                result["M0"] = float(parts[1])
                                result["Re"] = float(parts[2])
                                result["J/M^2"] = float(parts[3])
                                found_baryonic = True
                                break
                        except (ValueError, IndexError):
                            # If conversion fails, this wasn't the data line (e.g. empty or text)
                            continue
            
            if not found_baryonic:
                logging.warning(f"Could not find 'Baryonic Matter' data in {filename}")
                return None
                
        except FileNotFoundError:
            logging.warning(f"Output file not found: {out_file}")
            return None

        # 2. Parse HDF5
        try:
            with h5py.File(h5_file, "r") as hf:
                Omega_diffBM = hf["Omega_diffBM"][:]
                s_qp = hf["s_qp"][:]
                
                Omega_s = Omega_diffBM[:, 1] * CU_TO_KHZ
                
                # Interpolate
                interp = interp1d(s_qp, Omega_s, kind='linear', 
                                  bounds_error=False, fill_value="extrapolate")
                
                result["Ωc"] = float(interp(0))
                result["Ωe"] = float(interp(0.5))
                
                # Find Max
                res = minimize_scalar(lambda s: -interp(s), 
                                      bounds=(s_qp[0], s_qp[-1]), method='bounded')
                result["Ωmax"] = float(interp(res.x))
                
        except (FileNotFoundError, KeyError, OSError) as e:
            logging.warning(f"HDF5 Issue for {filename}: {e}")
            return None

        return result
    
    def compare(self, ref_data: Dict, runner: SimulationRunner):
        logging.info("Analyzing results...")
        errors = {}
        
        for k_pair, data_dict in ref_data.items():
            l1, l2 = k_pair
            count = len(data_dict["εc"])
            
            # Arrays to store mined data
            mined = {
                "M": [], "Re": [], "J/M^2": [], 
                "Ωc": [], "Ωe": [], "Ωmax": []
            }
            
            valid_indices = []

            for i in range(count):
                ec = data_dict["εc"][i]
                r_rat = data_dict["rp/re"][i]
                fname = runner.generate_filename(l1, l2, ec, r_rat)
                
                res = self.parse_output(fname)
                if res:
                    for key in mined:
                        mined[key].append(res.get(key, np.nan))
                    valid_indices.append(i)

            if not valid_indices:
                continue

            errors[k_pair] = {}
            
            # Calculate errors only for valid points
            for col in mined:
                mine_arr = np.array(mined[col])
                
                if col == "J/M^2":
                    ref_J = np.array(data_dict["J"])[valid_indices]
                    ref_M = np.array(data_dict["M"])[valid_indices]
                    ref_arr = ref_J / (ref_M**2)
                else:
                    ref_arr = np.array(data_dict[col])[valid_indices]
                
                # Relative Error
                err = np.abs(mine_arr - ref_arr) / np.where(ref_arr != 0, np.abs(ref_arr), 1.0)
                errors[k_pair][col] = err

        return errors

def plot_errors(errors, output_file="errors_comparison.pdf"):
    if not errors:
        logging.warning("No data to plot.")
        return

    cols = ["M", "Re", "J/M^2", "Ωc", "Ωe", "Ωmax"]
    labels = ["M", r"$R_e$", r"J/M$^2$", r"$\Omega_c$", r"$\Omega_e$", r"$\Omega_{max}$"]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Flatten data for boxplots
    data_to_plot = []
    for col in cols:
        col_errs = np.concatenate([v[col] for v in errors.values() if col in v])
        data_to_plot.append(col_errs)

    # Boxplot
    boxprops = dict(facecolor='none', edgecolor='black', linewidth=1.5)
    ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, boxprops=boxprops, showfliers=False)

    # Scatter (Jitter)
    colors = plt.cm.tab20(np.linspace(0, 1, len(errors)))
    for color, (k_pair, k_data) in zip(colors, errors.items()):
        for i, col in enumerate(cols):
            if col in k_data:
                y = k_data[col]
                x = np.random.normal(loc=i+1, scale=0.08, size=len(y))
                label = fr"$\lambda$={k_pair}" if i == 0 else ""
                ax.scatter(x, y, color=color, alpha=0.7, label=label, s=20)

    ax.set_ylabel("Relative Error", fontsize=14)
    ax.set_yscale('log')
    plt.xticks(rotation=45, fontsize=12)
    
    # Legend handling
    handles, lbls = ax.get_legend_handles_labels()
    by_label = dict(zip(lbls, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    logging.info(f"Plot saved to {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Run RNS regression tests.")
    parser.add_argument("--exec", required=True, help="Path to RNS executable")
    parser.add_argument("--data", default="./data.d", help="Path to reference data file")
    parser.add_argument("--configs", default="./configs", help="Directory for config files")
    parser.add_argument("--output", default="./outfiles", help="Directory for simulation output")
    parser.add_argument("--eos-bm", default="../../eos/DD2noY.rns", help="Path to BM EoS file")
    parser.add_argument("--eos-dm", default="../../eos/colpi_0.25GeV_4500_24PI.rns", help="Path to DM EoS file")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers")
    parser.add_argument("--plot", default="errors.pdf", help="Output filename for plot")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity")

    args = parser.parse_args()
    setup_logging(args.verbose)

    # 1. Load Data
    try:
        ref_data = ReferenceDataLoader.load(args.data)
    except Exception as e:
        logging.critical(f"Failed to load data: {e}")
        return

    # 2. Setup Runner
    eos_paths = {'bm': args.eos_bm, 'dm': args.eos_dm}
    runner = SimulationRunner(args.exec, args.configs, args.output, eos_paths)

    # 3. Prepare Tasks
    tasks = []
    for (l1, l2), data_dict in ref_data.items():
        for i in range(len(data_dict["εc"])):
            tasks.append((l1, l2, data_dict["εc"][i], data_dict["rp/re"][i]))

    # 4. Run Simulations
    runner.run_all(tasks, max_workers=args.workers)

    # 5. Analyze Results
    analyzer = ResultAnalyzer(args.output)
    errors = analyzer.compare(ref_data, runner)

    # 6. Report Stats
    logging.info("-" * 40)
    logging.info("ERROR STATISTICS")
    cols = ["M", "Re", "J/M^2", "Ωc", "Ωe", "Ωmax"]
    
    all_errors = {c: [] for c in cols}
    for k_data in errors.values():
        for c in cols:
            if c in k_data:
                all_errors[c].extend(k_data[c])
    
    for c in cols:
        arr = np.array(all_errors[c])
        if len(arr) > 0:
            logging.info(f"{c:10} | Median: {np.median(arr):.4e} | Max: {np.max(arr):.4e}")
    logging.info("-" * 40)

    # 7. Plot
    plot_errors(errors, args.plot)

if __name__ == "__main__":
    main()