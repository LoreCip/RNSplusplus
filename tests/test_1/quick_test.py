#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import h5py
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# --- Physical Constants & Conversions ---
C_CGS = 2.99792458e10
G_CGS = 6.6730831e-8
SQRT_KAPPA = np.sqrt(1.0e-15 * C_CGS * C_CGS / G_CGS)
CU_TO_KHZ = (C_CGS / SQRT_KAPPA) / 1000.0

# --- Hardcoded Reference Data ({λ1, λ2} = {2.0, 0.5}) ---
REF = {
    "lambda1": 2.0,
    "lambda2": 0.5,
    "r_ratio": 0.3208,
    "ec":      0.3812,
    # Expected Results from inspirehep.net/literature/1861167
    "M":       2.20,
    "M0":      2.4009,
    "J":       4.4513,
    "Oc":      7.6442,
    "Omax":    15.2884,
    "Oe":      3.8221,
    "Re":      18.1234
}

CONFIG_TEMPLATE = """# Quick Test Config
EoS_BM_type         0
EoS_BM_file         ../../eos/DD2noY.rns
EoS_DM_type         0
EoS_DM_file         ../../eos/colpi_0.25GeV_4500_24PI.rns
accuracy            1e-10

BM_central_energy   {ec}
DM_fraction         0.0
BM_rotation_type    2
BM_A                0.7
BM_B                0.5
BM_lambda1          {l1}
BM_lambda2          {l2}
r_ratio_BM          {r_ratio}

2Doutput            1
output_name         {name}
verbose             1
"""

def run_model(exe):
    # Construct base filename
    name = f"quick_test_l1{REF['lambda1']}_l2{REF['lambda2']}_e{REF['ec']}_r{REF['r_ratio']}"
    config_file = name + ".d"
    out_file = name + ".out"
    h5_local_name = name + ".h5"

    # 1. Generate Configuration File
    print(f"--- Generating config: {config_file} ---")
    with open(config_file, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            ec=REF['ec'], l1=REF['lambda1'], l2=REF['lambda2'], 
            r_ratio=REF['r_ratio'], name=name
        ))

    # 2. Run Simulation
    print(f"--- Running {exe} ---")
    result = subprocess.run([exe, "-c", config_file], capture_output=True, text=True)
    
    # Save standard output for parsing
    with open(out_file, "w") as f:
        f.write(result.stdout)

    # 3. Handle File Management (Move HDF5 to current folder)
    # The code saves it in ./output/output_2d/ based on internal structure
    h5_remote_path = os.path.join("output", "output_2d", h5_local_name)
    
    if os.path.exists(h5_remote_path):
        print(f"--- Moving HDF5 from {h5_remote_path} to current directory ---")
        shutil.move(h5_remote_path, h5_local_name)
    else:
        print(f"--- Warning: Could not find HDF5 file at {h5_remote_path} ---")

    # 4. Parse Output Data
    data_mine = {}
    with open(out_file, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Locate Baryonic Matter data block (ignoring target setting lines)
            if "Baryonic Matter:" in line and "Target" not in line:
                # Based on RNS output format, the data line is usually 2 lines after the header
                if i + 2 < len(lines):
                    parts = lines[i+2].split()
                    try:
                        data_mine["M"] = float(parts[0])
                        data_mine["M0"] = float(parts[1])
                        data_mine["Re"] = float(parts[2])
                        # The code outputs J/M^2; we convert to J for direct comparison
                        data_mine["J"] = float(parts[3]) * (data_mine["M"]**2)
                    except (ValueError, IndexError):
                        continue
                break

    # Parse HDF5 for Frequency Profiles
    if os.path.exists(h5_local_name):
        try:
            with h5py.File(h5_local_name, "r") as hf:
                # Omega_diffBM shape is typically (SDIV, MDIV)
                O_diff = hf["Omega_diffBM"][:, 1] * CU_TO_KHZ
                s_qp = hf["s_qp"][:]
                
                # Interpolate to find specific points
                interp = interp1d(s_qp, O_diff, kind='linear', fill_value="extrapolate")
                data_mine["Oc"] = float(interp(0))
                data_mine["Oe"] = float(interp(0.5))
                
                # Find the maximum Omega
                res = minimize_scalar(lambda s: -interp(s), bounds=(0, 0.5), method='bounded')
                data_mine["Omax"] = float(interp(res.x))
        except Exception as e:
            print(f"Error parsing HDF5 data: {e}")
    else:
        print("Error: HDF5 file missing. Frequency comparison will be skipped.")

    # 5. Compare Results with Reference
    print(f"\n{'Quantity':<10} | {'Reference':<12} | {'Simulation':<12} | {'Rel. Error':<12}")
    print("-" * 60)
    
    keys_to_compare = [
        ("M", "Mass (M)"), ("M0", "Mass0 (M0)"), ("Re", "Radius (Re)"), ("J", "Ang. Mom (J)"), 
        ("Oc", "Ω_center"), ("Omax", "Ω_max"), ("Oe", "Ω_equator")
    ]

    for ref_key, label in keys_to_compare:
        val_ref = REF[ref_key]
        val_mine = data_mine.get(ref_key, np.nan)
        
        if not np.isnan(val_mine):
            err = abs(val_mine - val_ref) / abs(val_ref) if val_ref != 0 else 0
            print(f"{label:<10} | {val_ref:<12.4f} | {val_mine:<12.4f} | {err:.4e}")
        else:
            print(f"{label:<10} | {val_ref:<12.4f} | {'FAILED':<12} | {'N/A':<12}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quick single-model test for RNSDiff_one.")
    parser.add_argument("--exec", default="../../RNSDiff_one", help="Path to the C executable")
    args = parser.parse_args()
    
    if not os.path.exists(args.exec):
        print(f"Error: Executable not found at {args.exec}")
    else:
        run_model(args.exec)