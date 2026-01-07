# RNS Regression Test Suite

This code provides automated testing for the **RNSDiff_one** numerical relativity code. It is designed to ensure the code's accuracy by comparing its results against reference data for rotating dark matter-admixed neutron stars.

---

The primary objective of this test is to verify that the code correctly solves the field equations for rotating stars in General Relativity. The reference data used for comparison is taken from the publication: 

> **"Models of binary neutron star remnants with tabulated equations of state"** > [inspirehep.net/literature/1861167](http://inspirehep.net/literature/1861167)

### What is being tested?
The script automates the comparison of several fundamental physical quantities:
* **Global Properties**: Total gravitational mass $M$, baryonic mass $M_0$, and the equatorial radius $R_e$.
* **Rotation Dynamics**: The angular momentum $J$ (normalized as $J/M^2$) and the angular velocity profile $\Omega$.
* **Differential Rotation**: Specifically, values for $\Omega$ at the center ($\Omega_c$), the equator ($\Omega_e$), and its maximum value ($\Omega_{max}$) based on the Uryu rotation law.

For each model, the script calculates the **relative error**:
$$\text{Error} = \left| \frac{\text{Simulated Value} - \text{Reference Value}}{\text{Reference Value}} \right|$$
This allows us to monitor the numerical precision across different physical regimes, such as varying central energy densities $\varepsilon_c$ and axis ratios $r_p/r_e$.

---

### Usage

Run the comparison script by pointing it to your executable:

```bash
python3 full_comparison.py --exec ../RNSDiff_one --data data.d
```

#### Special Instructions for $\lambda_1 = 1.6$ Models

Some models (specifically where $\lambda_1 = 1.6$) are numerically sensitive. If you encounter convergence issues:

1.  **Tinker with `r_step`**: You may need to manually adjust the `r_step` parameter in the configuration file until the model converges.
2.  **Use Initial Guesses**: Once a single model works, use its HDF5 output as a starting guess for the remaining models in that sequence.
3.  **Config Field**: To use an initial guess, add the following field to your config file:
    ```text
    id_file    path/to/your/initial_guess.h5
    ```

---

### Workflow

- **Initialization**: The script reads data.d and generates unique configuration files in configs/.
- **Execution**: It runs the C code in parallel (using ProcessPoolExecutor) and tracks progress via a status bar.
- **Parsing**: It extracts global data from .out logs and detailed profiles from .h5 files.
- **Analysis**: It calculates error statistics (Median/Max) for all physical quantities.
- **Visualization**: It generates a boxplot (errors.pdf) to visualize the distribution of errors.

### Output: 

- **Boxplot**: Shows the statistical spread of errors.
- **Jittered Scatter**: Each point representsa single simulation, color-coded by its $\lambda$ (rotation parameter) set.
- **Outliers**: Any point significantly higher than the rest indicates a potential issue in that specific physical configuration.
