# MoS2 Phase Diagram

Python script for the Supporting Information of A.S. Rosen, J.M. Notestein, R.Q. Snurr, "Comprehensive Phase Diagrams of MoS2 Edge Sites using Dispersion-Corrected DFT Free Energy Calculations" in J. Phys. Chem. C (2018). https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.8b02524

## Instructions
To run mos2_data.py, open the file and change T (temperature in Kelvins), f_HS2 (fugacity of H2S in bar), and f_H2 (fugacity of H2 in bar) to the desired parameters. Partial pressures can be used in place of fugacities if ideal gas conditions can be assumed. Make sure the excel_path variable points to the included Excel sheet. Then, run the Python script (e.g. via `python mos2_data.py` from the command line). 

## Output
It will print out the most stable Mo-edge and S-edge configuration at the specified conditions. It will also save four .csv files in a folder: `F_Mo_edge.csv` (Helmholtz free energy of the Mo-edge), `F_S_edge.csv` (Helmholtz free energy of the S-edge), `phi_Mo_edge.csv` (grand potential of the Mo-edge), `phi_S_edge.csv` (grand potential of the S-edge). All energies are in eV. The .csv files have 7 rows and 7 columns, representing theta_S = {0, 0.17, 0.33, 0.5, 0.67, 0.83, 1} and theta_H = {0, 0.33, 0.66, 1.0, 1.33, 1.66, 2.0}, respectively. Only theta_S = 1 has entries for theta_H = 1.33 - 2.0.
