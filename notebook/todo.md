Repo: https://github.com/BenjMy/tmpCATHY_ERT_DA_real/tree/feature/data-paths-registry-readme

- [x] Add a manifest file for each run to check which version of CATHY DA ERT (including the branch) is used and which version of CATHY and pyCATHY as well (including the branch). And other usefuls things as you consider for reproducibility.
- [ ] feat/withIrr-on-main or feature/eta-assimilation-improvements or twin-experiment? Help me to choice the last updated DA ERT branch.  Should I rebase?
- [ ] By the way can I close origin/feat/material-vgp-control? and origin/mesh-vtk-properties? 
- [ ] I want a data_path.toml where I can control all the path for inputs, outputs, plots, codes lib etc ...
- [ ] Register all the run in a csv files where I can keep track what has been done or not
- [ ] Investigate if it is possible to get ligher DA folder (.pkl files are heavy)
- [ ] Create a new branch on the github and add a README that explains all the code architecture, from inputs type, how to run, scenario construction, outputs data and plots. 
- [ ] Explain synthetic ERT DA pipeline in the README
- [ ] Explain how to use scenario.toml file in the README
- [ ] Check this warming:   warnings.warn("truth forward on the SAME mesh as the DA: mesh discretisation error is not tested "
- [ ] ERT frame: legacy (farming_DA_ERT3d.py); signs=(1, 1, 1), shift=[2.0, 2.0, 0.001]
/home/ben/Nextcloud/BenCSIC/Codes/BenjMy/pycathy_wrapper/pyCATHY/ERT/alignment.py:116: UserWarning: [ERT] pg forward mesh extends only 0.50 m beyond the electrodes (array span 1.10 m): boundary effects can bias the forward apparent resistivities of weak-signal quadrupoles (see nc_model.unreliable_quadrupoles: homogeneous test, rhoa/rho must be 1)
  warnings.warn(f"[{label}] pg forward mesh extends only {margin:.2f} m beyond the 

use .patch file that you commit as me as author (BenjMy)



