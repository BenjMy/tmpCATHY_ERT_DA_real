Repo: https://github.com/BenjMy/tmpCATHY_ERT_DA_real/tree/feature/data-paths-registry-readme

- [x] Add a manifest file for each run to check which version of CATHY DA ERT (including the branch) is used and which version of CATHY and pyCATHY as well (including the branch). And other usefuls things as you consider for reproducibility.
- [x] feat/withIrr-on-main or feature/eta-assimilation-improvements or twin-experiment? Help me to choice the last updated DA ERT branch.  Should I rebase?
  (checked 2026-11-24 on BenjMy/pycathy_wrapper; the branches live there, not in this repo)
  - DA branch to use: `feature/eta-assimilation-improvements` (4 commits not in main: set_atmbc_components / rain-ETp
    perturbation / mass-balance ETa operator, 2 mbeconv importer fixes, plot_state_timeseries). It is 7 commits behind main
    (merge-base 2026-10-26, before the mesh-vtk-properties work).
  - Yes, rebase it on main: `git rebase origin/main` is clean (no conflict, py_compile OK), then fast-forward/PR into main.
  - `feat/withIrr-on-main` is the most recently updated (2026-11-16) but it is the irrigation work (irrigation.py,
    cathy_tools, version_config), not DA; it is 0 behind main (main already merged into it). It touches the same files
    as eta (cathy_tools.py, cathy_plots.py) but merges with eta without conflict (checked): merge them in separate PRs.
  - `withIrrigation` is the older version of withIrr-on-main (68 behind): delete after withIrr-on-main is merged.
  - `twin-experiment` does not exist on origin of pycathy_wrapper, tmpCATHY_ERT_DA_real or mfc2_EA-EM_TDR: local-only?
    Push it or tell me where it is; not assessed.
- [x] By the way can I close origin/feat/material-vgp-control? and origin/mesh-vtk-properties? 
  - Yes, both. Each is an ancestor of origin/main (0 commits ahead; `git merge-base --is-ancestor` = YES;
    mesh-vtk-properties was merged by PR #33). Nothing is lost by deleting them:
    `git push origin --delete feat/material-vgp-control mesh-vtk-properties` (and `git fetch --prune` elsewhere).
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



