# tDCSsim_Aberra2023

Code associated with Aberra AS, Wang R, Grill WM, Peterchev AV. (2022). "Multi-scale model of neuronal polarization by conventional and high-definition transcranial direct current stimulation in realistic head geometry". (in preparation)

This repository contains MATLAB code for reproducing figures and the neuron repositioning algorithm used for placing neurons within their layers while minimizing morphological intersections with the GM/CSF boundary. 

To be able to run all relevant functions:

1) Run `addPaths_tDCSsim.m`

2) To plot figures using data simulated for the manuscript, first download the data from the following Zenodo repository:

[https://doi.org/10.5281/zenodo.10275614](https://doi.org/10.5281/zenodo.10275614)

3) Place contents of `tDCSsim_Aberra2023_data` in `tDCSsim_Aberra2023` (top level directory):
    `nrn_sim_data/`
    `cell_data/` 
    `output_data`
    `simnibs_data` 
    OR, optionally (to save save space on your main HD), `nrn_sim_data` can be placed in a different folder (e.g. on an external HD). In this case, the path to the parent directory should be saved in a text file named `data_dir.txt` placed in the top level. For example, if the data is saved in `/Volumes/MyExternalHD/tDCSsim_Aberra2023_data/nrn_sim_data`, `tDCSsim_Aberra2023/data_dir.txt` should contain the following text:

>`/Volumes/MyExternalHD/tDCSsim_Aberra2023_data`

4) You can now run the functions for generating figures, found in the `gen_figures/` folder

5) Alternatively, you can run the neuron repositioning algorithm on all positions for any cell type in the population using `repositionNeuronLayer.m`. Default arguments are for L2/3 PC clone #1. Note: Due to computational cost parallelization is recommended. A parallel for loop can be used to take advantage of all available CPUs by setting the optional argument `par_on` to 1. 

## Dependencies:
MATLAB R2022a
Toolboxes:
Statistics and Machine Learning Toolbox
Parallel Computing Toolbox (Optional, required for parallelized for loops)
