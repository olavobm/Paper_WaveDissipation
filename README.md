# Paper_WaveDissipation
Code for reproducing paper on sea-swell wave dissipation over rocky shores.

The repositoy Paper_WaveDissipation has code to analyze data and generate figures for the scientific paper "Observations of Wave Energy Dissipation by Bottom Friction on Rocky Shores", accepted for publication in the Journal of Physical Oceanography (JPO).
<!---
[Marques et al. 2025](https://journals.ametsoc.org/view/journals/atot/....).
-->

At the top directory level, this GitHub repository has:
* paper_directory.m: a function that returns the paper directory path.
* run_reproducepaper.m: a high-level function to run all scripts that process data and make figures.
* proc/: a directory with Matlab code to process data.
* figures/: a directory with Matlab code to generate the figures.


The data has been archived in this [Zenodo repository](https://doi.org/10.5281/zenodo.13242438). To reproduce the analysis, save the data repository (a folder named data) in the same level as the content in the Paper_EffectiveDepth GitHub repository.
To reproduce the figures, you must also have [cmocean](https://github.com/chadagreene/cmocean) for a few colormaps.

For using the repository, **first edit the directory path in paper_directory.m** to reflect the appropriate paper repository path in your machine.
After editing paper_directory.m, you may run the full analysis with run_all.m. This script simply runs the other high-level scripts ~/data_proc/run_alldataprocessing.m and ~/figures/run_allfigures.m.

