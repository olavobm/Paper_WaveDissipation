# Paper_WaveDissipation

The repository Paper_WaveDissipation has code to analyze data and generate figures for the science paper "Observations of Wave Energy Dissipation by Bottom Friction on Rocky Shores" [(Marques et al. 2025)](https://journals.ametsoc.org/view/journals/phoc/aop/JPO-D-24-0144.1/JPO-D-24-0144.1.xml), accepted for publication in the Journal of Physical Oceanography (JPO).

In this paper, we analyzed measurements collected on rocky shores of the Monterey Peninsula as part of the [ROXSI (ROcky shores: eXperiments and SImulations)](https://roxsi.ucsd.edu/) 2022 experiment. In particular, we observed significant dissipation of sea-swell waves outside of the surfzone, which indicates large wave dissipation by bottom friction. The dissipation observed between pairs of instruments provides estimates of the wave friction factor $f_e$. The observed $f_e$ typically varies between 1 and 10, and are amongst the largest friction factors observed in environments with rough seabeds, including coral reefs and rocky shores.

At the top directory level, this GitHub repository has:
* paper_directory.m: a function that returns the paper directory path.
* run_reproducepaper.m: a high-level function to run all scripts that process data and make figures.
* proc/: a directory with Matlab code to process data.
* figures/: a directory with Matlab code to generate the figures.

<!---
The data has been archived in this [Zenodo repository](https://doi.org/.../zenodo....).
-->

To reproduce the paper, save the data repository (i.e., a folder named data) in the same level as the content in the Paper_WaveDissipation GitHub repository.
To fully reproduce the figures, you must also have two additional toolboxes:
* [m_map](https://www-old.eoas.ubc.ca/~rich/map.html) to make a map in Fig. 1.
* [cmocean](https://github.com/chadagreene/cmocean) for a few colormaps.

For using the repository, **first edit the directory path in paper_directory.m** so that it provides the paper repository path in your machine.
After editing paper_directory.m, you may run the full analysis with the high-level script run_reproducepaper.m.

