# Public data, code, and notebooks for paper analyzing inner bars and nuclear rings

This git repository contains data files, Python code, and Jupyter
notebooks (for Python and R) which can be used to reproduce figures and
analyses from the paper "The Frequency and Sizes of Inner Bars and
Nuclear Rings in Barred Galaxies and Their Dependence on Galaxy
Properties" (Erwin, 2023, *Monthly Notices of the Royal Astronomical
Society*, in press; [arXiv:2312.12893](http://arxiv.org/abs/2312.12893)).


![Trends](./composite_fig_for_github.png)

Left: fractions of S0-Sd galaxies with inner bars (black squares) or nuclear rings
(open blue circles), counting all barred galaxies (top) or barred + unbarred galaxies 
(bottom; red circles = fraction of galaxies with bars).
Right: Trends in bar semi-major axis versus galaxy stellar mass.

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10417791.svg)](https://doi.org/10.5281/zenodo.10417791) -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10252782.svg)](https://doi.org/10.5281/zenodo.10252782)



## Data

The following data files are included:

   * `table_mainsample.dat` -- General galaxy and bar data for full (155-galaxy) sample
    of barred galaxies (the full version of Table 1 in the paper)
    
   * `table_innerbars.dat` -- Measurements of inner bars in double-barred galaxies 
   (the full version of Table 2 in the paper)

   * `table_nuclearrings.dat` -- Measurements of nuclear rings (the full version of Table 3 
   in the paper)

   * `table_unbarred_info.dat` -- Minimal data for unbarred galaxies in parent sample

   * `s4gbars_table.dat` -- data for 1322 individual galaxies in the Parent Disc Sample
   of Erwin (2018; used for Figs. 9 and 10)


## Dependencies

The Python code and notebooks require the following external Python modules and packages,
all of which are available on PyPI and can be installed via `pip`:

   * [Numpy](https://www.numpy.org), [Scipy](https://www.scipy.org), 
   [matplotlib](https://matplotlib.org)


## Jupyter Notebooks

There are three Jupyter notebooks:

   * `figures_for_paper.ipynb` -- Python notebook; generates the figures for the paper

   * `fits_for_paper.ipynb` -- Python notebook; computes best-fit parameters, AIC
   values, uncertainties for various fits

   * `logistic_fits_r.ipynb` -- R notebook; performs logistic fits for DB and NR presence
   as functions of stellar mass and/or bar size


## Python Code

   * `datautils.py`, `dbnr_utils.py`, `fitting_barsizes.py`, `plotutils.py` -- miscellaneous 
   utility functions for reading data tables, performing fits, and generating figures.
   


## How to Generate Figures and Analyses from the Paper

1. Download this repository.

2. Edit paths in the notebooks so they point to the correct locations, if necessary.
See notes in the initial cells of the notebooks; the main variable you will probably
need to edit is `plotDir` in the second cell of `figures_for_paper.ipynb`,
which is where saved PDF figures should go. Also make sure to set `savePlots = True`
if you want the PDF files to actually be generated (the default is `False`, which
means the figures will appear in the notebook but won't be saved to disk).

3. Run the notebook `figures_for_paper.ipynb` to generate the figures.

4. Run the notebook `fits_for_paper.ipynb` if you wish to do the fits for inner-bar
and nuclear-ring sizes as functions of galaxy stellar mass and/or outer/main bar size.

5. Run the notebook `logistic_fits_r.ipynb` if you wish to do the logistic regression
for inner-bar and nuclear-ring presence as functions of galaxy stellar mass and/or 
outer/main bar size.



## Licensing

Code in this repository is released under the BSD 3-clause license.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" 
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />
Text and figures are licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
