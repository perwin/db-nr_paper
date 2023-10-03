# Public data, code, and notebooks for paper analyzing inner bars and nuclear rings

<!-- 
This git repository contains data files, Python code, and Python and R 
Jupyter notebooks which can be used to reproduce figures and analyses
from the paper "The Profiles of Bars in Barred Galaxies" (Erwin,
Debattista, & Anderson 2023, *Monthly Notices of the Royal Astronomical
Society*, in press).
 -->


![Trends](./composite_fig_for_github.png)

Left: fractions of S0-Sc galaxies with inner bars (black squares) or nuclear rings
(open blue circles), counting all
barred galaxies (top) or barred + unbarred galaxies (bottom; red circles = fraction
of galaxies with bars).
Right: Trends in bar semi-major axis versus galaxy stellar mass.

<!-- [![DOI](https://zenodo.org/badge/579425923.svg)](https://zenodo.org/badge/latestdoi/579425923) -->



## Data

The following data files are included:

   * `table_mainsample.dat` -- General galaxy and bar data for full (156-galaxy) sample
    of barred galaxies (the full version of Table 1 in the paper)
    
   * `table_innerbars.dat` -- Measurements of inner bars in double-barred galaxies (Table 2
   in the paper)

   * `table_nuclearrings.dat` -- Measurements of nuclear rings (Table 3 in the paper)


## Dependencies

The Python code and notebooks require the following external Python modules and packages,
all of which are available on PyPI and can be installed via `pip`:

   * [Numpy](https://www.numpy.org), [Scipy](https://www.scipy.org), 
   [matplotlib](https://matplotlib.org), [Astropy](https://www.astropy.org)


## Jupyter Notebooks

There are XXX Jupyter notebooks:

   * `XXX.ipynb` -- Python notebook; generates the figures for 
   the paper



## Python Code

<!-- 
   * `angle_utils.py`, `barprofile_utils.py`, `plotutils.py` -- miscellaneous utility functions
   (including statistics).
 -->
   


## How to Generate Figures and Analyses from the Paper

<!-- 
1. Download this repository.

2. Edit paths in the notebooks so they point to the correct locations, if necessary.
See notes in the initial cells of the notebooks; the main variable you will probably
need to edit is `plotDir` in the second cell of `barprofiles_figures_for_paper.ipynb`,
which is where saved PDF figures should go. Also make sure to set `savePlots = True`
if you want the PDF files to actually be generated (the default is `False`, which
means the figures will appear in the notebook but won't be saved to disk).

3. **Optionally:** Run the notebook
`barprofiles_R_logistic-regression.ipynb` to generate and save the
various logistic fits. This is "optional" in that the output files
already exist in this directory (they will be overwritten if the
notebook is run).

4. Run the notebook `barprofiles_figures_for_paper.ipynb` to generate the figures
(it will read the coefficients of the fits from the file generated by running the
previous notebook).
 -->


## Licensing

Code in this repository is released under the BSD 3-clause license.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" 
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />
Text and figures are licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
