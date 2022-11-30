# Formation of the Helix Nebula

Data and analysis scripts for publication "Evidence for the disruption of a planetary system during the formation of the Helix Nebula", Marshall, J.P. et al., 2022, AJ accepted.

AJ:

arXiv: https://arxiv.org/abs/2211.02251

NASA ADS: https://ui.adsabs.harvard.edu/abs/2022arXiv221102251M/abstract

## Contents

### Data
./data

./data/imaging/ -- contains imaging data used in this work (ALMA, Herschel, SOFIA, Spitzer). These are arranged by program within each sub-directory. Imaging data are stored as fits files. The original files can also be obtained from the observatories' respective public archives. 

./data/photometry -- contains photometric data, as presented in Table 1 of the paper. There are three text files, one for each component of the SED (compact, extended, total), and a csv file which is used in the modelling.

./data/spectroscopy -- contains spectroscopic data used in this work (ALMA, Spitzer). These are arranged by program within each sub-directory. Spectroscopic data are stored as text files. The Spitzer IRS spectrum used in this work can also be obtained from the CASSIS archive.

### Models
./models

./models/dust/ -- contains optical constants for astronomical silicate used in the radiative transfer modelling presented in Sections 3.2, 4.2, and Figure 4 of the paper.

./models/photosphere/ -- contains the white dwarf model photosphere used in this work, taken from the Spanish Virtual Observatory.

### Scripts
./scripts -- contains all the scripts used in data analysis for this work. The scripts have been named after which figure in the paper they will produce when run.

## Instructions

Scripts can be run in their holding directory, output will be placed in the directory above that (i.e. '../'). Currently the Herschel two component image analysis script (Figure 2) is written in IDL, whilst all the remaining scripts are written in Python. It is planned to provide Python scripts for all aspects of the analysis. Radiative transfer analysis (Figure 4) was done using 'artefact', which is available as a separate repository. 

