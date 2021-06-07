The package provides utility code developed to facilitate the distribution and usage of IR templates presented in Boquien & Salim (2021). It contains three routines:

- bosa-templates for extracting custom templates (average IR spectra)
- bosa-estimators for estimating LTIR and SFR (total, obscured plus unobscured)
- bosa-plotter 

# Installation

The package requires python 3.8 or later and depends on three external libraries, numpy, matplotlib and astropy. The installation follows the standard python way: `python setup.py install`. Note that the data file take >140 MB, so if you want to keep the files locally rather than having a system-wide copy, you can do `python setup.py develop --user`.

# Content

The package contains:
- This `README.md` file.
- The source with the `estimators`, `plotter`, and `templates` packages.
- Data files in `pickle` format in `src/*/data`.

# Usage
## bosa-templates

This routine produces average IR spectra (templates) corresponding to arbitrary TIR luminosities (single parameter mode), or  arbitrary LTIR/sSFR combinations (two parameter mode). These luminosities must be provided as arguments and must be in logarithmic units and in solar luminosities between 10⁹ and 10¹² Lsun (a warning will be displayed if beyond these bounds). For instance `bosa-templates LTIR 10 11 12` will produce the `LTIR.fits` with three spectra. Note that the script will display an exception if a file already exists in order not to accidentally overwrite files. In single parameter mode LTIR can be replaced with some other parameter:  Mdust, Mstar, SFR, sSFR, or Z. Spectra averaged in LTIR, Mdust, Mstar, and SFR are in absolute luminosity, whereas the ones in sSFR and Z are normalized to 1 Lsun. Parameters are in log for all the properties except for Z, which is in 12+log O/H. In two parameter mode (e.g., `bosa-templates LTIR 10 11 12 sSFR -10 -9.5`) code produces `LTIR-sSFR.fits` with 6 spectra in this case.

The output files have `wavelength` in units of nm and `nuLnu` in units of solar luminosities.

## bosa-estimators

This routine estimates LTIR and the SFR based on one of more (up to four) fluxes as an input. One column must contain redshifts and be named `redshift` and one must be named `sSFR` and it should contain log sSFR if known, or some placeholder values if not. Fluxes are specified using column names:
- `spitzer.irac.ch4`
- `spitzer.mips.24`
- `herschel.pacs.70`
- `herschel.pacs.100`
- `herschel.pacs.160`
- `herschel.spire.PSW`
- `herschel.spire.PMW`
- `herschel.spire.PLW`
- `jwst.miri.F770W`
- `jwst.miri.F1000W`
- `jwst.miri.F1130W`
- `jwst.miri.F1280W`
- `jwst.miri.F1500W`
- `jwst.miri.F1800W`
- `jwst.miri.F2100W`
- `jwst.miri.F2550W`
- `WISE3`
- `WISE4`

Output includes `LTIR` in Lsun and `SFR` in Msun/year calculated based on each filter individually and L(TIR) parameterized templates, also on L(TIR)+sSFR parameterized templates, and finally based on all flux points simultaneously. In order not to accidentally erase data, the script will refuse to run if columns of that name are already present in the file.

# Contact

Questions and suggestions are welcome and can be sent to Médéric Boquien, mederic.boquien@uantof.cl.
