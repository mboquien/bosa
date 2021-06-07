import pickle
from pkg_resources import resource_stream as rs
import sys

from astropy.cosmology import WMAP7
from astropy.table import Column, Table
from astropy.io.registry import IORegistryError
import numpy as np
from scipy.constants import c, parsec
import scipy.interpolate as interp


class Estimator:
    def __init__(self, param):
        """Estimator of a physical property p for the luminosity of one to four
        bands. It is based on a simple relation:
            p = Σα×log(luminosity)+β.
        The coefficients α and β have been calibrated first by fitting a large
        number of galaxies with WISE+Herschel bands to estimate their physical
        properties and then relating the luminosities to these properties
        through a multi-linear relation in log space. These α and β
        coefficients are computed for all possible combination of up to four
        WISE, JWST, and Herschel bands. They are loaded through pickle files
        that are provided with the package.
        The value of p is in log for all the physical properties except for Z.

        Parameters
        ----------
        param: str
            Name of the parameter to be estimated to load the appropriate
            pickle files. Normally it should be SFR or LTIR
        """
        self.alpha = pickle.load(rs('estimators', f'data/alpha-{param}.pickle'))
        self.beta = pickle.load(rs('estimators', f'data/beta-{param}.pickle'))

        self.alpha2nd = pickle.load(
            rs('estimators', f'data/alpha-2nd-{param}.pickle')
        )
        self.beta2nd = pickle.load(
            rs('estimators', f'data/beta-2nd-{param}.pickle')
        )

        # Set up interpolators so the estimators work at any redshift
        self.z = np.linspace(0.0, 4.0, 401)
        self.alpha = {
            band: interp.interp1d(self.z, self.alpha[band], axis=0)
            for band in self.alpha
        }
        self.beta = {
            band: interp.interp1d(self.z, self.beta[band], axis=0)
            for band in self.beta
        }
        self.alpha2nd = {
            band: interp.interp1d(self.z, self.alpha2nd[band], axis=0)
            for band in self.alpha2nd
        }
        self.beta2nd = {
            band: interp.interp1d(self.z, self.beta2nd[band], axis=0)
            for band in self.beta2nd
        }

        # Build the line of bands
        self.bands = tuple(
            band[0] for band in self.alpha.keys() if len(band) == 1
        )

    def estimate(self, observations, bands=None):
        """Estimate the physical properties based on certain observations and
        based a given set of bands. Estimates are given both when using all the
        input bands simultaneously (up to four), and using each of these bands
        individually.

        Parameters
        ----------
        observations: astropy.Table
            Table containing the redshift and the luminosities in each band for
            each object.
        bands: list
            List of bands to use for computing the estimates

        Returns
        -------
        estimates: dictionary
            Contains all the estimates of the relevant physical property. The
            key is of the form 'property[band]', with property being SFR or
            LTIR and 'band' being the band name. When all the bands are used,
            the key is 'allbands'.
        """
        # If bands are not indicated explicitly, we take all the bands from the
        # input file
        if bands is None:
            bands = [
                band for band in observations.colnames if band in self.bands
            ]

        if len(bands) < 1 or len(bands) > 4:
            raise Exception(
                "Between one and four bands must be provided. The input bands "
                f"are: {', '.join(bands)}."
            )

        # Find the corresponding set of bands. This is necessary because the
        # bands mustbe sorted following a certain order
        for key in self.alpha.keys():
            if set(key) == set(bands):
                break

        # The key should always be found, but better double-check
        if set(key) != set(bands):
            raise Exception(
                f"Combination of bands {', '.join(bands)} not found."
            )

        # Verify that the redshifts as all within the bounds
        z = observations['redshift']
        if np.any((z < 0.0) | (z > 4.0)):
            raise Exception("The redshifts must all be between 0 and 4.")

        estimates = {}
        estimates2nd = {}

        # Compute the estimate using all the bands, if it makes sense
        if len(bands) > 1:
            estimate = self.beta[key](z)
            for idx_band, alpha in enumerate(self.alpha[key](z).T):
                estimate += np.log10(observations[key[idx_band]].data) * alpha
                wul = np.where(observations[f'{key[idx_band]}_err'].data <= 0.0)
                if len(wul[0]) > 0:
                    estimate[wul] = np.nan
            estimates['allbands'] = 10 ** estimate

        # Compute the estimates using individual bands
        for band in bands:
            lumin = np.log10(observations[band].data)
            estimates[band] = 10 ** (
                lumin * np.squeeze(self.alpha[(band,)](z))
                + self.beta[(band,)](z)
            )
            alpha = np.squeeze(self.alpha2nd[(band,)](z))

            estimates2nd[band] = 10 ** (
                lumin * alpha[:, 0]
                + observations['sSFR'] * alpha[:, 1]
                + self.beta2nd[(band,)](z)
            )

        return estimates, estimates2nd

    @staticmethod
    def prepare_data(data):
        """Convert the input flux densities normally given in mJy into solar
        luminosities given in W. We assume a WMAP7 cosmology, as for the
        derivation of the estimators.

        Parameters
        ----------
        data: astropy.Table
           Contains the redshift and the flux density of each object in mJy

        Returns
        -------
        data: astropy.Table
           Contains the redshift and the luminosity of each object in solar
           luminosities
        """
        # Pivot wavelength in m of the filters used to compute the estimators
        wl = {
            'spitzer.irac.ch4': 7.915762440837297e-06,
            'spitzer.mips.24': 2.3593476810206795e-05,
            'herschel.pacs.70': 7.073074183063004e-05,
            'herschel.pacs.100': 0.00010065133201547037,
            'herschel.pacs.160': 0.00016097632066582225,
            'herschel.spire.PSW': 0.000248540338729174,
            'herschel.spire.PMW': 0.0003481595910721862,
            'herschel.spire.PLW': 0.0005007949402133246,
            'jwst.miri.F770W': 7.639333710593687e-06,
            'jwst.miri.F1000W': 9.953116516778535e-06,
            'jwst.miri.F1130W': 1.1308501361802308e-05,
            'jwst.miri.F1280W': 1.2810137914219364e-05,
            'jwst.miri.F1500W': 1.5063506136295945e-05,
            'jwst.miri.F1800W': 1.7983722293636374e-05,
            'jwst.miri.F2100W': 2.079500506251338e-05,
            'jwst.miri.F2550W': 2.5364001534544257e-05,
            'WISE3': 1.2069436891637823e-05,
            'WISE4': 2.219568031582113e-05
        }

        if 'distance' not in data.colnames:
            dist = WMAP7.luminosity_distance(data['redshift']).to('m').value
            w0 = np.where(data['redshift'] == 0.0)
            if len(w0[0]) > 0:
                data[w0] = 10.0 * parsec
        else:
            dist = data.distance.value * 1e6 * parsec

        cst = 1e-29 * c * 4.0 * np.pi / 3.828e26 * dist ** 2
        conv = {k: cst / wl[k] for k in wl}
        for band in wl:
            if band in data.colnames:
                data[band] *= conv[band]
                if (banderr := f'{band}_err') in data.colnames:
                    data[banderr] *= conv[band]
                else:
                    data[banderr] = data[band] * 0.1

        return data


def main():
    # Instantiate the SFR and LTIR estimators
    estimatorSFR = Estimator('SFR')
    estimatorLTIR = Estimator('LTIR')

    # Check that we have at least one input parameter
    if len(sys.argv) < 2:
        bands = '\n'.join(estimatorSFR.bands)
        raise Exception(
            "You must provide the name of the file. This file must contain a "
            "column 'redshift' with the redshift and columns for the flux in "
            f"mJy. The band names are: {bands}.\nNote that the names are case "
            "sensitive. Optionally you can also provide the distance in Mpc, "
            "which will take precedence over the redshift for computing the "
            "distance. If no  distance is given and the redshift is 0, a "
            "distance of 10 pc is assumed."
        )

    # If we have more than one parameter, we take the additional ones as bands
    fname = sys.argv[1]
    if len(sys.argv) > 2:
        bands = sys.argv[2:]
    else:
        bands = None

    # Attempt to read the input files either as ASCII or as FITS
    try:
        data = Table.read(fname)
    except IORegistryError:
        try:
            data = Table.read(fname, format='ascii', delimiter=r'\s')
        except IORegistryError:
            raise Exception(f"The file {fname} could not be parsed.")

    # We propare the data. This mostly consists in converting observed flux
    # densities in mJy to luminosities in solar luminosities
    Ldata = Estimator.prepare_data(data)

    # Estimate the SFR and LTIR
    SFR, SFR2nd = estimatorSFR.estimate(Ldata, bands)
    LTIR, LTIR2nd = estimatorLTIR.estimate(Ldata, bands)

    estimates = Ldata.copy()
    for k in SFR:
        estimates.add_column(Column(SFR[k], name=f'SFR[{k}]', unit='Msun/year'))
        estimates.add_column(Column(LTIR[k], name=f'LTIR[{k}]', unit='W'))
    for k in SFR2nd:
        estimates.add_column(
            Column(SFR2nd[k], name=f'SFR2nd[{k}]', unit='Msun/year')
        )
        estimates.add_column(Column(LTIR2nd[k], name=f'LTIR2nd[{k}]', unit='W'))

    try:
        estimates.write(f'output-{fname}')
    except IORegistryError:
        estimates.write(
            f'output-{fname}', delimiter='', format='ascii.fixed_width'
        )
