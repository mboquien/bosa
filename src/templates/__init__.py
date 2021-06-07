from itertools import product
import os.path
import pickle
import sys

from astropy.table import Column, Table
import numpy as np
from pkg_resources import resource_stream as rs


class Spectra:
    def __init__(self):
        """Constructs dust emission spectra parametrized on the
        logarithms of any of the following physical properties: LTIR (Lsun),
        Mdust (Msun), Mstar (Msun), SFR (Msun yr¯¹), sSFR (yr¯¹), and Z.
        The construction is based on the relations between the physical
        property p and the luminosity λLλ (Lsun) at wavelength λ:
        log(λLλ) = α×p+β.  The coefficients α and β have been calibrated
        from a large number of galaxies with WISE+Herschel bands.
        Spectra averaged in LTIR, Mdust, Mstar, and SFR are in absolute
        luminosity, whereas the ones in sSFR and Z are normalized to 1
        Lsun.  p is in log for all the properties except for Z, which is
        in 12+log O/H.
        """
        self.alpha = pickle.load(rs('templates', 'data/alpha-spec.pickle'))
        self.beta = pickle.load(rs('templates', 'data/beta-spec.pickle'))
        self.wl = pickle.load(rs('templates', 'data/wl-spec.pickle'))

    def estimate(self, param, value):
        """Constructs the IR spectrum corresponding to a certain value of a
        physical property p with a relation:
            log(λLλ) = α×log(p)+β.

        Parameters
        ----------
        param: str
            Physical parameter of the templates. It must be LTIR, Mdust, Mstar,
            SFR, sSFR, or Z.
        value: float
            Value of the physical property.

        Returns
        -------
        estimate: array
            Spectrum corresponding to a given value of param
        """
        return self.alpha[param] * value + self.beta[param]


class Spectra2nd:
    def __init__(self):
        """Constructs dust emission spectra parametrized on two physical
        properties. Supported at this moment are LTIR (Lsun) in
        combination with sSFR (yr¯¹).  The construction is based on a
        relation between physical properties p1 and p2 and the
        luminosity λLλ (Lsun) at wavelength λ: log(λLλ) = α×p1+β×p2+γ.
        The coefficients α, β, and γ have been calibrated by fitting a
        large number of galaxies with WISE+Herschel bands.
        """

        self.alpha = pickle.load(rs('templates', 'data/alpha-spec-2nd.pickle'))
        self.beta = pickle.load(rs('templates', 'data/beta-spec-2nd.pickle'))
        self.gamma = pickle.load(rs('templates', 'data/gamma-spec-2nd.pickle'))
        self.wl = pickle.load(rs('templates', 'data/wl-spec-2nd.pickle'))

    def estimate(self, param, value1, value2):
        """Constructs the spectrum corresponding to a certain value of a
        physical property p with a relation:
            log(λLλ) = α×log(p1)+β×log(p2)+γ.

        Parameters
        ----------
        param: str
            Physical parameter of the templates. It must be LTIR, Mdust, Mstar,
            SFR, sSFR, or Z.
        value: float
            Value of the physical property.

        Returns
        -------
        estimate: array
            Spectrum corresponding to a given value of param
        """
        return (
            self.alpha[param] * value1
            + self.beta[param] * value2
            + self.gamma[param]
        )


def main():
    # Initiate the template builder
    spectra = Spectra()
    spectra2nd = Spectra2nd()

    # Make sure we have the minimum number of arguments
    if len(sys.argv) < 3:
        raise Exception(
            "You must provide the parameter and then a list of values as "
            "argument. The parameters are LTIR (Lsun), Mdust (Msun), Mstar "
            "(Msun), SFR (Msun yr¯¹), sSFR (yr¯¹), and Z. All the values are "
            "in log except for Z, which is in 12+log O/H."
        )

    # Make sure the parameter is valid
    if sys.argv[1] not in spectra.alpha:
        raise Exception(
            f"The parameter {sys.argv[1]} is not recognized. The parameters "
            "are LTIR (Lsun), Mdust (Msun), Mstar (Msun), SFR (Msun yr¯¹), "
            "sSFR (yr¯¹), and Z. All the values are in log except for Z, which "
            "is in 12+log O/H."
        )

    param = sys.argv[1]

    params = {}
    for arg in sys.argv[1:]:
        if arg in ['LTIR', 'Mdust', 'Mstar', 'SFR', 'sSFR', 'Z']:
            current = arg
            params[current] = []
        else:
            try:
                params[current].append(float(arg))
            except ValueError:
                raise Exception(f"Could not convert the input argument {arg}.")

    # Make sure we can interpret correctly the values
    try:
        for param in params:
            params[param] = np.array(params[param], dtype=np.float)
    except ValueError:
        raise Exception(
            f"Could not convert the input arguments {params[param]} to float."
        )

    # Check we do not overwrite an existing gile
    if os.path.exists(fname := f'{"-".join(params.keys())}.fits'):
        raise Exception(
            f"The file {fname} already exists. Please delete it before running "
            "the script again."
        )

    newtable = Table()
    newtable.add_column(Column(spectra.wl, name='wavelength', unit='nm'))
    if len(params) == 1:
        param = list(params.keys())[0]
        for value in params[param]:
            spectrum = Column(
                10 ** spectra.estimate(param, value),
                name=f'nuLnu[{param}={value}]',
                unit='Lsun'
            )
            newtable.add_column(spectrum)
    elif len(params) == 2:
        param1, param2 = list(params.keys())
        values1, values2 = (params[param1], params[param2])
        for value1, value2 in product(values1, values2):
            spectrum = Column(
                10 ** spectra2nd.estimate((param1, param2), value1, value2),
                name=f'nuLnu[{param1}={value1}, {param2}={value2}]',
                unit='Lsun'
            )
            newtable.add_column(spectrum)
    else:
        raise Exception("There can only be one or two parameters.")
    newtable.write(fname)
