from setuptools import find_packages, setup

entry_points = {
    'console_scripts': ['bosa-estimators=estimators:main',
                        'bosa-templates=templates:main',
                        'bosa-plotter=plotter:main']
}

setup(name="bosa",
      version="2021.06",
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=['astropy', 'matplotlib', 'numpy'],
      setup_requires=['astropy', 'matplotlib', 'numpy'],
      entry_points=entry_points,
      package_data={'': ['data/*']},
      author="Médéric Boquien",
      author_email="mederic.boquien@uantof.cl",
      url="http://pages.iu.edu/~salims/gswlc/",
      description="Spectral templates and estimators from Boquien & Salim (2021)",
      license="GPL",
      keywords="astrophysics, galaxy, infrared"
)
