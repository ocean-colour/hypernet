
# Standard imports
import glob, os
from setuptools import setup, find_packages


# Begin setup
setup_keywords = dict()
setup_keywords['name'] = 'hypernet'
setup_keywords['description'] = 'Analysis of WATER HYPERNET hyperspectral radiometry'
setup_keywords['author'] = 'J. Xavier Prochaska'
setup_keywords['author_email'] = 'jxp@ucsc.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/ocean-colour/hypernet'
setup_keywords['version'] = '0.0.dev0'
# Use README.md as long_description.
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.12'
setup_keywords['install_requires'] = [
    'numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn', 'cartopy',
    'scikit-learn', 'tqdm', 'IPython', 'pytest', 'tabulate',
    # netCDF / data handling (HYPSTAR products are netCDF)
    'xarray', 'netcdf4', 'h5netcdf', 'h5py', 'cftime', 'pyarrow',
    # Data retrieval from the HYPERNETS data portal
    'requests']
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages()

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

setup(**setup_keywords)
