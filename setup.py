import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '1.2.0.14'
PACKAGE_NAME = 'epiinfo'
AUTHOR = 'John Copeland'
AUTHOR_EMAIL = 'ita3@cdc.gov'
URL = 'https://github.com/Epi-Info/Epi-Info-Python-Package'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Epi Info: Import and analyze data'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "Imports data from different formats, including Epi Info sync files, and performs Epi Info statistical analyses."

INSTALL_REQUIRES = [
      'scipy',
      'ijson',
      'pycryptodome'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type="text/markdown",
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
