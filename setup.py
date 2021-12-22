import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'epiinfo'
AUTHOR = 'John Copeland (ported code to Pyton)'
AUTHOR_EMAIL = 'zfj4@cdc.gov'
URL = 'https://github.com/Epi-Info/Epi-Info-Python-Package'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Epi Info'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'scipy'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
