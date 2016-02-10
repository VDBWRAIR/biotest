from setuptools import setup, find_packages

import biotest

setup(
    name = biotest.__projectname__,
    version = biotest.__release__,
    packages = find_packages(),
    author = biotest.__authors__,
    author_email = biotest.__authoremails__,
    description = biotest.__description__,
    license = "GPLv2",
    keywords = biotest.__keywords__,
)
