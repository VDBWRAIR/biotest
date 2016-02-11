from setuptools import setup, find_packages

# Do this otherwise exception because biotest/__init__.py gets loaded
# which tries to import dependencies before they are installed by setup below
import shutil
shutil.copy('biotest/__meta__.py', './')
import __meta__

setup(
    name = __meta__.__projectname__,
    version = __meta__.__release__,
    packages = find_packages(),
    author = __meta__.__authors__,
    author_email = __meta__.__authoremails__,
    description = __meta__.__description__,
    license = "GPLv2",
    keywords = __meta__.__keywords__,
    install_requires = [
        'mock',
        'sh',
        'toolz',
        'hypothesis',
        'biopython',
    ],
)
