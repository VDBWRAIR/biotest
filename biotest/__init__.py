__version__ = '0.1.0'
__release__ = __version__ + '-dev'
__authors__ = 'Tyghe Vallard, Michael panciera'
__authoremails__ = 'vallardt@gmail.com, michael.panciera.work@gmail.com'
__description__ = 'Testing framework for python bioinformatics'
__projectname__ = 'biotest'
__keywords__ = "bioinformatics, unittest, framework, test"

from biotest.biofile import FileMocker, MockableFile

# Put this here so others can use it
try:
    import __builtin__ as builtins
except:
    import builtins as builtins
