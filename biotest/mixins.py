'''classes for unittest.TestCase to inherit'''
from os.path import *
from testfixtures import  TempDirectory
testname = splitext(basename(__file__))[0]
testdirbasepath = join(common.TESTDIR, 'testoutput', testname)
class RunsInFixtureDir(object):
    def setUp(self):
        # Start new tempdir
        if exists(testdirbasepath):
            shutil.rmtree(testdirbasepath)
        self.tdir = TempDirectory(path=testdirbasepath)
        os.makedirs(testdirbasepath)
