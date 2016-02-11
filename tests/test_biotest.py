from biotest import BioTest, MockableFile

class TestAssertFilesEqual(BioTestCase):
    def test_lines_not_sorted_not_equal(self):
        f1 = MockableFile('one\ntwo\nthree')
        f2 = MockableFile('three\ntwo\none')
        self.assertFilesEqual(f1, f2)
