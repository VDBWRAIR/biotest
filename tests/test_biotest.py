from nose.tools import assert_raises
from biotest import BioTestCase, MockableFile, MockSeqRecord, seq_record_strategy

def assertHasAttr(obj, attr):
    assert hasattr(obj, attr), "{0} does not have attribute {1}".format(obj, attr)

def test_mockseqrecord_sanity():
    x = MockSeqRecord()
    assertHasAttr(x, 'seq')
    assertHasAttr(x, 'letter_annotations')
    assert_raises(AttributeError, getattr, x, 'foo')

class TestAssertFilesEqual(BioTestCase):
    def test_lines_not_sorted_not_equal(self):
        f1 = MockableFile('foo.txt', contents='one\ntwo\nthree')
        f2 = MockableFile('bar.txt', contents='three\ntwo\none')
        self.assertRaises(AssertionError, self.assertFilesEqual, f1, f2)

    def test_lines_not_stripped_not_equal(self):
        f1 = MockableFile('foo.txt', contents='one\ntwo\nthree')
        f2 = MockableFile('bar.txt', contents='three\n\ntwo\n\none')
        self.assertRaises(AssertionError, self.assertFilesEqual, f1, f2)

    def test_empty_files_raises_error(self):
        f1 = MockableFile('foo.txt')
        f2 = MockableFile('bar.txt')
        self.assertRaises(AssertionError, self.assertFilesEqual, f1, f2)


class TestAssertSeqRecordEqual(BioTestCase):
    @seq_record_strategy
    def test_seqrecord_is_identical_passes(self, rec):
        self.assertSeqRecordEqual(rec, rec)

    def test_make_sure_sequence_records_equal(self):
        s = MockSeqRecord()
        s.id = 'id'
        s.description = 'd'
        s.seq = 'ATGC'
        s.name = 'name'
        s.letter_annotations['phred_quality'] = [40,40,40,40]
        self.assertSeqRecordEqual(s, s) # Will test that s is equal to itself
