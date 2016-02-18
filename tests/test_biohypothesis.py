from biotest import biohypothesis, BioTestCase, seqrec
from hypothesis import strategies as st
from hypothesis import given

class TestMakeSeqRecord(BioTestCase):
    @given(st.text(), biohypothesis.seq_qual_strategy_factory())
    def test_makes_record_with_quality_data(self, id, seq_qual):
        seq, qual = seq_qual
        x = biohypothesis.make_seqrec(id, seq, qual)
        self.assertEqual(id, x.id)
        self.assertEqual(seq, str(x.seq))
        self.assertEqual(qual, x.letter_annotations['phred_quality'])

    @given(st.text(), st.text(alphabet='ATGCN'))
    def test_makes_record_without_quality_data(self, id, seq):
        x = biohypothesis.make_seqrec(id, seq)
        self.assertEqual(id, x.id)
        self.assertEqual(seq, str(x.seq))

class TestSeqRecStrategyFactory(BioTestCase):
    @given(biohypothesis.seq_rec_strategy_factory())
    def test_generates_seqecs(self, rec):
        self.assertTrue(hasattr(rec, 'seq'))

    @given(seqrec())
    def test_generates_seqecs_import(self, rec):
        self.assertTrue(hasattr(rec, 'seq'))

    @biohypothesis.seq_record_strategy()
    def test_decorator_generates_seqrecs(self, rec):
        self.assertTrue(hasattr(rec, 'seq'))

class TestInterleavedStrategyFactory(BioTestCase):
    @given(st.lists(biohypothesis.interleaved_strategy_factory(), max_size=20))
    def test_ensure_ids_same_data_diff(self, interleave):
        ids = list()
        for f,r in interleave:
            self.assertEqual(f.id, r.id)
            ids.append(f)
        # Just make sure we are making different id's for each pair(for sanity)
        self.assertNotEqual(1, set(ids))
