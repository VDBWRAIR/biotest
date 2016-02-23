import string

from biotest import biohypothesis, BioTestCase, seqrec
from hypothesis import strategies as st
from hypothesis import given, assume
from nose.plugins.attrib import attr

# This makes sure valid id for SeqRecord as SeqRecord doesn't test
# for these things(yay hypothesis)
#assume(not set(['','\n','>',' ','\r']).intersection(set(id)))
ascii_text = st.text(alphabet=string.printable)

@attr('py27+')
class TestMakeSeqRecord(BioTestCase):
    @given(ascii_text, biohypothesis.seq_qual_strategy_factory())
    def test_makes_record_with_quality_data(self, id, seq_qual):
        assume(not set(['','\n','>',' ','\r']).intersection(set(id)))
        seq, qual = seq_qual
        x = biohypothesis.make_seqrec(id, seq, qual)
        self.assertEqual(id, x.id)
        self.assertEqual(seq, str(x.seq))
        self.assertEqual(qual, x.letter_annotations['phred_quality'])
        qualstr = ''.join([chr(j+33) for j in qual])
        self.assertEqual('@{0}\n{1}\n+\n{2}\n'.format(id, seq, qualstr), x.format('fastq'))

    @given(ascii_text, st.text(alphabet='ATGCN', min_size=1))
    def test_makes_record_without_quality_data(self, id, seq):
        assume(not set(['','\n','>',' ','\r']).intersection(set(id)))
        x = biohypothesis.make_seqrec(id, seq)
        self.assertEqual(id, x.id)
        self.assertEqual(seq, str(x.seq))
        self.assertNotIn('phred_quality', x.letter_annotations)
        self.assertEqual('>{0}\n{1}\n'.format(id, seq), x.format('fasta'))

@attr('py27+')
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

@attr('py27+')
class TestInterleavedStrategyFactory(BioTestCase):
    @given(st.lists(biohypothesis.interleaved_strategy_factory(), max_size=20))
    def test_ensure_ids_same_data_diff(self, interleave):
        ids = list()
        for f,r in interleave:
            self.assertEqual(f.id, r.id)
            ids.append(f)
        # Just make sure we are making different id's for each pair(for sanity)
        self.assertNotEqual(1, set(ids))

@attr('py27+')
class TestVCFStrategy(BioTestCase):
    @given(biohypothesis.vcf_dict_strategy_factory('chr1', 1, 'A'))
    def test_ensure_useful_dict(self, vcfrec):
        self.assertEqual(1, vcfrec['pos'])
        self.assertEqual('chr1', vcfrec['chrom'])
        self.assertEqual('A', vcfrec['ref'])
        if type(vcfrec['AO']) == int:
             self.assertGreater(vcfrec['DP'], vcfrec['AO'])
        else:
            self.assertGreater(vcfrec['DP'], sum(vcfrec['AO']))

    @given(biohypothesis.ref_with_vcf_dicts_strategy_factory())
    def test_ensure_useful_records(self, seq_vcfs):
        seq, vcfs = list(seq_vcfs[0]), list(seq_vcfs[1])
        self.assertGreaterEqual(len(seq), len(vcfs))
        # Assert all vcf ref seq chunks are same as on actual reference sequence
        # at specified position
        for vcf in vcfs:
            r = vcf['ref']
            p = vcf['pos']
            refseq = ''.join(seq[p-1:p+len(r)-1])
            self.assertEqual(refseq, r)
