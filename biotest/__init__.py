import mock

from . import biofile, biotest, biohypothesis

FileMocker = biofile.FileMocker
MockableFile = biofile.MockableFile
BioTestCase = biotest.BioTestCase
seq_record_strategy = biohypothesis.seq_record_strategy
seqrec = biohypothesis.seq_rec_strategy_factory
interleaved_seqrec = biohypothesis.interleaved_strategy_factory

# Simple mocked SeqRecord class
from functools import partial
from Bio.SeqRecord import SeqRecord
MockSeqRecord = partial(mock.MagicMock, SeqRecord)

# Put this here so others can use it
try:
    import __builtin__ as builtins
except:
    import builtins as builtins
