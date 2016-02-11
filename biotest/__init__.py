import mock

from biofile import FileMocker, MockableFile
from biotest import BioTestCase
from biohypothesis import seq_record_strategy

# Simple mocked SeqRecord class
from functools import partial
from Bio.SeqRecord import SeqRecord
MockSeqRecord = partial(mock.MagicMock, SeqRecord)

# Put this here so others can use it
try:
    import __builtin__ as builtins
except:
    import builtins as builtins
