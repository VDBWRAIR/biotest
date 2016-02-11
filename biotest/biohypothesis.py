from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from hypothesis import strategies as st
from hypothesis import given
from functools import wraps

##############
# Hypothesis #
##############
make_seqrec = lambda id, seq, quals: \
                SeqRecord(Seq(seq, IUPAC.ambiguous_dna), id=str(id), description='', letter_annotations={'phred_quality':quals})

def seq_rec_strategy_factory(min_length=1, max_length=250, min_qual=0, max_qual=40, alphabet='ATGCN'):
    '''
    Factory to generate Strategy that build sequence records
    '''
    return st.integers(min_value=min_length, max_value=max_length).flatmap(
        lambda n:
            st.builds(
                make_seqrec,
                st.integers(),
                st.text(alphabet=alphabet, min_size=n, max_size=n),
                st.lists(st.integers(min_value=min_qual, max_value=max_qual), min_size=n, max_size=n)
            )
    )

def seq_record_strategy(*args, **kwargs):
    '''
    Decorator to wrap test function to use SeqRecordStrategy

    min_length - Minimum sequence length
    max_length - Maximum sequence length
    alphabet - Choices for nucleotides
    min_qual - Minimum quality
    max_qual - Maximum quality
    '''
    def wrapper(func):
        seq_strategy = seq_rec_strategy_factory(*args, **kwargs)
        @given(seq_strategy)
        @wraps(func)
        def wrap_f(self, rec):
            return func(self, rec)
        return wrap_f
    return wrapper
