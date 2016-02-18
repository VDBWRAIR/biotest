from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from hypothesis import strategies as st
from hypothesis import given
from functools import wraps

##############
# Hypothesis #
##############
def seq_qual_strategy_factory(min_length=1, max_length=250, min_qual=0, max_qual=40, alphabet='ATGCN'):
    return st.integers(min_value=min_length, max_value=max_length).flatmap(
        lambda n:
            st.tuples(
                st.text(alphabet=alphabet, min_size=n, max_size=n),
                st.lists(st.integers(min_value=min_qual, max_value=max_qual), min_size=n, max_size=n)
            )
    )

def make_seqrec(id, seq, quals=None):
    '''
    Generate simple Bio.SeqRecord.SeqRecord objects

    seq - can be string or tuple of (dna_string, [quals])
    quals - None or sequence of quals(ints)
    '''
    if isinstance(seq, tuple) and not quals:
        seq, quals = seq
        assert not quals, 'Cannot supply quals when seq is a seq_qual tuple'
    if quals:
        quals = {'phred_quality':quals}
    return SeqRecord(
        Seq(seq, IUPAC.ambiguous_dna), id=id, description='', letter_annotations=quals
    )

def seq_rec_strategy_factory(min_length=1, max_length=250, min_qual=0, max_qual=40, alphabet='ATGCN', idstrat=st.text()):
    '''
    Factory to generate Strategy that build sequence records

    min_length - Minimum sequence length
    max_length - Maximum sequence length
    alphabet - Choices for nucleotides
    min_qual - Minimum quality
    max_qual - Maximum quality
    idstrat - strategy that generates sequence ids
    '''
    return st.integers(min_value=min_length, max_value=max_length).flatmap(
        lambda n:
            st.builds(
                make_seqrec,
                idstrat,
                st.text(alphabet=alphabet, min_size=n, max_size=n),
                st.lists(st.integers(min_value=min_qual, max_value=max_qual), min_size=n, max_size=n)
            )
    )

def seq_record_strategy(*args, **kwargs):
    '''
    Decorator to wrap test function to use SeqRecordStrategy

    '''
    def wrapper(func):
        seq_strategy = seq_rec_strategy_factory(*args, **kwargs)
        @given(seq_strategy)
        @wraps(func)
        def wrap_f(self, rec):
            return func(self, rec)
        return wrap_f
    return wrapper

#interleaved_strategy_factory = st.uuids().map(str).flatmap(
#    lambda id:
#        st.tuples(seqrecstrat(st.shared(st.just(id), key=id)), seqrecstrat(st.shared(st.just(id), key=id))))
#

interleaved_strategy_factory = st.uuids().map(str).flatmap(
    lambda id:
        st.tuples(
            seq_rec_strategy_factory(5, 20, idstrat=st.shared(st.just(id), key=id)), 
            seq_rec_strategy_factory(5, 20, idstrat=st.shared(st.just(id), key=id))))


'''
reads_and_indices = st.integers(min_value=1,max_value=10).flatmap(
    lambda n:
        st.tuples(
            st.integers(min_value=0, max_value=40),
            st.lists(rec, min_size=n, max_size=n).map(
                lambda x:
                    list(interleave(x, x))
            ),
            st.lists(rec, min_size=n, max_size=n),
            st.lists(rec, min_size=n, max_size=n),
       )
)

def make_io_matrix(seqs):
    (_min, reads, i1, i2) = seqs

    def input_file(cmd):
        fn = 'geninput.fastq'
        SeqIO.write(reads, fn, 'fastq')
        return cmd.bake(fn)
    #output gets the cmd object, changes it, runs it, returns input
    def output_stdout(cmd):
        sio = StringIO()
        cmd(_out=sio)
        sio.seek(0)
        return sio
    def output_file(cmd):
        fn = 'genoutput.fastq'
        cmd(o=fn)
        return open(fn)
    inputs = st.one_of(*map(st.just, [input_file, lambda cmd:
      cmd.bake(_in='\n'.join((read.format('fastq')) for read in reads))]))
    outputs = st.one_of(*map(st.just, [output_stdout, output_file]))
    ret =   st.tuples(*(map(st.just, seqs) + [inputs, outputs]))
    return ret


io_matrix = reads_and_indices.flatmap(make_io_matrix)
'''
