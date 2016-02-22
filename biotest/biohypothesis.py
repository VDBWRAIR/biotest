from __future__ import print_function
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import operator
from fn.iters import accumulate
from toolz.itertoolz import partition
import string
from functools import partial
from biotest.compat import imap, ifilter, takewhile, izip 
from hypothesis import strategies as st
from hypothesis import given, assume
from functools import wraps
from toolz import compose
from fn import Stream
from pyparsing import ParseException, Literal,  Word, \
   alphas, nums, quotedString, delimitedList, removeQuotes
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

def interleaved_strategy_factory():
    '''
    Generate interleaved fastq that guarantees ids are same for pairs
    *_kwargs are supplied to gen seq_rec_strategy_factory
    to customize forward and reverse reads
    '''
    strategy = st.uuids().map(str).flatmap(
        lambda id:
            st.tuples(
                seq_rec_strategy_factory(5, 20, idstrat=st.shared(st.just(id), key=id)),
                seq_rec_strategy_factory(5, 20, idstrat=st.shared(st.just(id), key=id))))
    return strategy

def rolling_sum(elem_min, elem_max, length):
    return st.lists(
        st.integers(min_value=elem_min,max_value=elem_max), max_size=length, min_size=length
    ).map(lambda xs: accumulate(xs, operator.add))

@st.composite
def ref_with_vcf_dicts_strategy_factory(draw):
    '''
    first a generator for increasing size numbers and take slices with them.
    end up with rolling chunks of the reference of parameterizable (but random) size.
    sample some of those and those be come the REF field.
    ALT field is whatever
    POS field is limited by size of reference
    chrom, like POS, is same accross all reference
    AO < DP
    '''
    seq = draw(st.text(alphabet='ACGT', min_size=10, max_size=20))
    size = len(seq)
    ranges = draw(rolling_sum(1, 3, size/2).map(lambda xs: ifilter(lambda x: x < size, xs)) )#.filter(_not(bool)))
    # Stream lets you re-use a generator without draining it.
    pairs = Stream() << partition(2, ranges)
    POSs = Stream() << imap(operator.itemgetter(0), pairs)
    # VCF files start at index 1; python starts at 0
    pairs_offset_1 = imap(lambda x: (x[0] - 1, x[1] - 1), pairs)
    #grab the pieces of the reference
    chunks = map(lambda x: seq[x[0]:x[1]], pairs_offset_1)
    #random chromosome name
    chrom = draw(st.text(string.ascii_letters))
    vcfs = map(compose(draw, partial(vcf_dict_strategy_factory, chrom)), POSs, chunks)
    #TODO: ranges must be non-empty. Assuming vcfs for now.
    assume(len(vcfs) > 0)
    return (seq, vcfs)

@st.composite
def vcf_dict_strategy_factory(draw, chrom, pos, ref):
    '''a generator that returns a single
    VCF dict at a certain position or w/e for testing `call_base`'''
    alts = draw(st.lists(st.text(alphabet='ACGT', min_size=0, max_size=6)), min_size=1, max_size=3)
    #an AO (alternate base count) of 0 doesn't make sense
    ao = draw(st.integers(min_value=1))
    dp = ao + draw(st.integers(min_value=1))
    fields = ['alt', 'ref', 'pos', 'chrom', 'DP', 'AO']
    values = [alts, ref, pos, chrom, dp, ao]
    if None in values:
        raise ValueError("Found none value, did the mapping fail?")
    return dict(zip(fields, values))

#TODO: now consensus should throw error if AO > DP
def parse_header(lines):
    def schema2strategies(schema):
        types = {
            'Integer' : st.integers(),
            'Float' : st.floats(),
            'String' : st.text(max_size=10)
        }
        strategy = types[schema['Type']]
        return schema['ID'], strategy
    header =  takewhile(lambda x: x[0] == '#', lines)
    isMedata = lambda x: x.startswith("##INFO") or x.startswith("##FORMAT")
    metadata = ifilter(isMedata, header)
    result = imap(compose(schema2strategies, parse_header_line), metadata)
    return st.fixed_dictionaries(dict(result))

#strategy_from_header = compose(st.fixed_dictionaries, parse_header)
#TODO: A way to add arbitrary contraints (< 1), and contraints based on relations (AO < DP). use `assume` ?
# some of the
def parse_header_line(lineString):
    lineName = Literal("##").suppress() + (Literal("FORMAT") | Literal("INFO"))
    sentence = quotedString(r'"').setParseAction(removeQuotes)
    def make_kv(key, valParser):
        return Literal(key) + Literal("=").suppress() + valParser
    keyVal = make_kv("ID", Word(alphas + nums + '.')) | make_kv("Type", (Literal("Float") | Literal("String") | Literal("Integer"))) | make_kv("Description", sentence) | make_kv("Number",  Word(alphas + nums))
    fields = delimitedList(keyVal, ",")
    line = lineName + Literal("=<").suppress() + fields + Literal(">").suppress()
    pairs = lambda xs: [] if len(xs) == 0 else [(xs[0], xs[1])] + pairs(xs[2:])
    try:
        res = line.leaveWhitespace().parseString(lineString)
    except ParseException as e:
        print(lineString)
        raise e
    metadata_type = res[0] # e.g. INFO/FORMAT
    schema = dict(pairs(res[1:]))
    return schema

    #return metadata_type, schema
    #return update_in(d['INFO'], ['Type'], compose(types.__getitem__, str.lower))
    #res = line.leaveWhitespace().parseString(r'##INFO=<ID=technology.ILLUMINA,Number=A,Type=Float,Description="Fraction of observations supporting the alternate observed in reads from ILLUMINA">')

# could create a generator for CIGAR strings at some point
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


#st.integers(min_value=1, max_value=3).flatmap(lambda ao:  #not possible with jsut map, need flatmap
#  st.integers(min_value=0, max_value=2).map(lambda dp:
#     (ao, ao+dp))).example()
'''
