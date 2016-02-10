from __future__ import print_function
from . import unittest
from operator import attrgetter as attr
from functools import partial
import sh
from toolz import compose
from Bio import SeqIO
import sys
import os
import mock

THISD = os.path.dirname(os.path.abspath(__file__))
here = partial(os.path.join, THISD)

TESTDIR = os.path.dirname(os.path.abspath(__file__))
PROJDIR = os.path.dirname(TESTDIR)




# file mocking? http://www.ichimonji10.name/blog/6/
if sys.version[0] == '2':
        import __builtin__ as builtins  # pylint:disable=import-error
else:
        import builtins  # pylint:disable=import-error
def mock_file(func, read_data, *args, **kwargs):
    with mock.patch.object(builtins, 'open', mock.mock_open(read_data=read_data)): #, create = True) as m:
        with open('_') as handle:
            '''mock_open doesn't iterate properly, so work-around read it here'''
            return func(handle.read().split('\n'), *args, **kwargs)

class BioTest(unittest.TestCase):
    #TODO: avoid transforming filenames to files all the time
    ''' i.e. @as_files  --> converts filenames to files and asserts they exist'''
    def assertFilesEqual(self, fn1, fn2, sort=False, strip=False):
        no_newlines = partial(filter, bool)
        strip_lines = partial(map, str.strip)
        strip_all = compose(no_newlines, strip_lines)
        ''' http://stackoverflow.com/a/3943697 '''
        if not hasattr(fn1, 'read'):
            fh1 = open(fn1)
            fh2 = open(fn2)
        else: fh1, fh2 = fn1, fn2
        tolines = strip_all if strip else list
        lines1, lines2 = tolines(fh1), tolines(fh2)

        if sort:
            lines1, lines2 = sorted(lines1), sorted(lines2)
        self.assertFalse(len(lines1) == 0)
        return self.assertMultiLineEqual('\n'.join(lines1), '\n'.join(lines2))

    def print_file_diff(self, fh1, fh2):
        '''show the diff of two files when a test fails for easy debugging'''
        sh.diff(fh1, fh2)

    #TODO: accept a mock string (and use BytesIO to open)'''
    def seqs_equal(self, fn1, fn2, format):
        open_sorted = compose(partial(sorted, key=compose(str, attr('seq'))), partial(SeqIO.parse, format=format))
        fq1, fq2 = map(open_sorted, [fn1, fn2])
        self.assertFalse(len(fq1) == 0)
        map(self.assert_seq_recs_equal, fq1, fq2)


    def get_record_set(self, fh1):
        pass

    #TODO: allow for ignoring ids, names, etc.
    #TODO: Use FastqGeneralIterator: https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py#L799
    def assert_seq_recs_equal(self, seq1, seq2, quality=False):
        '''This is necessary because the __eq__ in SeqRecord is weird.'''
        _fields = ['id', 'name', 'description']
        seqstr = compose(str, attr('seq'))
        self.assertEquals(seq1._per_letter_annotations, seq2._per_letter_annotations)
        self.assertEquals(seqstr(seq1), seqstr(seq2))
        for field in _fields:
            f1, f2 = getattr(seq1, field), getattr(seq2, field)
            #TODO: include field name in message?
            self.assertEquals(f1, f2) #, msg=)
        #self.seqrecs =  list(SeqIO.parse(BytesIO(self.reads), format='fastq'))

