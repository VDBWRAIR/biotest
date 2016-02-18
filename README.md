# biotest

[![Build Status](https://travis-ci.org/VDBWRAIR/biotest.svg?branch=master)](https://travis-ci.org/VDBWRAIR/biotest)

Test framework to help make testing bioinformatics easier

## Features

* More easily mock out open calls
* Generic test base classes you can use to provide easy functionality to assert
  common things for sequence type data
* More to come...

## How To Use

Add this project as a dependency in your requirements.txt or setup.py tests_require

### FileMocker

This is a handy base class for any generic File Mocking you may want to do.

You can easily tell FileMocker the file paths you expect that will need to be
mocked or you can just give it a string and it will mock any open calls and that'
string will have the contents you give it.

#### Mock multiple files

```python
from biotest import FileMocker, builtins
import mock

# FileMocker will use the dictionary later when you call open to retrieve
# the correct contents
mock_contents = {'foo.txt': 'foo', 'bar.txt': 'bar'}
with mock.patch(FileMocker(mock_contents)) as mock_open:
    x = open('foo.txt').read() # Produces 'foo'
    x = open('bar.txt').read() # Produces 'bar'

# However, at this time even if a path is relative/absolute to the same file,
# it will be a different instance
# Let's assume that the CWD is /
# Both paths below should represent the same file, however, at this time that
# functonality is not implemented so you will have to watch out for that
mock_contents = {'path/foo.txt': 'foo', '/path/foo.txt': 'bar'}
with mock.patch(FileMocker(mock_contents)) as mock_open:
    x = open('path/foo.txt').read() # Produces 'foo'
    x = open('/path/foo.txt').read() # Produces 'bar' not 'foo'
```

#### Mock single call

You can easily just give it a string if you don't really care about different
filenames having different values

```python
from biotest import FileMocker, builtins
import mock

with mock.patch(FileMocker('foo')) as mock_open:
    x = open('foo.txt').read() # Produces 'foo'
    # you can just keep calling it, same result
    x = open('foo.txt').read() # Produces 'foo'
    # Doesn't even matter what path you use
    x = open('turkey.txt').read() # Produces 'foo'
```

#### Reading a file you have not setup generates exception as expected

```python
from biotest import FileMocker, builtins
import mock

with mock.patch(FileMocker()) as mock_open:
    # Cannot open a file that does not exist
    x = open('foo.txt') # Generates IOError

with mock.patch(FileMocker()) as mock_open:
    # You can open a file that does not exist if in write mode
    x = open('foo.txt', 'w')
```

### MockSeqRecord

For convienience there is a MockSeqRecord class that simply allows you to utilize
the mock.MagicMock class that is spec'd around Bio.SeqRecord.SeqRecord.
That is, it ensures that attributes are correct when you try to access/set them.

```python
from biotest import MockSeqRecord
x = MockSeqRecord()
x.id = 'id'
x.description = 'description
x.seq = 'ATGC # YAY, we don't have to make a Seq instance!
x.foo # Will raise AttributeError as it does not exist in SeqRecord class
```

### BioTestCase

There is a BioTestCase that inherits directly from unittest.TestCase that gives you
some nice functionality to test sequence type data. All you need to do is have your
test classes inherit from it.

```python
from biotest import BioTestCase, MockableFile, MockSeqRecord

class TestSomething(BioTestCase):
    def test_make_sure_files_equal(self):
        f1 = MockableFile('foo.txt', contents='foo')
        f2 = MockableFile('bar.txt', contents='bar')
        self.assertFilesEqual(f1, f2) # Will generate AssertionError since contents are not equal

    def test_make_sure_sequence_records_equal(self):
        s = MockSeqRecord()
        s.id = 'id'
        s.description = 'd'
        s.seq = 'ATGC'
        s.name = 'name'
        s.letter_annotations['phred_quality'] = [40,40,40,40]
        self.assertSeqRecordEqual(s, s) # Will test that s is equal to itself
```

### Hypothesis testing

Hypothesis testing is a new-ish way of testing that allows you to "frame" your tests
such that you are not locking the functionality of your code with unittests.

Read more at https://hypothesis.readthedocs.org

#### seqrecord hypothesis

You can use the seqrec strategy to generate SeqRecord objects with the hypothesis
package.

```python
from biotest import BioTestCase, seqrec
class TestSeqRecord(BioTestCase):
    @given(seqrec())
    def test_something_with_seqrecord(self, record):
        self.assertSeqRecordEqual(record, record)
```

#### seqrecord decorator
To make it a bit easier you can use the test decorator as follows to do the same
thing
```python
from biotest import seq_record_strategy, BioTestCase

class TestSeqRecord(BioTestCase):
    @seq_record_strategy()
    def test_something_with_seqrecord(self, record):
        self.assertSeqRecordEqual(record, record)
```

You can customize the records that get generated with either way by using any of
the following args:

- min_length
  Default: 1
- max_length
  Default: 250
- min_qual
  Default: 0
- max_qual
  Default: 40
- alphabet
  Default: ATGCN

```python
class TestSeqRecord(BioTestCase):
    @seq_record_strategy(min_length=10, max_length=50, min_qual=20, max_qual=30, alphabet='ATGC')
    def test_something_with_seqrecord(self, record):
        self.assertSeqRecordEqual(record, record)
```

#### Interleaved sequence records

Sometimes you may want to test interleaved sequence records
This strategy just gives you a tuple of forward, reverse
Essentially two records from calling `seqrec`, but ensuring the ids are the same
for the forward and reverse records.

```python
from biotest import BioTestCase, interleaved_seqrec
class TestSeqRecord(BioTestCase):
    @given(interleaved_seqrec)
    def test_showing_off_interleaved(self, seqrec)
        f, r = seqrec
        self.assertTrue(f.id, r.id)
```
