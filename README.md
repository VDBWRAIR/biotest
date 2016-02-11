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
