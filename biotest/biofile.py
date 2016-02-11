import mock
import tempfile

class MockableFile(object):
    '''
    Actual class that mocks File object(what is returned from open call)

    Should have the same API as the File class
    '''
    def __init__(self, name, mode='r', bufsize=-1, contents=None):
        if contents:
            mode = 'w+'
        self.tfile = tempfile.SpooledTemporaryFile(0, mode, bufsize)
        if contents:
            self.tfile.write(contents)
            self.tfile.seek(0)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.tfile.close()
    
    def __getattr__(self, attr):
        # Delegate all attributes to the NamedTemporaryFile
        return getattr(self.tfile, attr)

class FileMocker(object):
    '''
    Provide easily mockable object for the open call

    >>> import mock
    >>> with mock.patch('__builtin__.open', MockableFile()) as mock_open:
    ...    mock_open.set_contents('foo')
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo', x
    >>> with mock.patch('__builtin__.open', MockableFile('foo')) as mock_open:
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo', x
    >>> with mock.patch('__builtin__.open', MockableFile()) as mock_open:
    ...    mock_open.set_contents({'foo.txt': 'foo', 'bar.txt': 'bar'})
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo'
    ...    x = open('bar.txt').read()
    ...    assert x == 'bar'
    >>> with mock.patch('__builtin__.open', MockableFile({'foo.txt': 'foo', 'bar.txt': 'bar'})) as mock_open:
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo'
    ...    x = open('bar.txt').read()
    ...    assert x == 'bar'
    '''
    def __init__(self, contents=None):
        self._files = {}
        if contents:
            self.set_contents(contents)

    def __call__(self, path, mode='r', buffering=-1):
        ''' Mocks the open call '''
        if path in self._files:
            contents = self._files.get(path)
        elif '' in self._files:
            contents = self._files['']
        else:
            raise IOError(
                2, "No such file or directory: '{0}'".format(path)
            )
        return MockableFile(
            path, mode, buffering, contents=contents
        )

    def set_contents(self, contents):
        if isinstance(contents, basestring):
            self._files[''] = contents
        elif isinstance(contents, dict):
            self._files.update(contents)
        else:
            raise ValueError('Contents must be a string or dict. {0} was supplied'.format(type(contents)))

