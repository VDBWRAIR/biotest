import mock

class MockableFile(mock.MagicMock):
    '''
    Provide easily mockable object for the open call

    >>> import mock
    >>> print globals()['__builtins__']

    >>> with mock.patch.object(__builtins__, 'open', MockableFile) as mock_open:
    ...    mock_open.set_contents('foo')
    ...    print open
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo'

    #>>> with mock.patch(open, MockableFile('foo') as mock_open:
    #...    x = open('foo.txt').read()
    #...    assert x == 'foo'
    #>>> with mock.patch(open, MockableFile) as mock_open:
    #...    mock_open.set_contents('foo.txt', 'foo')
    #...    mock_open.set_contents('bar.txt', 'bar')
    #...    x = open('foo.txt').read()
    #...    assert x == 'foo'
    #...    x = open('bar.txt').read()
    #...    assert x == 'bar'
    '''
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        pass

    def __exit__(self):
        pass

    def __call__(self, *args, **kwargs):
        return self

    def readlines(self):
        pass

    def readline(self):
        pass

    def read(self):
        pass

    def write(self):
        pass

    @classmethod
    def set_contents(cls, *args, **kwargs):
        pass

