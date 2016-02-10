import mock

class MockableFile(mock.MagicMock):
    '''
    Provide easily mockable object for the open call

    >>> import mock
    >>> with mock.patch('__builtin__.open', MockableFile) as mock_open:
    ...    mock_open.set_contents('foo')
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo', x
    >>> with mock.patch('__builtin__.open', MockableFile('foo')) as mock_open:
    ...    x = open('foo.txt').read()
    ...    assert x == 'foo', x

    #>>> with mock.patch('__builtin__.open', MockableFile) as mock_open:
    #...    mock_open.set_contents('foo.txt', 'foo')
    #...    mock_open.set_contents('bar.txt', 'bar')
    #...    x = open('foo.txt').read()
    #...    assert x == 'foo'
    #...    x = open('bar.txt').read()
    #...    assert x == 'bar'
    '''
    __file_contents = None
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
        return self.__file_contents

    def write(self):
        pass

    @classmethod
    def set_contents(cls, *args, **kwargs):
        # Not correct way to do this
        cls.__file_contents = args[0]
        #cls.__file_contents = 'foo'

