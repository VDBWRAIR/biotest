import mock

class MockableFile(mock.MagicMock):
    '''
    Provide easily mockable object for the open call

    import mock
    from biotest import MockableFile

    with mock.patch(open, MockableFile) as mock_open:
        mock_open.set_contents('foo')
        x = open('foo.txt').read()
        assert x == 'foo'

    with mock.patch(open, MockableFile('foo') as mock_open:
        x = open('foo.txt').read()
        assert x == 'foo'

    with mock.patch(open, MockableFile) as mock_open:
        mock_open.set_contents('foo.txt', 'foo')
        mock_open.set_contents('bar.txt', 'bar')
        x = open('foo.txt').read()
        assert x == 'foo'
        x = open('bar.txt').read()
        assert x == 'bar'
    '''
