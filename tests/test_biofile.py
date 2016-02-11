import mock

from biotest import FileMocker, MockableFile, builtins, BioTestCase

class TestFileMockerPatches(BioTestCase):
    def test_can_be_used_with_patch(self):
        with mock.patch.object(builtins, 'open', FileMocker('foo')) as mo:
            x = open('anything.txt')
            self.assertIsInstance(x, MockableFile)
            self.assertEqual(x.read(), 'foo')

class TestFileMockerStoresContents(BioTestCase):
    def test_set_single_string_contents_via_set_contents(self):
        x = FileMocker()
        x.set_contents('foo')
        self.assertEqual(x('anything.txt').read(), 'foo')

    def test_set_single_string_contents_via_init(self):
        x = FileMocker('foo')
        self.assertEqual(x('anything.txt').read(), 'foo')

    def test_set_multiple_file_via_set_contents(self):
        x = FileMocker()
        x.set_contents({'foo.txt': 'foo', 'bar.txt': 'bar'})
        self.assertEqual(x('foo.txt').read(), 'foo')
        self.assertEqual(x('bar.txt').read(), 'bar')

    def test_set_multiple_file_via_init(self):
        x = FileMocker({'foo.txt': 'foo', 'bar.txt': 'bar'})
        self.assertEqual(x('foo.txt').read(), 'foo')
        self.assertEqual(x('bar.txt').read(), 'bar')

    def test_raises_ioerror_if_no_file(self):
        x = FileMocker({'foo.txt': 'bar'})
        self.assertRaises(IOError, x, 'bar.txt')

    def test_file_used_multiple_times(self):
        x = FileMocker({'foo.txt': 'foo'})
        self.assertEqual(x('foo.txt').read(), 'foo')
        self.assertEqual(x('foo.txt').read(), 'foo')

    def test_support_readline(self):
        x = FileMocker({'foo.txt': 'foo'})('foo.txt')
        self.assertEqual(x.readline(), 'foo')
        self.assertEqual(x.readline(), '')

    def test_supports_filepaths(self):
        x = FileMocker({'/path/foo.txt': 'foo', 'path/bar.txt': 'bar'})
        self.assertEqual(x('/path/foo.txt').read(), 'foo')
        self.assertEqual(x('path/bar.txt').read(), 'bar')

    def test_returns_mockablefile_if_missing_name_but_write(self):
        x = FileMocker()
        try:
            x('foo.txt', 'w')
        except IOError:
            self.fail("raised IOError even though 'w' mode")

class TestMockableFile(BioTestCase):
    def test_can_set_contents_in_init(self):
        x = MockableFile('foo.txt', contents='foo')
        self.assertEqual(x.read(), 'foo')

    def test_reads_and_writes(self):
        x = MockableFile('foo.txt', 'w+')
        x.write('foo')
        x.seek(0)
        self.assertEqual(x.read(), 'foo')

    def test_support_context_manager(self):
        with MockableFile('foo.txt', 'w+') as fh:
            fh.write('foo')
            fh.seek(0)
            self.assertEqual(fh.read(), 'foo')
