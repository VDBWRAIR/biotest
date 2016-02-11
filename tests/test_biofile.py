from . import unittest
import mock

from biotest import FileMocker, MockableFile

class TestFileMockerPatches(unittest.TestCase):
    def test_can_be_used_with_patch(self):
        with mock.patch('__builtin__.open', FileMocker('foo')) as mo:
            x = open('anything.txt')
            self.assertIsInstance(x, MockableFile)
            self.assertEqual(x.read(), 'foo')

class TestFileMockerStoresContents(unittest.TestCase):
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
        print "Testing foo.txt"
        self.assertEqual(x('foo.txt').read(), 'foo')
        print "Testing bar.txt"
        self.assertEqual(x('bar.txt').read(), 'bar')

    def test_raises_ioerror_if_no_file(self):
        x = FileMocker({'foo.txt': 'bar'})
        print x._files
        self.assertRaises(IOError, x, 'bar.txt')

class TestMockableFile(unittest.TestCase):
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
