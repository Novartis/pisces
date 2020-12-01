from unittest import TestCase
from pisces.cli import create_parser


class TestRun(TestCase):

    def test_defaults_config(self):
        "Test that default configuration is set."
        args = create_parser(args='submit -m metadata.csv')
        assert 'human' in args.conf.keys()
