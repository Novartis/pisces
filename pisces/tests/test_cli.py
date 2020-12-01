from unittest import TestCase
from pisces.cli import create_parser


class TestRun(TestCase):
    def setUp(self):
        self.parser = create_parser()

    def test_defaults_config(self):
        "Test that default configuration is set."
        args = self.parser.parse_args(['run', '-fq1', 'data/test1.fastq'])
        assert 'human' in args.conf.keys()
