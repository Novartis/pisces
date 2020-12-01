from unittest import TestCase
from pisces.cli import create_parser


class TestRun(TestCase):

    def test_defaults_config(self):
        "Test that default configuration is set."
        parser = create_parser(['index'])
        args, unknown_args = parser.parse_known_args()
        assert 'human' in args.conf.keys()
