import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1]))
from run_trgtime_pipeline import TraceProgress


class ProgressTest(unittest.TestCase):
    def test_partial_line_does_not_acknowledge(self):
        progress = TraceProgress()
        progress.feed('MESSAGE,0,2,1\nMETA,0,0,42,42,1,1\nTRG,0,0,42,0,0,123,0\nSLICE,0,0,42,0,123,0')
        self.assertFalse(progress.complete(1))
        progress.feed('\n')
        self.assertTrue(progress.complete(1))
        self.assertFalse(progress.complete(2))

    def test_zero_trigger_meta_is_complete(self):
        progress = TraceProgress()
        progress.feed('MESSAGE,0,1,1\nMETA,0,0,42,42,0,1\n')
        self.assertTrue(progress.complete(1))

    def test_message_header_is_not_completion(self):
        progress = TraceProgress()
        progress.feed('MESSAGE,0,2,1\n')
        self.assertFalse(progress.complete(1))

    def test_missing_trigger_is_not_completion(self):
        progress = TraceProgress()
        progress.feed('MESSAGE,0,2,1\nMETA,0,0,42,42,1,1\nSLICE,0,0,42,0,123,0\n')
        self.assertFalse(progress.complete(1))
        with self.assertRaisesRegex(RuntimeError, 'incomplete'):
            progress.feed('MESSAGE,1,1,1\n')


if __name__ == '__main__':
    unittest.main()
