from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).parents[1]))
from run_trgtime_pipeline import count_input_frames


def record(magic, length, header):
    return struct.pack('<QIHH', magic, length, header, 0) + bytes(length - 16)


class InputCountTest(unittest.TestCase):
    def count(self, data):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / 'input.dat'
            path.write_bytes(data)
            return count_input_frames(path)

    def test_frames_with_file_header_and_trailer(self):
        data = record(0x004b4e53454c4946, 304, 304)
        data += record(0x004d5246454d4954, 128, 24)
        data += record(0x004d5246454d4954, 24, 24)
        data += record(0x004c5254454c4946, 304, 304)
        self.assertEqual(self.count(data), 2)

    def test_filter_prefix_skips_header_only(self):
        prefix = struct.pack('<QIHH', 0x0049474f4c544c46, 9999, 56, 0) + bytes(40)
        self.assertEqual(self.count(prefix + record(0x004d5246454d4954, 24, 24)), 1)

    def test_truncated_or_unknown_records_are_rejected(self):
        tf = record(0x004d5246454d4954, 128, 24)
        for data in (tf[:-1], tf + b'x', record(123, 24, 24), b''):
            with self.subTest(data=data[:16]):
                with self.assertRaises(ValueError):
                    self.count(data)

    def test_invalid_header_length_is_rejected(self):
        with self.assertRaises(ValueError):
            self.count(record(0x004d5246454d4954, 24, 25))


if __name__ == '__main__':
    unittest.main()
