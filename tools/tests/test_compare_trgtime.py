import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("compare_trgtime", Path(__file__).parents[1] / "compare_trgtime.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

TRACE = """# nestdaq-trgtime-v1,time_unit=4ns
MESSAGE,0,2,1
META,0,0,42,42,1,1
TRG,0,0,42,0,0,123,2852126720
SLICE,0,0,42,0,123,2852126720
END,1,1,1,1
"""


class CompareTest(unittest.TestCase):
    def compare(self, cpu=TRACE, gpu=TRACE):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "cpu.csv", Path(directory) / "gpu.csv"
            a.write_text(cpu)
            b.write_text(gpu)
            return module.compare(a, b)

    def test_match(self):
        self.assertEqual(self.compare(), (1, 1, 1, 1))

    def test_mismatch_in_time_type_or_frame(self):
        for before, after in [("123", "124"), ("2852126720", "0"), ("42", "43")]:
            with self.subTest(before=before), self.assertRaisesRegex(ValueError, "Mismatch"):
                self.compare(gpu=TRACE.replace(before, after))

    def test_truncated_identical_traces_are_not_a_match(self):
        trace = TRACE[:TRACE.index("END")]
        with self.assertRaisesRegex(ValueError, "missing END"):
            self.compare(trace, trace)

    def test_empty_runs_are_inconclusive(self):
        trace = module.HEADER + "\nEND,0,0,0,0\n"
        with self.assertRaisesRegex(ValueError, "inconclusive"):
            self.compare(trace, trace)

    def test_bad_counts_are_rejected(self):
        trace = TRACE.replace("END,1,1,1,1", "END,1,1,2,1")
        with self.assertRaisesRegex(ValueError, "counts disagree"):
            self.compare(trace, trace)

    def test_extra_records_after_end_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "invalid record"):
            self.compare(TRACE + "MESSAGE,1,0,0\n", TRACE + "MESSAGE,1,0,0\n")


if __name__ == "__main__":
    unittest.main()
