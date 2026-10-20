#!/usr/bin/env python3
"""Compare completed receiver-side TrgTime traces without ignoring order."""
import argparse
import csv
from itertools import zip_longest
from pathlib import Path

HEADER = "# nestdaq-trgtime-v1,time_unit=4ns"
WIDTHS = {"MESSAGE": 4, "META": 7, "TRG": 8, "SLICE": 7, "END": 5}


def records(path):
    counts = dict.fromkeys(("MESSAGE", "META", "TRG", "SLICE"), 0)
    ended = False
    with Path(path).open(newline="") as stream:
        if stream.readline().rstrip("\r\n") != HEADER:
            raise ValueError(f"{path}: missing/unsupported trace header")
        for line, row in enumerate(csv.reader(stream), 2):
            kind = row[0] if row else ""
            if ended or kind not in WIDTHS or len(row) != WIDTHS[kind]:
                raise ValueError(f"{path}:{line}: invalid record: {row}")
            if any(not value.isascii() or not value.isdecimal() for value in row[1:]):
                raise ValueError(f"{path}:{line}: expected unsigned integers")
            values = tuple(map(int, row[1:]))
            if kind == "END":
                if values != tuple(counts.values()):
                    raise ValueError(f"{path}:{line}: END counts disagree with records")
                if not counts["TRG"] or not counts["SLICE"]:
                    raise ValueError(f"{path}: no received triggers/slices; comparison is inconclusive")
                ended = True
            else:
                if kind == "MESSAGE" and values[0] != counts["MESSAGE"]:
                    raise ValueError(f"{path}:{line}: nonconsecutive message index")
                if kind != "MESSAGE" and values[0] != counts["MESSAGE"] - 1:
                    raise ValueError(f"{path}:{line}: wrong message association")
                counts[kind] += 1
            yield (kind, *values)
    if not ended:
        raise ValueError(f"{path}: missing END (run interrupted or not stopped cleanly)")


def compare(cpu, gpu):
    end = None
    for count, (left, right) in enumerate(zip_longest(records(cpu), records(gpu)), 1):
        if left != right:
            raise ValueError(f"Mismatch at record {count}:\n  CPU: {left}\n  GPU: {right}")
        end = left
    return end[1:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpu", type=Path)
    parser.add_argument("gpu", type=Path)
    args = parser.parse_args()
    try:
        messages, metas, triggers, slices = compare(args.cpu, args.gpu)
    except (OSError, ValueError) as error:
        parser.exit(1, f"FAIL: {error}\n")
    print(f"MATCH: messages={messages}, metas={metas}, TrgTime={triggers}, slices={slices}")


if __name__ == "__main__":
    main()
