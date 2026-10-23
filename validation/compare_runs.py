"""Compare two regression runs file by file.

Usage: python compare_runs.py <reference_dir> <new_dir>

Every .npy, .npz, and .csv file under the reference must exist in the new
run with identical contents (NaNs compare equal). Pickled interpolators and
plots are skipped: they are derived from the arrays compared here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

reference, new = (Path(arg) for arg in sys.argv[1:3])
checked, problems = 0, []


def same(first, second):
    first, second = np.asarray(first), np.asarray(second)
    if first.shape != second.shape or first.dtype.kind != second.dtype.kind:
        return False
    if first.dtype.kind in "fc":
        return np.array_equal(first, second, equal_nan=True)
    return np.array_equal(first, second)


for path in sorted(reference.rglob("*")):
    if path.suffix not in {".npy", ".npz", ".csv"} or path.is_dir():
        continue
    relative = path.relative_to(reference)
    other = new / relative
    if not other.exists():
        problems.append(f"missing: {relative}")
        continue
    checked += 1
    if path.suffix == ".npy":
        ok = same(np.load(path), np.load(other))
    elif path.suffix == ".csv":
        ok = pd.read_csv(path).equals(pd.read_csv(other))
    else:
        with np.load(path, allow_pickle=False) as a, \
                np.load(other, allow_pickle=False) as b:
            # The new code may add keys (e.g. a dimension record); every
            # reference key must match exactly
            missing = set(a.files) - set(b.files)
            ok = not missing and all(same(a[key], b[key]) for key in a.files)
    if not ok:
        problems.append(f"differs: {relative}")

print(f"Compared {checked} files")
print("\n".join(problems) if problems else "All identical")
sys.exit(1 if problems else 0)
