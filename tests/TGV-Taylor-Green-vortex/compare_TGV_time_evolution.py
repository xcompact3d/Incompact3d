import sys
import argparse

from numpy import loadtxt, testing


desc = "Compare last output to reference Taylor-Green Vortex output."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-i', '--input', type=str)
parser.add_argument('-r', '--reference', type=str)
args = parser.parse_args()

values = loadtxt(args.input)
reference = loadtxt(args.reference)

all_columns_match = True
for column_idx in range(values.shape[1]):
    print(f"Checking column no {column_idx + 1} of {values.shape[1]}... ", end="")
    try:
        testing.assert_allclose(
            values[:, column_idx], reference[:, column_idx]
        )
        print("OK")
    except AssertionError:
        print("MISMATCH")
        all_columns_match = False

if all_columns_match:
    print("SUCCESS -- output matches reference values.")
else:
    print("MISMATCH -- some columns do not match reference values, see above.")
    sys.exit(1)
