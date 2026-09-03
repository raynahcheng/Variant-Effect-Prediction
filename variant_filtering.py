import argparse

import pandas as pd


def filter_significant(input_path: str, output_path: str, threshold: float = 0.05) -> pd.DataFrame:
    df = pd.read_csv(input_path)

    before = len(df)
    df = df[df["P-Value"] < threshold]
    after = len(df)

    print(f"Loaded {before} rows from {input_path}")
    print(f"Kept {after} rows with P-Value < {threshold}")

    df.to_csv(output_path, index=False)
    print(f"Saved filtered data to {output_path}")

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter variants by P-Value threshold.")
    parser.add_argument("--input", default="data/GRCh38.csv", help="Path to input CSV")
    parser.add_argument("--output", default="data/GRCh38_significant.csv", help="Path to write filtered CSV")
    parser.add_argument("--threshold", type=float, default=0.05, help="P-Value threshold (exclusive)")
    args = parser.parse_args()

    filter_significant(args.input, args.output, args.threshold)
