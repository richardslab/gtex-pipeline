#!/usr/bin/env python3
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Extract variant-gene pairs from list of associations.')
parser.add_argument('--input_pairs', help="output from FastQTL.")
parser.add_argument('--prefix', help='Prefix for output file: <prefix>.extracted_qtls.txt')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

print('Loading input')
ref_pairs_df = pd.read_csv(args.input_pairs, sep='\t', usecols=['variant_id', 'gene_id'], dtype=str)
ref_pairs_df[["sid_chr", "sid_pos", "sid_ref", "sid_alt"]] = ref_pairs_df['variant_id'].str.split(":", expand=True)
ref_pairs_df = ref_pairs_df.rename(columns={"variant_id": "sid", "gene_id": "pid"}).reindex(
    columns=["sid", "pid", "sid_chr", "sid_pos"])[["sid", "pid", "sid_chr", "sid_pos"]]

with open(os.path.join(args.output_dir, args.prefix + '.extracted_qtls.txt'), 'wt') as f:
    ref_pairs_df.to_csv(f, sep='\t', na_rep='NA', float_format='%.6g', index=False)
