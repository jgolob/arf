#!/usr/bin/env python3
"""
Add is_type column

1. sequences in args.types
2. sequences with ATCC in description
3. sequences that start with NR_ in version
"""

import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'seq_info',
        help='full seq_info file with description column')
    p.add_argument(
        'types',
        type=argparse.FileType('r'),
        help='txt file of version numbers that are type strains')

    p.add_argument(
        'out',
        help='list of all version downloaded')

    args = p.parse_args()

    info = pandas.read_csv(args.seq_info, dtype=str)
    columns = info.columns.tolist()
    if 'is_type' in info.columns:
        info = info.drop('is_type', axis=1)
    else:
        columns.append('is_type')

    info.loc[:, 'is_type'] = False

    ncbi = set(i.strip() for i in args.types if i)
    ncbi_types = info['version'].isin(ncbi)
    info.loc[ncbi_types, 'is_type'] = True

    atcc_types = info['description'].apply(lambda x: 'ATCC' in x)
    info.loc[atcc_types, 'is_type'] = True

    nr_types = info['version'].apply(lambda x: x.startswith('NR_'))
    info.loc[nr_types, 'is_type'] = True

    info.to_csv(args.out, index=False, columns=columns)


if __name__ == '__main__':
    main()
