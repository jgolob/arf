#!/usr/bin/env python3
"""
output all valid tax_ids in database
"""

import argparse
import pandas
import sqlalchemy


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'feather', help='full seq_info file with description column')
    p.add_argument(
        'url',
        help='Database string URI or filename.')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    url = args.url
    engine = sqlalchemy.create_engine(url)
    conn = engine.connect()
    meta = sqlalchemy.MetaData(engine, reflect=True)
    nodes = meta.tables['nodes']
    s = sqlalchemy.select([nodes.c.tax_id]).where(nodes.c.is_valid)
    named = set(i[0] for i in conn.execute(s))
    info['is_valid'] = False
    info.loc[info['tax_id'].isin(named), 'is_valid'] = True
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
