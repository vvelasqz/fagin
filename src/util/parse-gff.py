#!/usr/bin/env python3

import argparse
import sys


class Entry:
    def __init__(self, row):
        self.row = row
        self.attr = {t:v for t,v in (s.split("=") for s in row[8].split(";"))}

    def print_(self, tags=None, split=False, use_ids=False):
        if(tags):
            try:
                if split:
                    nine = "\t".join([self.attr[k] for k in tags])
                else:
                    if use_ids and "Name" in tags and not "Name" in self.attr:
                        self.attr["Name"] = self.attr["ID"]
                    if len(tags) == 1:
                        nine = self.attr[tags[0]]
                    else:
                        nine = ";".join(["%s=%s" % (k, self.attr[k]) for k in tags])
            except KeyError:
                err("Input error: requested tag missing")
        else:
            nine = self.row[8]

        out = "\t".join(self.row[0:8] + [nine])

        print(out)


def err(msg):
    sys.exit(msg)


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'gfffile',
        help="GFF input file",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '-s', '--select',
        help="select where 3rd column matches one of these (1+)",
        nargs='+'
    )
    parser.add_argument(
        '-r', '--reduce',
        help="reduce the attribute column (9th) to these tags",
        nargs='+'
    )
    parser.add_argument(
        '-p', '--split',
        help="split attribute column by tag (in REDUCE order)",
        action="store_true",
        default=False
    )
    parser.add_argument(
        '-m', '--mapids',
        help="map parent ID to parent Name",
        action="store_true",
        default=False
    )
    parser.add_argument(
        '-d', '--use-id-if-unnamed',
        help="If no Name field is present, use ID instead",
        action="store_true",
        default=False
    )
    args = parser.parse_args()

    return(args)

def rowgen(gff):
    entries = []
    for line in gff.readlines():
        if(line[0] == "#"):
            continue

        row = line.rstrip().split("\t")

        if(len(row) != 9):
            err("Bad GFF, must be TAB-delimited with 9 columns")

        entries.append(Entry(row))

    return entries

def mapids(entries):
    idmap = dict()
    for entry in entries:
        try:
            ID = entry.attr['ID']
        except KeyError:
            err("Bad GFF, 9th column must have ID tag")

        try:
            Name = entry.attr['Name']
        except KeyError:
            Name = entry.attr['ID']

        idmap[ID] = Name
    return(idmap)


def parent_id2name(entries, idmap):
    for entry in entries:
        try:
            entry.attr['Parent'] = idmap[entry.attr['Parent']]
        except KeyError:
            pass
        except TypeError:
            pass


if __name__ == '__main__':
    args = parser()

    entries = rowgen(args.gfffile)

    if(args.mapids):
        idmap = mapids(entries)
        parent_id2name(entries, idmap)

    attr_join = "\t" if args.split else None

    if(args.select):
        entries = [e for e in entries if e.row[2] in args.select]

    for entry in entries:
        entry.print_(args.reduce, attr_join, args.use_id_if_unnamed)
