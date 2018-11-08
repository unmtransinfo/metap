#!/usr/bin/env python3

"""
    Usage: uni2mgi.py INPUT OUTPUT

"""

from docopt import docopt
from lxml import objectify


def main(arguments):
    root = objectify.parse(arguments["INPUT"])
    root = root.getroot()
    with open(arguments["OUTPUT"], "w") as ofs:
        ofs.write('ACCESSION\tMGI\n')
        for i in range(len(root.entry)):
            entry = root.entry[i]
            mgi = get_mgi(entry)
            if mgi is None:
                continue
            ofs.write("{0}\t{1}\n".format(entry.accession, mgi))


def get_mgi(entry):
    if entry.dbReference is not None and len(entry.dbReference) > 0:
        for dbref in entry.dbReference:
            if(dbref.attrib["type"] == "MGI"):
                return dbref.attrib["id"]
    return None


if __name__ == "__main__":
    argv = docopt(__doc__, version="1.0")
    main(argv)