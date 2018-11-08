#!/usr/bin/env python3

"""
    Usage: uniprot.py INPUT OUTPUT

"""

from docopt import docopt
from lxml import objectify


def main(arguments):
    root = objectify.parse(arguments["INPUT"])
    root = root.getroot()
    with open(arguments["OUTPUT"], "w") as ofs:
        ofs.write('SWISSPROT\tACCESSION\tGENE\tNAME\tORGANISM\tTAX_ID\tENSEMBL_GENE_ID\tENTREZ_GENE_ID\n')
        for i in range(len(root.entry)):
            entry = root.entry[i]
            gene = "NA"
            if len([el.tag for el in entry.iterchildren(tag='{http://uniprot.org/uniprot}gene')]) > 0:
                gene = entry.gene.name
            ofs.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(entry.name, entry.accession, gene,
                                                    entry.protein.recommendedName.fullName,
                                                    get_organism_name(entry),
                                                    get_tax_id(entry),
                                                    get_ensembl_gene_id(entry),
                                                    get_entrez_geneid(entry)))

def get_organism_name(entry):
    if entry.organism.name is not None and len(entry.organism.name) > 0:
        for name in entry.organism.name:
            if name.attrib["type"] == "scientific":
                return name.text
    return None


def get_tax_id(entry):
    if entry.organism.dbReference is not None and len(entry.organism.dbReference) > 0:
        tax = entry.organism.dbReference[0]
        if entry.organism.dbReference.attrib["type"] == "NCBI Taxonomy":
            return entry.organism.dbReference.attrib["id"]
    return None

def get_ensembl_gene_id(entry):
    if entry.dbReference is not None and len(entry.dbReference) > 0:
        for dbref in entry.dbReference:
            if(dbref.attrib["type"] == "Ensembl"):
                if dbref.property is not None and len(dbref.property) > 0:
                    for prop in dbref.property:
                        if prop.attrib["type"] == "gene ID":
                            return prop.attrib["value"]
    return None
    
def get_entrez_geneid(entry):
    if entry.dbReference is not None and len(entry.dbReference) > 0:
        for dbref in entry.dbReference:
            if(dbref.attrib["type"] == "GeneID"):
                return dbref.attrib["id"]
    return None


if __name__ == "__main__":
    argv = docopt(__doc__, version="1.0")
    main(argv)