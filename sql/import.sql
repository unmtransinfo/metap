COPY protein(swissprot, accession, symbol, name, species, tax_id, ensembl_gene_id, ncbi_gene_id)
FROM '/Users/oleg/workspace/metap/data/uniprot/human.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY protein(swissprot, accession, symbol, name, species, tax_id, ensembl_gene_id, ncbi_gene_id)
FROM '/Users/oleg/workspace/metap/data/uniprot/mouse.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY protein(swissprot, accession, symbol, name, species, tax_id, ensembl_gene_id, ncbi_gene_id)
FROM '/Users/oleg/workspace/metap/data/uniprot/rat.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

CREATE INDEX gene_ncbi_idx on protein(ncbi_gene_id);

COPY go(go_id, name, definition, aspect)
FROM '/Users/oleg/workspace/metap/data/go/go.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY goa(protein_id, go_id, evidence, assigned_by)
FROM '/Users/oleg/workspace/metap/data/go/goa_human.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY reactome(reactome_id, name, species)
FROM '/Users/oleg/workspace/metap/data/reactome/pathways.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY reactomea(protein_id, reactome_id, evidence)
FROM '/Users/oleg/workspace/metap/data/reactome/reactome_annon.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY homology(ncbi_gene_id, homologene_group_id, tax_id, symbol, protein_gi, ref_seq, protein_id)
FROM '/Users/oleg/workspace/metap/data/homology/homologene.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);

COPY mouse_pheno(marker_accession_id, mp_term_id, mp_term_name, marker_symbol, significant, p_value, phenotype_sex, life_stage_acc, life_stage_name, statistical_method, parameter_name, effect_size, zygosity, protein_id)
FROM PROGRAM 'gunzip -c /Users/oleg/workspace/metap/data/impc/all_api_data.tsv.gz'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);
