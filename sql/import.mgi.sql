COPY protein2mgi(protein_id, mgi_id)
FROM '/home/oleg/workspace/metap/data/uniprot/protein2mgi.tsv'
WITH (
	FORMAT CSV,
	HEADER true,
	NULL 'None',
	DELIMITER E'\t',
	QUOTE '"'
);