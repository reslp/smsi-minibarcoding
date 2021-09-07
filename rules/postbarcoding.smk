rule download_bold_sequences:
	output:
		reference_sequences = "results/{name}/postprocessing/" + config["bold_taxon"] + ".fas"
	params:
		taxon=config["bold_taxon"]
	shell:
		"""
		wget http://v3.boldsystems.org/index.php/API_Public/sequence?taxon={params.taxon} -O {output.reference_sequences}
		"""

rule create_blast_database:
	input:
		reference_sequences = rules.download_bold_sequences.output.reference_sequences 
	output:
		db = "results/{name}/postprocessing/blastdb/{name}_db"
	params:
		name = "{name}"
	singularity: "docker://reslp/ncbi-blast:2.9.0"
	shell:
		"""
		makeblastdb -in {input.reference_sequences} -dbtype 'nucl' -hash_index -out {output.db}
		"""

rule blast_against_bold:
	input:
		db = rules.create_blast_database.output.db,
		checkpoint = rules.reference_free_barcoding.output
	output:
		"results/{name}/postprocessing/{name}_blast_results.xml"
	params:
		name = "{name}"
	singularity: "docker://reslp/ncbi-blast:2.9.0"
	shell:
		"""
		blastn -db {input.db} -query results/{params.name}/reffreebarcoding/all_barcodes.fa -outfmt 5 > {output}
		"""

#rules reorient_sequences:
#	input:
#		blast_results = rules.blast_against_bold.output
#	output:
#		5P_barcodes = "results/{name}/postprocessing/reoriented_barcodes/{name}_all_barcodes_5P.fasta"
#	params:
#		name = "{name}"
#	singularity: "docker://reslp/biopython_plus:1.77"
#	shell:
#		"""
#		bin/reorient_sequences.py
#		"""
	
