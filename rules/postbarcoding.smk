rule download_bold_sequences:
	output:
		reference_sequences = "results/{name}/postprocessing/" + config["post_barcoding"]["bold_taxon"] + ".fas"
	params:
		taxon=config["post_barcoding"]["bold_taxon"]
	shell:
		"""
		wget http://v3.boldsystems.org/index.php/API_Public/sequence?taxon={params.taxon} -O {output.reference_sequences}
		"""

rule create_blast_database:
	input:
		reference_sequences = rules.download_bold_sequences.output.reference_sequences 
	output:
		db = "results/{name}/postprocessing/blastdb/{name}_db.done"
	params:
		name = "{name}",
		wd = os.getcwd()
	singularity: "docker://reslp/ncbi-blast:2.9.0"
	shell:
		"""
		makeblastdb -in {input.reference_sequences} -dbtype 'nucl' -hash_index -out {params.wd}/results/{params.name}/postprocessing/blastdb/{params.name}_db 
		touch {output.db}
		"""

rule blast_consensus_barcode_against_bold:
	input:
		db = rules.create_blast_database.output.db,
		checkpoint = rules.reference_free_barcoding.output
	output:
		"results/{name}/postprocessing/{name}_blast_results.xml"
	params:
		name = "{name}",
		wd = os.getcwd()
	threads: config["post_barcoding"]["blast_threads"]
	singularity: "docker://reslp/ncbi-blast:2.9.0"
	shell:
		"""
		blastn -db {params.wd}/results/{params.name}/postprocessing/blastdb/{params.name}_db -query results/{params.name}/reffreebarcoding/all_barcodes.fa -outfmt 5 -num_threads {threads} > {output}
		"""

rule blast_individual_barcode_against_bold:
	input:
		db = rules.create_blast_database.output.db,
		checkpoint = rules.reference_free_barcoding.output,
		barcode_file = "results/{name}/reffreebarcoding/temp_alls_500/{barcode}_all.fa"
	output:
		"results/{name}/postprocessing/individual_barcodes/blast_results/{barcode}_blast_results.xml"
	params:
		name = "{name}",
		wd = os.getcwd()
	threads: config["post_barcoding"]["blast_threads"]
	singularity: "docker://reslp/ncbi-blast:2.9.0"
	shell:
		"""
		blastn -db {params.wd}/results/{params.name}/postprocessing/blastdb/{params.name}_db -query {input.barcode_file} -outfmt 5 -num_threads {threads} > {output}
		"""

rule reorient_individual_barcodes:
	input:
		blast_results = rules.blast_individual_barcode_against_bold.output,
		barcode_file = "results/{name}/reffreebarcoding/temp_alls_500/{barcode}_all.fa"
	output:
		P5_barcodes = "results/{name}/postprocessing/individual_barcodes/reoriented/{barcode}_all_barcodes_5P.fasta"
	params:
		name = "{name}",
		wd = os.getcwd()
	singularity: "docker://reslp/biopython_plus:1.77"
	shell:
		"""
		bin/reorient_sequences.py --blast_results {input.blast_results} --sequence_file {input.barcode_file} > {output.P5_barcodes} 
		"""

rule reorient_consensus_barcodes:
	input:
		blast_results = rules.blast_consensus_barcode_against_bold.output
	output:
		P5_barcodes = "results/{name}/postprocessing/reoriented_barcodes/{name}_all_barcodes_5P.fasta"
	params:
		name = "{name}",
		wd = os.getcwd()
	singularity: "docker://reslp/biopython_plus:1.77"
	shell:
		"""
		bin/reorient_sequences.py --blast_results {input.blast_results} --sequence_file {params.wd}/results/{params.name}/reffreebarcoding/all_barcodes.fa > {output.P5_barcodes} 
		"""
	

