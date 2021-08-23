configfile: "data/config.yaml"

#RF, = glob_wildcards(config["fastq_files"]+"/{readfile}.fastq.gz")


		
rule create_demfile:
	input:
		get_demfile		
	output:
		"results/{name}/demfile.txt"
	shell:
		"""
		#tail -n +2 {input} | awk -F ";" '{{printf $1","; printf substr($3,1,13)",";  printf substr($5,1,13)","; printf substr($3,14,38)","; printf substr($5,14,38)"\\n";}}' > {output}
		cp {input} {output}
		"""

rule fastqtofasta:
	input:
		get_reads_dir
	output:
		fasta_path = directory("results/{name}/fasta_pass/")
	shell:
		"""
		for file in $(ls {input}/*.gz); do
			filename="${{file%.*}}"
			filename=$(basename $filename)	
			gzfile=$(basename $file)
			zcat {input}/$gzfile | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | sed 's/\\t/\\n/'> {output.fasta_path}/$filename.fa
		done
		"""

rule combine_fasta:
	input:
		rules.fastqtofasta.output.fasta_path 
	output: 
		"results/{name}/fasta_combined.fasta"
	shell:
		"""
		files=$(ls {input}/*.fa)
		cat $files > {output}
		"""

rule reference_free_barcoding:
	input:
		fasta = rules.combine_fasta.output,
		demfile = get_demfile
	output:
		"results/{name}/reference_free_barcoding.done"
	params:
		len = config["barcoding_settings"]["length"],
		depth = config["barcoding_settings"]["depth"],
		name = "{name}"
	singularity:
		"docker://reslp/minibarcoder:5e1dc3b"
	threads: config["barcoding_settings"]["threads"]
	shell:
		"""
		rm -rf results/reffreebarcoding
		cp {input.demfile} results/{params.name}/demfile.txt
		python /software/miniBarcoder/miniBarcoder.py -f {input.fasta} -d results/{params.name}/demfile.txt -o results/{params.name}/reffreebarcoding -D {params.depth} -t {threads} -l {params.len}
		touch {output}
		"""
