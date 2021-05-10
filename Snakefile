configfile: "data/config.yaml"

RF, = glob_wildcards("data/Jacky_171220_guppy_basecalls/{readfile}.fastq.gz")


rule all:
	input:
		"results/fasta_combined.fasta",
		#"results/demfile.txt",
		"results/reference_free_barcoding.done",
		"results/report.pdf"
		
rule create_demfile:
	input:
		config["demfile"]		
	output:
		"results/demfile.txt"
	shell:
		"""
		#tail -n +2 {input} | awk -F ";" '{{printf $1","; printf substr($3,1,13)",";  printf substr($5,1,13)","; printf substr($3,14,38)","; printf substr($5,14,38)"\\n";}}' > {output}
		cp {input} {output}
		"""

rule fastqtofasta:
	input:
		"data/Jacky_171220_guppy_basecalls/{readfile}.fastq.gz"
	output:
		fasta_file = "results/fasta_pass/{readfile}.fa"
	shell:
		"""
		zcat {input} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | sed 's/\\t/\\n/'> {output}
		"""

rule combine_fasta:
	input:
		expand("results/fasta_pass/{rf}.fa", rf=RF) 
	output: 
		"results/fasta_combined.fasta"
	shell:
		"""
		cat {input} > {output}
		"""

rule reference_free_barcoding:
	input:
		fasta = rules.combine_fasta.output,
		demfile = config["demfile"]
	output:
		"results/reference_free_barcoding.done"
	params:
		len = config["barcoding_settings"]["length"],
		depth = config["barcoding_settings"]["depth"]
	singularity:
		"docker://reslp/minibarcoder:5e1dc3b"
	shell:
		"""
		rm -rf results/reffreebarcoding
		cp {input.demfile} results/
		python /software/miniBarcoder/miniBarcoder.py -f {input.fasta} -d {input.demfile} -o results/reffreebarcoding -D {params.depth} -t 8 -l {params.len}
		touch {output}
		"""

rule get_summary_for_plotting:
	input:
		checkpoint = rules.reference_free_barcoding.output,
		fasta = rules.combine_fasta.output
	output:
		summary = "results/barcode_summary.txt",
		total_seqs = "results/total_seqs.txt"
	shell:
		"""
		cat results/reffreebarcoding/all_barcodes.fa | grep ">" | awk -F ";" '{{printf substr($1,2,length($1)-5)","; printf $2","; printf $3"\\n";}}' > {output.summary}
		cat {input.fasta} | grep ">" -c > {output.total_seqs}
		"""
rule report:
	input:
		rules.get_summary_for_plotting.output.summary,
		rules.get_summary_for_plotting.output.total_seqs
	output:
		"results/report.pdf"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/rmarkdown:3.6.3"
	shell:
		"""
		#cp -rf /usr/local/texlive tmp/texlive
		#cp -rf /opt/texlive/ tmp/opt_texlive
		cp -f {params.wd}/data/report.Rmd {params.wd}/results/report.Rmd
		cd results
		Rscript --vanilla -e "library(rmarkdown); rmarkdown::render('report.Rmd')" $(cat fasta_combined.fasta | grep ">" | wc -l)
		"""
