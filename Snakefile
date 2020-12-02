RF, = glob_wildcards("data/fast5_pass/{readfile}.fast5")

rule all:
	input:
		"results/fastq_files/flappie.done"
		#"results/demfile.txt",
		#"results/reference_free_barcoding.done",
		#"results/report.pdf"
		

rule flappie:
	input:
		dir = "data/fast5_pass"
	output:
		checkpoint = "results/fastq_files/flappie.done"
	log:
		stdout = "log/flappie.out",
		stderr = "log/flappie.err"
	params:
		wd = os.getcwd(),
		model = "r941_native"
	singularity: "docker://reslp/flappie:4de542f"
	threads: 10
	shell:
		"""
			cd results/fastq_files/
			export OPENBLAS_NUM_THREADS={threads}
			for file in $(ls {params.wd}/{input.dir}); do
				filename=$(echo $(basename $file) | awk -F'.' '{{print $1}}')	
				if [ ! -f $filename.fastq.gz ]; then
					echo -e "\\n$(date) - Processing file $filename"			
					echo -e "$(date) - converting fast5"
					if [ -d $filename ]; then rm -rf $filename/*; else mkdir $filename; fi
                                	multi_to_single_fast5 -i {params.wd}/{input.dir}/$file -t {threads} -s $filename/
					echo -e "$(date) - basecalling"
					flappie --model={params.model} $filename/0 | gzip > $filename.fastq.gz.tmp
					mv $filename.fastq.gz.tmp $filename.flappie.fastq.gz
					rm -rf $filename
				fi
				echo -e "$(date) - Done"
			done 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
			touch {params.wd}/{output.checkpoint}
	"""

			

rule create_demfile:
	input:
		"data/Primer_combination_MinION_Teil1.csv"
	output:
		"results/demfile.txt"
	shell:
		"""
		tail -n +2 {input} | awk -F ";" '{{printf $1","; printf substr($4,1,13)",";  printf substr($6,1,13)","; printf substr($4,14,38)","; printf substr($6,14,38)"\\n";}}' > {output}
		"""

rule fastqtofasta:
	input:
		"data/fastq_pass/{readfile}.fastq.gz"
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
		demfile = rules.create_demfile.output
	output:
		"results/reference_free_barcoding.done"
	params:
		len = "600",
		depth = 10000
	singularity:
		"docker://reslp/minibarcoder:5e1dc3b"
	shell:
		"""
		rm -rf results/reffreebarcoding
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
		Rscript --vanilla -e "library(rmarkdown); rmarkdown::render('report.Rmd')"
		"""
