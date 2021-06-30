rule get_summary_for_plotting:
	input:
		checkpoint = rules.reference_free_barcoding.output,
		fasta = rules.combine_fasta.output
	output:
		summary = "results/{name}/barcode_summary.txt",
		total_seqs = "results/{name}/total_seqs.txt"
	params:
		name = "{name}"
	shell:
		"""
		cat results/{params.name}/reffreebarcoding/all_barcodes.fa | grep ">" | awk -F ";" '{{printf substr($1,2,length($1)-5)","; printf $2","; printf $3"\\n";}}' > {output.summary}
		cat {input.fasta} | grep ">" -c > {output.total_seqs}
		"""
rule report:
	input:
		rules.get_summary_for_plotting.output.summary,
		rules.get_summary_for_plotting.output.total_seqs
	output:
		"results/{name}/report.html"
	params:
		wd = os.getcwd(),
		name = "{name}"
	singularity:
		"docker://reslp/rmarkdown:4.0.3"
	shell:
		"""
		#cp -rf /usr/local/texlive tmp/texlive
		#cp -rf /opt/texlive/ tmp/opt_texlive
		cp -f {params.wd}/data/report.Rmd {params.wd}/results/{params.name}/report.Rmd
		cd results/{params.name}
		Rscript --vanilla -e "library(rmarkdown); rmarkdown::render('report.Rmd')" $(cat fasta_combined.fasta | grep ">" | wc -l)
		"""
