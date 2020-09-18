
# smsi-minibarcode

This workflow analyses is a snakemake/sigularity wrapper around the [miniBarcoder](https://github.com/asrivathsan/miniBarcoder) tool. It will take bascalled Minion seqeuence data and extract barcode sequences, count data and provide an overview report about the barcoding experiment.

I will run this on Sauron in this folder:

	/cl_tmp/reslph/projects/minibarcoding
	
	
## Manual steps:

**Data preparation:**

The demultiplex file needs to be created from Jackys info about the samples:

	$ tail -n +2 Primer_combination_MinION_Teil1.csv | awk -F ";" '{printf $1","; printf substr($4,1,13)",";  printf substr($6,1,13)","; printf substr($4,14,38)","; printf substr($6,14,38)"\n";}' > demfile.txt
	

MinION FASTQ files need to be converted to FASTA as well. Here an example for a single file:

	$ zcat FAO01594_pass_6519cb82_1.fastq.gz | paste - - - - | cut -f 1,2 | sed 's/^@/>/' > FAO01594_pass_6519cb82_1.fasta
	

**Analysis:**

These two files are needed to create first, reference-free consensus barcodes:

	$ singularity shell docker://reslp/minibarcoder:5e1dc3b
	

This is how I ran the barcoder inside the singularity container:

	$ python /software/miniBarcoder/miniBarcoder.py -f FAO01594_pass_6519cb82_1.fasta -d demfile.txt -o results/reffreebarcoding -D 10000 -t 8 -l 600
	

The main results file from this is called:

	all_barcodes.fa 

in the directory specified above.



From this file and the original sequence FASTA I had to create two files which are needed for creating the summary report:

	$ cat results/reffreebarcoding/all_barcodes.fa | grep ">" | awk -F ";" '{printf substr($1,2,length($1)-5)","; printf $2","; printf $3"\n";}' > barcode_summary.txt
	$ cat FAO01594_pass_6519cb82_1.fasta | grep ">" -c > total_seqs.txt
	



Now we can create a report a PDF report using a Docker container. The files report.Rmd, total_seqs.txt and barcode_summary.txt need to be in the wd.

	$ docker run --rm -it -v $(pwd):/data reslp/rmarkdown:3.6.3 Rscript -e "rmarkdown::render('./report.Rmd')"
	

This command should create a filde called report.pdf which gives a summary of the barcoding experiment.



## Automation:

I also created a snakemake pipeline to automate all the necessary steps. For it to work the folder with FASTQ files from Minion (fastq_pass) should be placed in data.

The snakemake command for this is:

	$ snakemake -p --use-singularity --singularity-args "-B $(pwd)/tmp:/usertmp -B $(pwd)/tmp/texlive:/usr/local/texlive -B $(pwd)/tmp/opt_texlive:/opt/texlive" --cores 8
	

There is also a submission script to do this on Sauron:

	$ gsub minibarcode.sge

