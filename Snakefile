configfile: "data/config.yaml"

RF, = glob_wildcards(config["fastq_files"]+"/{readfile}.fastq.gz")

def get_all_reports(wildcards):
	names = sample_data["name"].to_list()
	names = ["results/"+name+"/report.html" for name in names]
	return names

rule all:
	input:
		get_all_reports
	
include: "rules/functions.smk"
include: "rules/barcoding.smk"
include: "rules/report.smk"	
