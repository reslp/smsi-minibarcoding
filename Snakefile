configfile: "data/config.yaml"
import pandas as pd
samples = pd.read_csv(config["data_file"], sep="\t")["name"].to_list()

RF, = glob_wildcards(config["fastq_files"]+"/{readfile}.fastq.gz")

def get_all_reports(wildcards):
	names = sample_data["name"].to_list()
	names = ["results/"+name+"/report.html" for name in names]
	return names

rule all:
	input:
		get_all_reports,
		expand("results/{n}/postprocessing/{n}_blast_results.xml", n=samples)
			
include: "rules/functions.smk"
include: "rules/barcoding.smk"
include: "rules/report.smk"	
include: "rules/postbarcoding.smk"
