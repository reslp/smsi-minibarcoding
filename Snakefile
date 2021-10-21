configfile: "data/config.yaml"
import pandas as pd
samples = pd.read_csv(config["data_file"], sep="\t")["name"].to_list()
sample_data =  pd.read_csv(config["data_file"], sep="\t")

RF, = glob_wildcards(config["fastq_files"]+"/{readfile}.fastq.gz")

def get_all_reports(wildcards):
	names = sample_data["name"].to_list()
	names = ["results/"+name+"/report.html" for name in names]
	return names

#def get_all_barcodes(wildcards):
#	print(wildcards)
#	demfile = sample_data.loc[sample_data["name"] == wildcards.n]["demfile"].tolist()
#	file = open(demfile, "r")
#	barcodes = []
#	for line in file:
#		barcodes.append(line.split(",")[0])
#	print(barcodes)
#	return barcodes

def get_all_barcodes(wildcards):
	barcodes = []
	samples = pd.read_csv(config["data_file"], sep="\t")["name"].to_list()
	for sample in samples:
		demfile = sample_data.loc[sample_data["name"] == sample]["demfile"].tolist()
		file = open(demfile[0], "r")
		for line in file:
			if os.path.exists("results/" + sample + "/reffreebarcoding/temp_alls_500/" + line.split(",")[0] + "_all.fa"):
				barcodes.append("results/" + sample + "/postprocessing/individual_barcodes/reoriented/"+line.split(",")[0]+"_all_barcodes_5P.fasta")
	#print(barcodes)
	return barcodes


rule all:
	input:
		get_all_reports,
		expand("results/{n}/postprocessing/{n}_blast_results.xml", n=samples),
		expand("results/{n}/postprocessing/reoriented_barcodes/{n}_all_barcodes_5P.fasta", n=samples),			
		get_all_barcodes
#		expand("results/{n}/postprocessing/individual_barcodes/reoriented/{bc}_all_barcodes_5P.fasta", n=samples, bc = get_all_barcodes)

rule barcoding:
	input:
		get_all_reports

rule postprocessing:
	input:
		expand("results/{n}/postprocessing/{n}_blast_results.xml", n=samples),
		expand("results/{n}/postprocessing/reoriented_barcodes/{n}_all_barcodes_5P.fasta", n=samples),			
		get_all_barcodes

include: "rules/functions.smk"
include: "rules/barcoding.smk"
include: "rules/report.smk"	
include: "rules/postbarcoding.smk"
