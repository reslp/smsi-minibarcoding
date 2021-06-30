import pandas as pd

sample_data = pd.read_csv(config["data_file"], sep ="\t").set_index("name", drop=False)

def get_reads_dir(wildcards):
	return sample_data.loc[wildcards.name, ["reads_dir"]].astype(str)

def get_demfile(wildcards):
	return sample_data.loc[wildcards.name, ["demfile"]].to_list()

def get_all_reports(wildcards):
	names = sample_data["name"].to_list()
	names = ["results/"+name+"report.html" for name in names]
	return names
			
