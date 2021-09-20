import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import json

wd = "/nfs/proj/REPO-TRIAL/montserrat"

### Number of patients interacting with each drug ###

print("Computing patient-drug statistics and plots")

pat_drug_matrix = pd.read_csv("%s/Data/Drugs/patients-drugs_table_-1_0_1_v5.csv" % wd, index_col=0)
pat_drug_matrix = pat_drug_matrix.replace({0:np.nan})
pat_drug_matrix_counts = pat_drug_matrix.count()
pat_drug_matrix_counts = pat_drug_matrix_counts.sort_values()

hist1 = pat_drug_matrix_counts.plot(kind="hist", bins=80, logy=True)
fig1 = hist1.get_figure()
fig1.savefig("%s/Data/Drugs/Plots/Distr_n_of_patients_drug_v5.pdf" % wd)
plt.show()

pat_drug_matrix_counts.to_csv("%s/Data/Drugs/drug_counts_v5.txt" % wd, sep="\t", header=True)


### Number of drugs interacting with each patient ###
pat_drug_matrix_drugcounts = pat_drug_matrix.count(axis=1).sort_values()

pat_drug_matrix_drugcounts

hist2 = pat_drug_matrix_drugcounts.plot(kind="hist", bins=80, 
	title="Amount of drugs/patient")
plt.axvline(x=pat_drug_matrix_drugcounts.mean(), c="r")
fig2 = hist2.get_figure()
fig2.savefig("%s/Data/Drugs/Plots/Distr_n_of_drugs_patient_v5.pdf" % wd)
plt.show()


pat_drug_matrix_drugcounts.to_csv("%s/Data/Drugs/pat_counts_v5.txt" % wd, sep="\t", header=True)

print("Finished!\n")
### Number of drugs interacting with each disease ###

print("Starting with drug-disease associations")
'''
with open("%s/Data/Drugs/drugs_v2.json" % wd, "r") as f:
	patients_drugs = json.load(f)

patients_count = {}

with open("%s/Patients_count.txt" % wd, "r") as fd:
	for line in fd:
		line = line.rstrip()
		dis, number = line.split("\t")
		patients_count[dis] = number

info = pd.read_csv("%s/Patients_information_table_v2.txt" % wd)

drugs_ids = pd.read_csv("%s/Mapping/drugs_mapping_lincs.txt" % wd, sep="\t")
drugnames = drugs_ids["Drug_lincs"].to_list()

d = {}

for n,patient in enumerate(patients_drugs):
	disease = info[info["Patient"] == patient]["Disease"].item()

	pre_mimickers = set(patients_drugs[patient]["mimickers"])
	mimickers = []

	for drug in pre_mimickers:
		if drug in drugnames:
			mimickers.append(drug)

	pre_reversers = set(patients_drugs[patient]["reversers"])
	reversers = []

	for drug in pre_reversers:
		if drug in drugnames:
			reversers.append(drug)
	if disease not in d.keys():
		d[disease] = {"mimickers": mimickers, "reversers": reversers}
	else:
		d[disease]["mimickers"].extend(mimickers)
		d[disease]["reversers"].extend(reversers)
	print("Patient %i of %i" % (n+1, len(patients_drugs)))

d2 = {}

for m, disease in enumerate(d):
	pat_number = int(patients_count[disease])
	count_mimickers = Counter(d[disease]["mimickers"])
	for key in count_mimickers:
		count_mimickers[key] /= pat_number
	if disease not in d2.keys():
		d2[disease] = {"mimickers": count_mimickers.most_common()}
	else:
		d2[disease]["mimickers"] = count_mimickers.most_common()
	count_reversers = Counter(d[disease]["reversers"])
	for key in count_reversers:
		count_reversers[key] /= pat_number
	d2[disease]["reversers"] = count_reversers.most_common()
	print("Disease %i of %i" % (m+1, len(d)))

with open("%s/Data/Drugs/drugs-diseases-associations_v2.json") as fw:
	json.dumps(d2, fw)

'''

info = pd.read_csv("%s/Patients_information_table_v2.txt" % wd)
diseases = info["Disease"].tolist()

d = {}


for disease in set(diseases):
	df = pat_drug_matrix[pat_drug_matrix.index.str.contains(disease, regex=False, case=False)]
	df_T = df.transpose()
	mimickers = []
	reversers = []
	for n, column in enumerate(df_T.columns):
		pre_mimickers = []
		pre_reversers = []
		mim = df_T[df_T[column] > 0].index.values
		for drug in mim:
			pre_mimickers.append(drug)
		rev = df_T[df_T[column] < 0].index.values
		for drug in rev:
			pre_reversers.append(drug)
		mimickers = mimickers + pre_mimickers
		reversers = reversers + pre_reversers
	m_counts = Counter(mimickers)
	for key in m_counts:
		m_counts[key] /= (n+1)/100
	r_counts = Counter(reversers)
	for key in r_counts:
		r_counts[key] /= (n+1)/100
	d[disease] = {"mimickers": m_counts.most_common(), "reversers": r_counts.most_common(), "patients": n+1}

with open("%s/Data/Drugs/drugs-disease-assoc_v5.json" % wd, "w") as fw:
	json.dump(d, fw)

##########################################################################

with open("%s/Data/Drugs/drugs-disease-assoc_v5.json" % wd, "r") as fw:
	d = json.load(fw)

info = pd.read_csv("%s/Patients_information_table_v2.txt" % wd)
diseases = info["Disease"].tolist()

inter_d = {}

diseases = list(set(diseases))

for disease1 in diseases:
	mimickers1 = [drug[0] for drug in d[disease1]["mimickers"][:10]]
	reversers1 = [drug[0] for drug in d[disease1]["reversers"][:10]]
	inter_d[disease1] = {}
	for disease2 in diseases:
		if disease1 == disease2:
			continue
		mimickers2 = [drug[0] for drug in d[disease2]["mimickers"][:10]]
		reversers2 = [drug[0] for drug in d[disease2]["reversers"][:10]]
		mim_mim = set(mimickers1).intersection(mimickers2)
		rev_rev = set(reversers1).intersection(reversers2)
		mim_rev = set(mimickers1).intersection(reversers2)
		rev_mim = set(reversers1).intersection(mimickers2)
		inter_d[disease1][disease2] = {"mim-mim":mim_mim, "rev_rev":rev_rev, "mim-rev": mim_rev, "rev_mim": rev_mim}

