import numpy as np
import pandas as pd
import json

wd = "/nfs/proj/REPO-TRIAL/montserrat"


print("We start loading the data from LINCS")

with open("%s/Data/Drugs/drugs_v3.json" % wd, "r") as f:
	patients_drugs = json.load(f)

'''
namedrugs = []

for patient in patients_drugs:
	mimickers = patients_drugs[patient]["mimickers"]
	reversers = patients_drugs[patient]["reversers"]
	for drug in mimickers:
		namedrugs.append(str(drug))
	for drug in reversers:
		namedrugs.append(str(drug))

print("Drugs names extracted")

namedrugs_unique = list(set(namedrugs))
namedrugs_file = open("%s/Data/Drugs/drugs_names_list2.txt" % wd, "w")

for elem in namedrugs_unique:
	namedrugs_file.write(elem+"\n")

namedrugs_file.close()

print("Starting the drugnames mapping")

drugs_id = pd.read_csv("%s/Mapping/drug_id_name.txt" % wd, sep="\t", names=["drug_id", "drug_name"])
drugs_synonyms = pd.read_csv("%s/Mapping/drug_id_synonyms.txt" % wd, sep="\t", names=["drug_id", "drug_syn"])
drugs_synonyms = drugs_synonyms.dropna()

ids_s = {}

f = open("%s/Mapping/drugs_mapping_lincs2.txt" % wd, "w")

f.write("Drug_lincs\tDrugBankID")

for drug in namedrugs_unique:
	id_n = drugs_id[drugs_id["drug_name"].str.lower() == drug]["drug_id"].tolist()
	if len(id_n) > 0:
		f.write("%s\t%s\n" % (drug, id_n))
	else:
		l = drugs_synonyms[drugs_synonyms["drug_syn"].str.contains(drug, regex=False, case=False)]
		ids = l["drug_id"].tolist()
		l = l["drug_syn"].tolist()
		for elem in l:
			names = elem.split(", ")
			if drug in names:
				f.write("%s\t%s\n" % (drug, ids))

f.close()

print("Mapping done!")'''

drugs_ids = pd.read_csv("%s/Mapping/drugs_mapping_lincs2.txt" % wd, sep="\t")
drugnames = drugs_ids["Drug_lincs"].to_list()
drug_dict = {}

for idx, drug in enumerate(drugnames):
	drug_dict[drug] = idx

drug_matrix = np.empty((0,len(drugnames)), int)
patients_names = []

print("We begin with creating the patient-drug matrix")

for n, patient in enumerate(patients_drugs):
	patients_names.append(patient)
	row = np.zeros(len(drugnames),int)
	mimickers = set(patients_drugs[patient]["mimickers"])
	for drug in mimickers:
		try:
			drug_idx = drug_dict[str(drug)]
			row[drug_idx] = 1
		except:
			pass
	reversers = set(patients_drugs[patient]["reversers"])
	for drug in reversers:
		try:
			drug_idx = drug_dict[str(drug)]
			row[drug_idx] = -1
		except:
			pass
	drug_matrix = np.vstack((drug_matrix, row))
	print("Patient n°%i out of %i completed." % (n+1, len(list(patients_drugs))))

print("Patient-drug matrix completed!")

df = pd.DataFrame(data=drug_matrix, index=patients_names, columns=drugnames)
df.to_csv("/nfs/proj/REPO-TRIAL/montserrat/Data/Drugs/patients-drugs_table_-1_0_1_v5.csv")
print("Table finished!")
