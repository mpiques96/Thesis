import pandas as pd
import numpy as np

drugs_info = pd.read_csv("../Data/Drugs/drugs_names.txt", names=["drug"])
drugs_id = pd.read_csv("../Mapping/drug_id_name.txt", sep="\t", names=["drug_id", "drug_name"])
drugs_synonyms = pd.read_csv("../Mapping/drug_id_synonyms.txt", sep="\t", names=["drug_id", "drug_syn"])
drugs_synonyms = drugs_synonyms.dropna()

drugs = drugs_info["drug"].tolist()
ids_n = {}
ids_s = {}

for drug in drugs:
	id_n = drugs_id[drugs_id["drug_name"].str.lower() == drug]["drug_id"].tolist()
	if len(id_n) > 0:
		ids_n[drug] = id_n
	else:
		id_s = drugs_synonyms[drugs_synonyms["drug_syn"].str.contains(drug, regex=False, case=False)]["drug_id"].tolist()
		if len(id_s) > 0:
			ids_s[drug] = id_s

for drug in ids_s.keys():
	l = drugs_synonyms[drugs_synonyms["drug_syn"].str.contains(drug, regex=False, case=False)]
	ids = l["drug_id"].tolist()
	l = l["drug_syn"].tolist()
	for elem in l:
		names = elem.split(", ")
		if drug in names:
			ids_n[drug] = ids

df = pd.DataFrame.from_dict(ids_n, orient="index")
df.reset_index(level=0, inplace=True)
df = df.rename(columns={0: "Drugbank_id", "index": "Drugname"})

df.to_csv("../Data/Drugs/drugs_mapping.txt", sep="\t", index=False)