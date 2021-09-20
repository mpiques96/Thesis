import pandas as pd
import numpy as np

df = pd.read_csv("../Results_all/Patient_patient_similarity_network.txt",sep="\t")

dis_dict = {"Asthma_central":"Asthma",
       "Asthma_peripheral":"Asthma",
       "Atopic_dermatitis_non_lesional":"Atopic_dermatitis",
       "Endometriosis_early_mild":"Endometriosis",
       "Endometriosis_early_severe":"Endometriosis",
       "Endometriosis_mid_mild":"Endometriosis",
       "Endometriosis_mid_severe":"Endometriosis",
       "Endometriosis_prol_mild":"Endometriosis",
       "Endometriosis_prol_severe":"Endometriosis",
       "Ependymoma_1":"Ependymoma",
       "Ependymoma_2":"Ependymoma",
       "Erythematotelangiectatic_rosacea":"Rosacea",
       "Human_Papillomavirus_Cervix_1":"Human_Papillomavirus_Cervix",
       "Human_Papillomavirus_Cervix_2":"Human_Papillomavirus_Cervix",
       "Obesity_fat":"Obesity",
       "Obesity_skin":"Obesity",
       "Papulopustular_rosacea":"Rosacea",
       "Phymatous_rosacea":"Rosacea",
       "Stage_1_Huntingtons_disease":"Huntingtons_disease",
       "Stage_2_Huntingtons_disease":"Huntingtons_disease",
       "Vitiligo_lesional":"Vitiligo",
       "Vitiligo_non_lesional":"Vitiligo",
       "Vitiligo_peri_lesional":"Vitiligo"}

dis1 = []
dis2 = []

print("Starting creating the disease columns")

for idx,row in df.iterrows():
	n1 = row[0].split("_GSE")[0]
	n2 = row[1].split("_GSE")[0]
	if n1 in dis_dict.keys():
		dis1.append(dis_dict[n1])
		if n2 in dis_dict.keys():
			dis2.append(dis_dict[n2])
		else:
			dis2.append(n2)
	elif n2 in dis_dict.keys():
		dis1.append(n1)
		dis2.append(dis_dict[n2])
	else:
		dis1.append(n1)
		dis2.append(n2)

# interaction = df["Interaction"].tolist()

# d = {"Disease1": dis1, "Disease2": dis2, "Interaction": interaction}

# df2 = pd.DataFrame(d)

# df2 = df2[df2['Disease1'] != df2['Disease2']]
# df2_count = df2.groupby(["Disease1","Disease2","Interaction"]).size().reset_index(name="edges")
# df2_count.to_csv("../Results_all/Disease_disease_count.txt",sep="\t",index=False)

# patients = pd.read_csv("../Results_all/Patients_information_table.txt")
# for key in dis_dict:
# 	patients.loc[patients["Disease"] == key,"Disease"] = dis_dict[key]
# patients_count = patients.groupby("Disease").size().reset_index(name="pat_number")
# patients_count.to_csv("../Results_all/Patients_count.txt",sep="\t",index=False)

# pat_a = []
# pat_b = []
# totals = []
# percs = []

# for idx,row in df2_count.iterrows():
# 	d1 = row[0]
# 	d2 = row[1]
# 	a = patients_count.loc[patients_count["Disease"] == d1, "pat_number"].values
# 	pat_a.append(a.item())
# 	b = patients_count.loc[patients_count["Disease"] == d2, "pat_number"].values
# 	pat_b.append(b.item())
# 	total = a * b
# 	perc = (row[3]/total)*100
# 	totals.append(total.item())
# 	percs.append(perc.item())

# df2_count["N°_patients_Disease1"] = pat_a
# df2_count["N°_patients_Disease2"] = pat_b
# df2_count["Total_possible_edges"] = totals
# df2_count["Percentages"] = percs


# df_all = pd.concat([df,df2],axis=1)
# df_all = df_all.iloc[:, :-1]

# df_all.to_csv("pat-pat-interaction.txt",index=False)



#####################################
#### NEW disease-disease network ####
#####################################

# Repetim el mateix fins línea 49

df_dis = df
df_dis.columns = ["Patient_1", "Patient_2", "Interaction"]
df_dis["Disease_1"] = dis1
df_dis["Disease_2"] = dis2
df_dis = df_dis[df_dis['Disease_1'] != df_dis['Disease_2']]

print("Disease columns added!")

diseases1 = set(df_dis["Disease_1"].tolist())
diseases2 = set(df_dis["Disease_2"].tolist())


df_counts = pd.DataFrame(columns = ["Disease_1","Pats_1","Disease_2","Pats_2","Interaction"])

print("Starting to count patients involved in each disease-disease association")
for idx,disease1 in enumerate(diseases1):
	for disease2 in diseases2:
		df_temp = df_dis[(df_dis["Disease_1"] == disease1) & (df_dis["Disease_2"] == disease2)]
		
		# Positive interaction
		df_temp_pos = df_temp[df_temp["Interaction"] == 1]
		pat1_pos = len(set(df_temp_pos["Patient_1"].tolist()))
		pat2_pos = len(set(df_temp_pos["Patient_2"].tolist()))
		if (pat1_pos > 0) and (pat2_pos > 0):
			df_counts = df_counts.append({
				"Disease_1": disease1,
				"Pats_1": pat1_pos,
				"Disease_2": disease2,
				"Pats_2": pat2_pos,
				"Interaction": 1
				}, ignore_index=True)

		# Negative interaction
		df_temp_neg = df_temp[df_temp["Interaction"] == -1]
		pat1_neg = len(set(df_temp_neg["Patient_1"].tolist()))
		pat2_neg = len(set(df_temp_neg["Patient_2"].tolist()))
		if (pat1_neg >0 and pat2_neg > 0):
			df_counts = df_counts.append({
				"Disease_1": disease1,
				"Pats_1": pat1_neg,
				"Disease_2": disease2,
				"Pats_2": pat2_neg,
				"Interaction": -1
				}, ignore_index=True)
	print("Disease n°%i out of %i" % (idx, len(diseases1)))

patients_count = pd.read_csv("../Results_all/Patients_count.txt", sep="\t")
total_pats1 = []
total_pats2 = []
for idx,row in df_counts.iterrows():
	d1 = row[0]
	d2 = row[2]
	pat1_num = patients_count.loc[patients_count["Disease"] == d1, "pat_number"].values
	pat2_num = patients_count.loc[patients_count["Disease"] == d2, "pat_number"].values
	total_pats1.append(pat1_num.item())
	total_pats2.append(pat2_num.item())

df_counts["Total_pats1"] = total_pats1
df_counts["Total_pats2"] = total_pats2

print("Number of patients involved in disease-disease pairs completed!")

print("Starting with the selection of at least a certain percentage")
percs = [10, 20, 30, 40, 50, 60]
for i in percs:
	d_d = pd.DataFrame(columns = ["Disease1","Pats1","Perc1","Disease2","Pats2","Perc2","Interaction"])
	for idx, row in df_counts.iterrows():
		perc1 = (row[1]/row[5])*100
		perc2 = (row[3]/row[6])*100
		if (perc1 >= i) and (perc2 >= i):
			d_d = d_d.append({
				"Disease1": row[0],
				"Perc1": perc1,
				"Pats1": row[5],
				"Disease2": row[2],
				"Perc2": perc2,
				"Pats2": row[6],
				"Interaction": row[4]
				}, ignore_index=True)
	# d_d2 = d_d.drop_duplicates(subset=["Disease1","Disease2"], keep=False)
	# d_d2.to_csv("../Results_all/D-D-network/Disease_disease_similarity_network_%s.txt" % i,sep="\t",index=False)
	d_d.to_csv("../Results_all/D-D-network/Disease_disease_similarity_network_%s.txt" % i,sep="\t",index=False)
	print("%i percentage finished!" % i)