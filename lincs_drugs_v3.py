import requests
import json
from time import sleep

###### STEP 1 ######
# Convert gene names to UUIDs using the metadata API #

METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"


with open("/nfs/proj/REPO-TRIAL/montserrat/Data/Top_500_sDEGs_names.json") as f:
    list_patients = json.load(f)
l = list(list_patients)#[n-1:]

fp = open("/nfs/proj/REPO-TRIAL/montserrat/Data/Drugs/drugs_v3.json", "w")

n = 0
sleep_time = 2
for patient in l:
    input_gene_set = {
        "up_genes": list_patients[patient]["up"][:200],
        "down_genes": list_patients[patient]["down"][-200:]
    }
    all_genes = input_gene_set["up_genes"] + input_gene_set["down_genes"]

    payload = {
        "filter": {
            "where": {
                "meta.symbol": {
                    "inq": all_genes
                }
            },
            "fields": ["id", "meta.symbol"]
        }
    }
    for x in range(4):
        try:
            res = requests.post(METADATA_API + "entities/find", json=payload)
            str_error = None
        except Exception as str_error:
            pass
        if str_error:
            sleep(sleep_time)
            sleep_time *= 2
        else:
            break

    entities = res.json()

    for_enrichment = {
        "up_entities": [],
        "down_entities": []
    }

    for e in entities:
        symbol = e["meta"]["symbol"]
        if symbol in input_gene_set["up_genes"]:
            for_enrichment["up_entities"].append(e["id"])
        elif symbol in input_gene_set["down_genes"]:
            for_enrichment["down_entities"].append(e["id"])

    # print(json.dumps(for_enrichment, indent=2))


    ###### STEP 2 ######
    # Perform signature search #

    query = {
        **for_enrichment,
        "limit": 1000,
        "database": "l1000_cp"
    }

    for x in range(4):
        try:
            res = requests.post(DATA_API + "enrich/ranktwosided", json=query)
        except Exception as str_error:
            pass
        if str_error:
            sleep(sleep_time)
            sleep_time *= 2
        else:
            break

    results = res.json()

    # Optional, multiply z-down and direction-down with -1
    for i in results["results"]:
        i["z-down"] = -i["z-down"]
        i["direction-down"] = -i["direction-down"]

    # print(json.dumps(results, indent=2))


    ###### STEP 3 ######
    # Resolve signature UUIDs using the metadata API #

    sigids = {i["uuid"]: i for i in results["results"]}

    payload = {
        "filter": {
            "where": {
                "id": {
                    "inq": list(sigids.keys())
                }
            }
        }
    }
    for x in range(4):
        try:
            res = requests.post(METADATA_API + "signatures/find", json=payload)
        except Exception as str_error:
            pass
        if str_error:
            sleep(sleep_time)
            sleep_time *= 2
        else:
            break

    signatures = res.json()

    ## Merge the scores and the metadata
    for sig in signatures:
        uid = sig["id"]
        scores = sigids[uid]
        scores.pop("uuid")
        sig["scores"] = scores

    # print(json.dumps(signatures, indent=2))

    ##### STEP 4 #####
    # Save results for each patient #

    mimickers = []
    reversers = []

    for signature in signatures:
        z_sum = signature["scores"]["z-sum"]
        drug_name = signature["meta"]["pert_name"]
        if z_sum > 5:
            mimickers.append(drug_name)
        elif z_sum < -5:
            reversers.append(drug_name)

    drugs_results = {"mimickers": mimickers, "reversers": reversers}

    if n == 0:
        fp.write('{"%s": %s' % (patient, json.dumps(drugs_results)))
    else:
        fp.write(',"%s": %s' % (patient, json.dumps(drugs_results)))
    # plt.hist(z_scores, bins=20, range=(4, 20))
    # plt.title(patient)
    # plt.savefig("../Data/Drugs/Histograms/%s.png" % patient)
    # plt.close()

    n += 1
    print("Patient n°%i out of %i" % (n, len(list_patients)))


fp.write("}")
fp.close()