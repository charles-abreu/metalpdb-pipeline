
from glob import glob
import requests, os, sys, gzip
from Bio.PDB.MMCIFParser import MMCIFParser
from io import StringIO
import pymongo

MONGO_URL = "mongodb://localhost:27017/"
DB_NAME = "metalpdb2"
SITE_COLLECTION = "site"

PDB_STORE = "/bioinfo_data/WBIOCAT/pdb"

def unzip(file_name, outpu="string"): # string, bytes, file
  with gzip.open(file_name, 'rb') as f:
    file_content = f.read().decode("utf-8")
    return(StringIO(file_content))
  
def has_metals(model):
    
    metal_type_list = ["LI", "BE", "B", "NA", "MG", "AL", "SI" "K", "CA", "SC", "TI", "V", "CR", "MN",
        "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS" "RB", "SR", "Y", "ZR", "NB",
        "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "CS", "BA",
        "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER",
        "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
        "TL", "PB", "BI", "PO", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU",
        "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG","AS",
        "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC", "LV"]
    

    for atom in model.get_atoms():
        if atom.element.upper() in metal_type_list:
            return True
    
    return False
    

if __name__ == "__main__":

    client = pymongo.MongoClient(MONGO_URL)
    db = client[DB_NAME]
    collection = db[SITE_COLLECTION]
    
    dir_list = glob(f"{PDB_STORE}/*")
    dir_list.sort()

    already_done = []

    with open("dir_lsit.list") as in_file:
        already_done = [line.strip() for line in in_file.read().split("\n")]
    
    for d in dir_list:

        print(d)

        if d in already_done:
            continue

        pdb_list = glob(f"{d}/*.cif.gz")

        error_list = []
        with open("error_list.list") as in_file:
            error_list = [line.strip() for line in in_file.read().split("\n")]

        for i, pdb_file in enumerate(pdb_list):

            pdb_id = os.path.basename(pdb_file).split(".")[0]

            resultado = collection.find_one({"site_id": f"{pdb_id}_1"})
            
            if resultado:
                continue

            # check error list
            if pdb_file in error_list:
                continue

            # check file size
            tamanho_bytes = os.path.getsize(pdb_file)
            tamanho_mb = tamanho_bytes / (1024 * 1024)

            if tamanho_mb > 10:
                with open("size_list.list", "a") as out_file:
                        out_file.write(f"{pdb_file}\n")
                continue
            
            
            parser = MMCIFParser(QUIET= True)
            model = parser.get_structure(pdb_id, unzip(pdb_file))

            if has_metals(model):
                
                print(i, pdb_file)

                x = requests.get(f"http://localhost:50001/mfs/extract?pdb_id={pdb_id}&dssp=true&findgeo=true&insert=true")
                
                if x.status_code != 200:
                    print(f"{pdb_file} - ERROR")
                    with open("error_list.list", "a") as out_file:
                        out_file.write(f"{pdb_file}\n")
                
                    

        with open("dir_lsit.list", "a") as out_file:
            out_file.write(f"{d}\n")