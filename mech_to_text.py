import os
import sys
import re

import pandas as pd 
import numpy as np

from GECKO_to_SMILES import *

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Functions

def cell_to_SMILES(df_cell):
    
    GECKO_smiles = []
    
    for item in df_cell:
        
        try:
            nome = item.split()[1]
#            gecko_id = item.split()[0]
            GECKO_smiles.append(nome)
        except:
            break
            
    return GECKO_smiles


def get_smiles_list(filename):
    
    df_smiles = pd.read_csv(filename, header = None)
    species_list = df_smiles.apply(lambda x: cell_to_SMILES(x))[0]
    
    return species_list


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MAIN TASK

sim_dir = os.getcwd()

dictionary_files = []
mechanism_files = []

for r, d, f in os.walk(sim_dir):
    for file in f:
        if '.dict' in file:
            dictionary_files.append(os.path.join(r, file))
        elif '.mech' in file:
            mechanism_files.append(os.path.join(r, file))

# 1) Generating a list of all GECKO-A SMILES to us in RDKit: GECKO_smiles 
GECKO_smiles = []

# Open all "*.dict" files in working directory
#for file in os.listdir(sim_dir):
for file in dictionary_files:
#    if file.endswith(".dict"):
        # Read GECKO smiles per each file
        file_smiles_list = get_smiles_list(file)
        
        # Select only unique smiles using set() and list()
        GECKO_smiles = list(set(GECKO_smiles) | set(file_smiles_list))

print("Total number of GECKO species: ", len(GECKO_smiles))
#pd.Series(GECKO_smiles).to_csv("atmospheric_species.csv", index = False, header=False)


SMILESlist_init = GECKO_smiles

print('\n Remove fictional species : \n LOSS CARBON \n --mm###')

SMILESlist, GECKO_screened = translate_smiles(SMILESlist_init)
#pd.Series(SMILESlist).to_csv("SMILES_species.csv", index = False, header=False)


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

MFPlist = []

fileA = open('ERROR_RDKit.txt', 'w')
fileB = open('ERROR_GECKO.txt', 'w')

fileSMI = open('prova.smi', 'w')

for ind, specie in enumerate(SMILESlist):

    smile = specie

    try:
#        bi = {}
        m1 = Chem.MolFromSmiles(smile)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, radius = 2, nBits = 1024)
        array = np.zeros((0,), dtype = np.int8)
        DataStructs.ConvertToNumpyArray(fp1, array)

        MFPlist.append(smile)

        fileSMI.write("{}\n".format(smile))        

    except:
        fileA.write("{}\n".format(specie))
        fileB.write("{}\n".format(GECKO_screened[ind]))

fileA.close()
fileB.close()
fileSMI.close()

print('\n')
print('\n')
print('\n')
print('\n')

print('TOTAL SMILESlist: ', len(SMILESlist))
print('TOTAL MFPlist: ', len(MFPlist))

print('ERROR smiles: ', len(SMILESlist) - len(MFPlist))