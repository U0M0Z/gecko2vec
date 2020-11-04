
import sys
import os
import re

def translate_smiles(SMILESlist_init):

    SMILESlist = []
    SMILESlist_corrected = []
    GECKO_list_screened = []

    # (1) +++ main task +++: processing list of input SMILES from GECKO-A output and convert to valid SMILES +++++
    for ind,molecule in enumerate(SMILESlist_init):
        
        st = str(molecule)

        #Ignore some rows indicating issues like #lostcarbon in smiles file and conc file;    
        if not( re.match("#lostcarbon\n", st) or re.match("#mm", st)):  #then not an exception, so apply string processing:
            GECKO_list_screened.append(st)
            #remove GECKO-A notation of 'Cd' for a carbon with a double bond:
            st = st.replace('Cd','C')
            #safely convert certain specific GECKO-A component string notations (like CO, CO2) to SMILES:
            st = st.replace('CO2','O=C=O')
#            st = st.replace('CO2\n','O=C=O\n')
            st = st.replace('CO\n','[C-]#[O+]\n')
            st = st.replace('ONO2','O[N+]([O-])=O')
            st = st.replace('NO2','[N+](=O)[O-]')
            #replace GECKO-A 'CO' for aldehyde or ketone 'CHO' (or formaldehyde) within molecule by appropriate carbonyl group:
            st = st.replace('CO','C(=O)')
            st = st.replace('CH2O\n','C(=O)\n') #here a whole molecule substring as SMILES
            st = st.replace('CHO','C(=O)')
            #string beginning with - ether or peroxide notation:
            if st[0:1] == '-':
                st = st[1:]  #remove first character as it is a dash (-), which is not a valid SMILES notation as starting character;
            st = st.replace('(-O-)','(O)')
            st = st.replace('(-O1-)','(O1)')
            st = st.replace('(-O2-)','(O2)')
            st = st.replace('(-O3-)','(O3)')
            #radical species and peroxides
            st = st.replace('O.','[O]')
            st = st.replace('OO.','[O]O')
            st = st.replace('OOH', 'OO')
            st = st.replace('CO(OO.)', 'C(=O)[O]O')
            #remove explicit hydrogens stated:
            st = st.replace('H3','')
            st = st.replace('H2','')
            st = st.replace('H1','')
            st = st.replace('H','')
            st = st.replace('C.', '[C]')
            st = st.replace('-O1-', 'O1')

#           Special molecules
            st = st.replace('SO2', '[S](=O)=O')
            st = st.replace('sulf', 'OS(=O)(=O)O')
            st = st.replace('H2O2', 'O(O[H])[H]')
            st = st.replace('HNO3', '[N+](=O)([O-])O[H]')
            st = st.replace('NO3', '[N+](=O)([O-])O')
            st = st.replace('O3', '[O-][O+]=O')
            st = st.replace('HO2', 'O([H])[O]')
            st = st.replace('N2O5', '[N+](O[N+](=O)[O-])([O-])=O')

            st = st.replace('CH4', 'C')
            st = st.replace('O2', 'O=O')
            st = st.replace('NO4', '[N+](=O)([O-])OO')
            st = st.replace('HNO4', '[N+](=O)([O-])OO[H]')
            st = st.replace('C4', 'CCCC')


            #safe in modified SMILES list and input concentrations list:
            SMILESlist.append(st.strip())        #strip() removes leading and trailing whitespace, tabs, etc.
    
    for ind,molecule in enumerate(SMILESlist):
            
            st = str(molecule)

            #Ignore some rows indicating issues like #lostcarbon in smiles file and conc file;    
            if not( re.match("#lostcarbon\n", st) or re.match("#mm", st)):
                st = st.replace('--','-')
                SMILESlist_corrected.append(st.strip())


    return SMILESlist_corrected, GECKO_list_screened
#++++++++++++++++++++++++++++++++++++++++++++++++


# >>>> (2) >>>>  for tests, check with Indigo toolkit whether SMILES is valid:      
# this code block uses the epam Indigo API (it can be commentend out if Indigo is not installed);
#ii = len(SMILESlist)
#is_invalidSMILES = [True]*ii

#------------------------------------------------------------------------------------------