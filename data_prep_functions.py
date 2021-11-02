import re
from interpro_scraping import interpro_scraping_pandas
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import scipy.stats

def percent_containing(string_to_parse, length, search_string):
    return sum(1 for _ in re.finditer(r'\b%s\b' % re.escape(search_string), str(string_to_parse))) / length


def regularize_aa(aa_dict):
    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    keys = aa_dict.keys()
    if len(keys) != 20:
        for aa in aa_list:
            if aa not in keys:
                aa_dict[aa] = 0
    return aa_dict


def alphabetize_aa(aa_dict):
    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    aa_list.sort()
    return [aa_dict[aa] for aa in aa_list]


def pandas_column_list_expansion(data, column_of_interest):
    if column_of_interest == 'InterPro':
        quick_split = [y for x in data['InterPro'] for y in x]

    elif column_of_interest == 'Gene ontology IDs':
        quick_split = [y for x in data['Gene ontology IDs'] for y in str(x).split(';')]

    else:
        raise ValueError('That column has not been specified.'
                         ' Edit function or use InterPro or Gene ontology IDs')
    
    unique_list = []
    for i in quick_split:  # find every unique value
        if i not in unique_list:
            unique_list.append(i)

    if 'na' in unique_list:  # get rid of a na that appears
        unique_list.remove('na')

    for i in unique_list: # create a column for each unique value
        data[i] = 0  # put in the cleaned up dataframe

    # create booleans for unique values
    n_rows = data.shape[0]
    counter = 0
    for x in data[column_of_interest]:
        if column_of_interest == 'InterPro':
            value_list = x

        elif column_of_interest == 'Gene ontology IDs':
            value_list = str(x).split('; ')

        else:
            raise ValueError('That column has not been specified.'
                             ' Edit function or use InterPro or Gene ontology IDs')

        if value_list == ['']:
            counter += 1
            continue
        else:
            for value in value_list:
                data.loc[counter, value] = 1
        counter += 1

    return data

def clean_up_data(raw_data):
    raw_data = interpro_scraping_pandas(raw_data, 'Entry', 'InterPro')  # collect interpro numbers
    raw_data = raw_data.fillna(0) # fill nans

    # to be changed ... as of now if a value isnt 0 its 1 
    raw_data['DNA binding'] = raw_data['DNA binding'].apply(dna_binding_length)
    raw_data['dna_binding_length'] = raw_data['DNA binding'] / raw_data['Length'] 
    raw_data['calcium_binding_zones'] = raw_data['Calcium binding'].apply(calcium_zones)
    raw_data['calcium_binding_length'] = raw_data['Calcium binding'].apply(calcium_length) / raw_data['Length']
    raw_data['nucleotide_binding_zones'] = raw_data['Nucleotide binding'].apply(nucleotide_zones) 
    raw_data['nucleotide_binding_length'] = raw_data['Nucleotide binding'].apply(nucleotide_length) / raw_data['Length']
    

    raw_data['percent_glycosylated'] = [percent_containing(raw_data.Glycosylation[i], raw_data.Length[i], 'CARBOHYD') for i in range(len(raw_data.Length))]

    cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass', 'InterPro', 'dna_binding_length',
                           'Metal binding', 'calcium_binding_zones', 'calcium_binding_length', 'Nucleotide binding', "nucleotide_binding_zones",
                           "nucleotide_binding_length", 'Gene ontology IDs', 'percent_glycosylated']]

    cleaned_data = pandas_column_list_expansion(cleaned_data, 'InterPro')
    cleaned_data = pandas_column_list_expansion(cleaned_data, 'Gene ontology IDs')
    cleaned_data = cleaned_data.fillna(0)
    cleaned_data = metal_binding_analysis(cleaned_data)
    cleaned_data = nucleotide_types(cleaned_data)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        aromat = analyzed_seq.aromaticity() #
        instab = analyzed_seq.instability_index() # float
        flex = analyzed_seq.flexibility() # returns a list
        iso = analyzed_seq.isoelectric_point() 
        secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        secStruct_disorder = 1 - sum(secStruct)
        flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder]])
        

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight', 'aromaticity', 'instability_index', 'flexibility_mean', 'flexibility_std', 'flexibility_var',
                 'flexibility_max', 'flexibility_min', 'flexibility_median', 'isoelectric_point', 'secondary_structure_fraction_helix',
                 'secondary_structure_fraction_turn', 'secondary_structure_fraction_sheet', 'secondary_structure_fraction_disordered']
        
    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')
    cleaned_data = cleaned_data.fillna(0)
    return cleaned_data

def grab_max(df1, df2, column):
    if df1[column].max() >= df2[column].max():
        return df1[column].max()
    
    else:
        return df2[column].max()

def normalize_mass_length(df1, df2):

    max_length = grab_max(df1, df2, 'Length')
    max_mass = grab_max(df1, df2, 'Mass')
    max_mw = grab_max(df1, df2, 'molecular_weight')

    df1['length'] = df1['Length'] / max_length
    df2['length'] = df2['Length'] / max_length

    df1['mass'] = df1['Mass'] / max_mass
    df2['mass'] = df2['Mass'] / max_mass

    df1['molecular_weight'] = df1['molecular_weight'] / max_mw
    df2['molecular_weight'] = df2['molecular_weight'] / max_mw

    return df1, df2


def metal_binding_analysis(data):

    reg_list = [r"note\=\"(copper)", r"note\=\"(calcium)", r"note\=\"(zinc)", r"note\=\"(iron)", r"note\=\"(manganese)", r"note\=\"(nickel)", r"note\=\"(potassium)", r"note\=\"(sodium)", r"note\=\"(cobalt)"]
    metals = ['copper', 'calcium', 'zinc', 'iron', 'manganese', 'nickel', 'potassium', 'sodium', "cobalt"]


    for metal in metals:
        data['metal_binding_'+metal] = 0

    length = data['Metal binding'].shape[0]
    for i in range(length):
        current_string = data['Metal binding'][i]
        if current_string != 0:
            current_string = str(current_string)

            current_string = current_string.lower()
            for j in range(len(reg_list)):
                
                metal_found = re.findall(reg_list[j], current_string)

                data['metal_binding_'+metals[j]][i] = len(metal_found)

    data = data.drop(columns=['Metal binding'])
    return data    

def dna_binding_length(string_to_parse):
    if string_to_parse != 0:
        string_to_parse = str(string_to_parse)
        regex_list = re.findall(r"\d+", string_to_parse)
        if len(regex_list) == 2:
            return int(regex_list[1]) - int(regex_list[0])
    else: 
        return 0

def calcium_zones(string_to_parse):
    if string_to_parse != 0:
        string_to_parse = str(string_to_parse)
        return len(re.findall(r"(?:CA_BIND )(\d+)..(\d+)", str(string_to_parse)))
    else: 
        return 0

def calcium_length(string_to_parse):
    if string_to_parse != 0:
        total = 0 
        string_to_parse = str(string_to_parse)
        regex_list = re.findall(r"(?:CA_BIND )(\d+)..(\d+)", str(string_to_parse))
        for j , k in regex_list:
            total += (int(k)-int(j)) 
        return total

    else: 
        return 0

def nucleotide_zones(string_to_parse):
    if string_to_parse != 0:
        string_to_parse = str(string_to_parse)
        try:
            return len(re.findall(r"(?:NP_BIND )(\d+)..(\d+)", str(string_to_parse)))
        except:
            return 0
    else: 
        return 0

def nucleotide_length(string_to_parse):
    if string_to_parse != 0:
        total  = 0 
        string_to_parse = str(string_to_parse)
        regex_list = re.findall(r"(?:NP_BIND )(\d+)..(\d+)", str(string_to_parse))
        for j , k in regex_list:
            total += (int(k)-int(j)) 
        return total

    else: 
        return 0

def nucleotide_types(data):
    df = data.copy()

    unique_list = []
    update_list = []

    for i in range(len(df['Nucleotide binding'])):  # collect unique keys
        string_to_parse = df['Nucleotide binding'][i]
        if string_to_parse != 0:
            string_to_parse = str(string_to_parse)
            regex_list_type = re.findall(r"note\=\"(.*?)\"", str(string_to_parse))
            update_list.append(i)
            for j in regex_list_type:
                if j not in unique_list:
                    unique_list.append(j)

    # go back and add columns
    for i in unique_list:
        df['nucleotide_binding_'+i] = 0


    # go through the ones that were found to have hits
    for i in update_list:
        string_to_parse = df['Nucleotide binding'][i]
        string_to_parse = str(string_to_parse)
        regex_list_type = re.findall(r"note\=\"(.*?)\"", str(string_to_parse))

        for j in unique_list:
            df['nucleotide_binding_'+j][i] = regex_list_type.count(j) 
    
    df = df.drop(columns=['Nucleotide binding'])
    return df

def clean_up_data_no_interpro_go(raw_data):
    
    raw_data = raw_data.fillna(0) # fill nans

    # to be changed ... as of now if a value isnt 0 its 1 
    raw_data['DNA binding'] = raw_data['DNA binding'].apply(dna_binding_length)
    raw_data['dna_binding_length'] = raw_data['DNA binding'] / raw_data['Length'] 
    raw_data['calcium_binding_zones'] = raw_data['Calcium binding'].apply(calcium_zones)
    raw_data['calcium_binding_length'] = raw_data['Calcium binding'].apply(calcium_length) / raw_data['Length']
    raw_data['nucleotide_binding_zones'] = raw_data['Nucleotide binding'].apply(nucleotide_zones) 
    raw_data['nucleotide_binding_length'] = raw_data['Nucleotide binding'].apply(nucleotide_length) / raw_data['Length']
    

    raw_data['percent_glycosylated'] = [percent_containing(raw_data.Glycosylation[i], raw_data.Length[i], 'CARBOHYD') for i in range(len(raw_data.Length))]

    cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass', 'dna_binding_length',
                           'Metal binding', 'calcium_binding_zones', 'calcium_binding_length', 'Nucleotide binding', "nucleotide_binding_zones",
                           "nucleotide_binding_length", 'percent_glycosylated']]

    cleaned_data = cleaned_data.fillna(0)
    cleaned_data = metal_binding_analysis(cleaned_data)
    cleaned_data = nucleotide_types(cleaned_data)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        aromat = analyzed_seq.aromaticity() #
        instab = analyzed_seq.instability_index() # float
        flex = analyzed_seq.flexibility() # returns a list
        iso = analyzed_seq.isoelectric_point() 
        secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        secStruct_disorder = 1 - sum(secStruct)
        flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder]])
        

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight', 'aromaticity', 'instability_index', 'flexibility_mean', 'flexibility_std', 'flexibility_var',
                 'flexibility_max', 'flexibility_min', 'flexibility_median', 'isoelectric_point', 'secondary_structure_fraction_helix',
                 'secondary_structure_fraction_turn', 'secondary_structure_fraction_sheet', 'secondary_structure_fraction_disordered']
        
    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')

    return cleaned_data


def clean_up_data_biopy(raw_data):
    # least data only biopy
    raw_data = raw_data.fillna(0) # fill nans
    
    cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass']]

    
    cleaned_data = cleaned_data.fillna(0)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        aromat = analyzed_seq.aromaticity() #
        instab = analyzed_seq.instability_index() # float
        flex = analyzed_seq.flexibility() # returns a list
        iso = analyzed_seq.isoelectric_point() 
        secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        secStruct_disorder = 1 - sum(secStruct)
        flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        gravy = analyzed_seq.gravy()
        temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder, gravy]])
        

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight', 'aromaticity', 'instability_index', 'flexibility_mean', 'flexibility_std', 'flexibility_var',
                 'flexibility_max', 'flexibility_min', 'flexibility_median', 'isoelectric_point', 'secondary_structure_fraction_helix',
                 'secondary_structure_fraction_turn', 'secondary_structure_fraction_sheet', 'secondary_structure_fraction_disordered', 'gravy']
        
    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')
    cleaned_data = cleaned_data.fillna(0)
    return cleaned_data

def clean_up_data_biopy_no_ss_flex(raw_data):
    raw_data = raw_data.fillna(0) # fill nans
    
    cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass']]

    
    cleaned_data = cleaned_data.fillna(0)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        aromat = analyzed_seq.aromaticity() #
        instab = analyzed_seq.instability_index() # float
        flex = analyzed_seq.flexibility() # returns a list
        iso = analyzed_seq.isoelectric_point() 
        #secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        #secStruct_disorder = 1 - sum(secStruct)
        #flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        #temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder]])
        temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, iso]])

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight', 'aromaticity', 'instability_index',
                  'isoelectric_point']

    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')
    cleaned_data = cleaned_data.fillna(0)
    return cleaned_data


def clean_up_data_aa(raw_data):
    raw_data = raw_data.fillna(0) # fill nans
    
    cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass']]

    
    cleaned_data = cleaned_data.fillna(0)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        #aromat = analyzed_seq.aromaticity() #
        #instab = analyzed_seq.instability_index() # float
        #flex = analyzed_seq.flexibility() # returns a list
        #iso = analyzed_seq.isoelectric_point() 
        #secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        #secStruct_disorder = 1 - sum(secStruct)
        #flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        #temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder]])
        temp_df = pd.DataFrame([[seq, *aa_values, mw]])

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight']
                  
    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')
    cleaned_data = cleaned_data.fillna(0)
    return cleaned_data


### TO MERGE ####
def zf_length(string_to_parse):
    if string_to_parse != 0:
        total = 0 
        
        regex_list = re.findall(r"(?:ZN_FING )(\d+)..(\d+)", str(string_to_parse))
        
        for j , k in regex_list:
            total += (int(k)-int(j)) 
        return total

    else: 
        return 0

def zf_zones(string_to_parse):
    if string_to_parse != 0:
        
        return len(re.findall(r"(?:ZN_FING )(\d+)..(\d+)", str(string_to_parse)))
    else: 
        return 0

def regex_length(string_to_parse, regex_form):
    if string_to_parse != 0:
        total = 0 
        
        regex_list = re.findall(regex_form, str(string_to_parse))
        
        for j , k in regex_list:
            total += (int(k)-int(j)) 
        return total

    else: 
        return 0

def regex_zones(string_to_parse, regex_form):
    if string_to_parse != 0:
        return len(re.findall(regex_form, str(string_to_parse)))
    else: 
        return 0

def act_site_analysis(data):

    reg_list = [r"note\=\"(proton donor)", r"note\=\"(proton acceptor)", r"note\=\"(proton donor/acceptor)", r"note\=\"(charge relay system)"]
    metals = ['proton_donor', 'proton_acceptor', 'proton_donor_acceptor', 'charge_relay']

    #r"note\=\"(Proton donor/acceptor)", 'proton_donor_acceptor',
    for metal in metals:
        data['active_site_'+metal] = 0

    length = data['Active site'].shape[0]
    for i in range(length):
        current_string = data['Active site'][i]
        # print(current_string)
        if current_string != 0:
            current_string = str(current_string)
            # print(current_string)
            current_string = current_string.lower()
            #print(current_string)
            for j in range(len(reg_list)):
                
                metal_found = re.findall(reg_list[j], current_string)
                
                data['active_site_'+metals[j]][i] = len(metal_found)

    data = data.drop(columns=['Active site'])
    return data


def cofactor_analysis(data):

    reg_list = [r"name\=(fmn+)", r"name\=(ca+)", r"name\=(cu+)", r"name\=(mg+)", r"name\=(zn+)"]
    metals = ['fmn', 'ca', 'cu', 'mg', 'zn']

    #r"note\=\"(Proton donor/acceptor)", 'proton_donor_acceptor',
    for metal in metals:
        data['cofactor_'+metal] = 0

    length = data['Cofactor'].shape[0]
    for i in range(length):
        current_string = data['Cofactor'][i]
        
        if current_string != 0:
            current_string = str(current_string)
            current_string = current_string.lower()
            
            for j in range(len(reg_list)):
                
                metal_found = re.findall(reg_list[j], current_string)
                
                data['cofactor_'+metals[j]][i] = len(metal_found)

    data = data.drop(columns=['Cofactor'])
    return data


def clean_up_data_v2(data):
    raw_data = data.fillna(0) # fill nans

    # to be changed ... as of now if a value isnt 0 its 1 
    # many of these could be rewritten using new functions but they were left using the original function
    raw_data['DNA binding'] = raw_data['DNA binding'].apply(dna_binding_length)
    raw_data['dna_binding_length'] = raw_data['DNA binding'] / raw_data['Length'] 

    raw_data['calcium_binding_zones'] = raw_data['Calcium binding'].apply(calcium_zones)
    raw_data['calcium_binding_length'] = raw_data['Calcium binding'].apply(calcium_length) / raw_data['Length']
    raw_data = raw_data.drop(columns=['Calcium binding'])

    raw_data['nucleotide_binding_zones'] = raw_data['Nucleotide binding'].apply(nucleotide_zones) 
    raw_data['nucleotide_binding_length'] = raw_data['Nucleotide binding'].apply(nucleotide_length) / raw_data['Length']

    raw_data['zinc_finger_zones'] = raw_data['Zinc finger'].apply(regex_zones, regex_form=r"(?:ZN_FING )(\d+)..(\d+)")
    raw_data['zinc_finger_length'] = raw_data['Zinc finger'].apply(regex_length, regex_form=r"(?:ZN_FING )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Zinc finger'])

    raw_data['intramembrane_zones'] = raw_data['Intramembrane'].apply(regex_zones, regex_form=r"(?:INTRAMEM )(\d+)..(\d+)")
    raw_data['intramembrane_length'] = raw_data['Intramembrane'].apply(regex_length, regex_form=r"(?:INTRAMEM )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Intramembrane'])

    raw_data['transmembrane_zones'] = raw_data['Transmembrane'].apply(regex_zones, regex_form=r"(?:TRANSMEM )(\d+)..(\d+)")
    raw_data['transmembrane_length'] = raw_data['Transmembrane'].apply(regex_length, regex_form=r"(?:TRANSMEM )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Transmembrane'])

    raw_data['coiled_coil_zones'] = raw_data['Coiled coil'].apply(regex_zones, regex_form=r"(?:COILED )(\d+)..(\d+)")
    raw_data['coiled_coil_length'] = raw_data['Coiled coil'].apply(regex_length, regex_form=r"(?:COILED )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Coiled coil'])

    raw_data['lipidation_zones'] = raw_data['Lipidation'].apply(regex_zones, regex_form=r"(?:LIPID )(\d+)..(\d+)")
    raw_data = raw_data.drop(columns=['Lipidation'])

    raw_data['binding_site_zones'] = raw_data['Binding site'].apply(regex_zones, regex_form=r"(?:BINDING )(\d+)")
    raw_data = raw_data.drop(columns=['Binding site'])

    raw_data['topo_dom_extracellular_zones'] = raw_data['Topological domain'].apply(regex_zones, regex_form=r"(?:TOPO_DOM )(\d+)..(\d+)(?:\;  \/note\=\"Extracellular)") 
    raw_data['topo_dom_extracellular_length'] = raw_data['Topological domain'].apply(regex_zones, regex_form=r"(?:TOPO_DOM )(\d+)..(\d+)(?:\;  \/note\=\"Extracellular)") / raw_data['Length']
    raw_data['topo_dom_cytoplasmic_zones'] = raw_data['Topological domain'].apply(regex_zones, regex_form=r"(?:TOPO_DOM )(\d+)..(\d+)(?:\;  \/note\=\"Cytoplasmic)")
    raw_data['topo_dom_cytoplasmic_length'] = raw_data['Topological domain'].apply(regex_zones, regex_form=r"(?:TOPO_DOM )(\d+)..(\d+)(?:\;  \/note\=\"Cytoplasmic)") / raw_data['Length']
    raw_data['topo_dom_total_zones'] = raw_data['topo_dom_extracellular_zones'] + raw_data['topo_dom_cytoplasmic_zones']
    raw_data['topo_dom_total_length'] = raw_data['topo_dom_extracellular_length'] + raw_data['topo_dom_cytoplasmic_length']
    raw_data = raw_data.drop(columns=['Topological domain'])

    raw_data['mod_res_zones'] = raw_data['Modified residue'].apply(regex_zones, regex_form=r"(?:MOD_RES )(\d+)")
    raw_data = raw_data.drop(columns=['Modified residue'])

    raw_data['disulfide_zones'] = raw_data['Disulfide bond'].apply(regex_zones, regex_form=r"(?:DISULFID )(\d+)..(\d+)")
    # no need for disulfide length, effectively nonsensical
    #raw_data['disulfide_length'] = raw_data['Disulfide bond'].apply(regex_length, regex_form=r"(?:DISULFID )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Disulfide bond'])
    
    raw_data['cross_link_zones'] = raw_data['Cross-link'].apply(regex_zones, regex_form=r"(?:CROSSLNK )(\d+)")
    raw_data = raw_data.drop(columns=['Cross-link'])

    raw_data['percent_glycosylated'] = [percent_containing(raw_data.Glycosylation[i], raw_data.Length[i], 'CARBOHYD') for i in range(len(raw_data.Length))]

    raw_data['init_methionine'] = raw_data['Initiator methionine'].apply(regex_zones, regex_form=r"(?:INIT_MET )(\d+)")
    raw_data = raw_data.drop(columns=['Initiator methionine'])

    raw_data['signal_peptide_zones'] = raw_data['Signal peptide'].apply(regex_zones, regex_form=r"(?:SIGNAL )(\d+)..(\d+)")
    raw_data['signal_peptide_length'] = raw_data['Signal peptide'].apply(regex_length, regex_form=r"(?:SIGNAL )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Signal peptide'])

    raw_data['propeptide_zones'] = raw_data['Propeptide'].apply(regex_zones, regex_form=r"(?:PROPEP )(\d+)..(\d+)")
    raw_data['propeptide_length'] = raw_data['Propeptide'].apply(regex_length, regex_form=r"(?:PROPEP )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Propeptide'])

    raw_data['helix_zones'] = raw_data['Helix'].apply(regex_zones, regex_form=r"(?:HELIX )(\d+)..(\d+)")
    raw_data['helix_length'] = raw_data['Helix'].apply(regex_length, regex_form=r"(?:HELIX )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Helix'])

    raw_data['turn_zones'] = raw_data['Turn'].apply(regex_zones, regex_form=r"(?:TURN )(\d+)..(\d+)")
    raw_data['turn_length'] = raw_data['Turn'].apply(regex_length, regex_form=r"(?:TURN )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Turn'])

    raw_data['beta_strand_zones'] = raw_data['Beta strand'].apply(regex_zones, regex_form=r"(?:STRAND )(\d+)..(\d+)")
    raw_data['beta_strand_length'] = raw_data['Beta strand'].apply(regex_length, regex_form=r"(?:STRAND )(\d+)..(\d+)") / raw_data['Length']
    raw_data = raw_data.drop(columns=['Beta strand'])

    

    drop_list = ['Absorption', 'Activity regulation','Catalytic activity', 'ChEBI', 'ChEBI (Catalytic activity)', 'ChEBI (Cofactor)', 'ChEBI IDs', 'Chain', 'Compositional bias', 'Domain [CC]', 'Domain [FT]', 'EC number', 'Function [CC]', 'Gene ontology (GO)', 'Gene ontology (biological process)', 'Gene ontology (molecular function)', 'Gene ontology IDs', 'Induction', 'Interacts with', 'Kinetics', 'Motif', 'Pathway', 'Peptide', 'Post-translational modification', 'Protein families', 'Redox potential', 'Region', 'Repeat', 'Rhea Ids',  'Status', 'Site', 'Subcellular location [CC]', 'Subunit structure [CC]', 'Temperature dependence', 'Tissue specificity', 'Transit peptide', 'pH dependence']
    cleaned_data = raw_data.drop(columns=drop_list)
    # cleaned_data = raw_data[['Entry', 'Protein names', 'Sequence', 'Length', 'Mass', 'dna_binding_length',
    #                        'Metal binding', 'calcium_binding_zones', 'calcium_binding_length', 'Nucleotide binding', "nucleotide_binding_zones",
    #                        "nucleotide_binding_length", 'percent_glycosylated']]

    cleaned_data = cleaned_data.fillna(0)
    cleaned_data = metal_binding_analysis(cleaned_data)
    cleaned_data = nucleotide_types(cleaned_data)
    cleaned_data = cofactor_analysis(cleaned_data)
    cleaned_data = act_site_analysis(cleaned_data)

    sequences = raw_data['Sequence']
    seq_data = pd.DataFrame([])
    first_pass = True
    for seq in sequences:
        analyzed_seq = ProteinAnalysis(seq)
        aaCount = analyzed_seq.count_amino_acids()
        aa_percent = analyzed_seq.get_amino_acids_percent()
        reg_aa = regularize_aa(aa_percent)# dict
        aa_values = alphabetize_aa(reg_aa)
        mw = analyzed_seq.molecular_weight() #MW 
        aromat = analyzed_seq.aromaticity() #
        instab = analyzed_seq.instability_index() # float
        flex = analyzed_seq.flexibility() # returns a list
        iso = analyzed_seq.isoelectric_point() 
        secStruct = analyzed_seq.secondary_structure_fraction() # tuple of three floats (helix, turn, sheet)
        secStruct_disorder = 1 - sum(secStruct)
        flex_stat = (np.mean(flex), np.std(flex), np.var(flex), np.max(flex), np.min(flex), np.median(flex))
        gravy = analyzed_seq.gravy()
        temp_df = pd.DataFrame([[seq, *aa_values, mw, aromat, instab, *flex_stat, iso, *secStruct, secStruct_disorder, gravy]])
        
        

        if first_pass:
            seq_data = temp_df
            first_pass = False
        else:
            seq_data = seq_data.append(temp_df)
    
    aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    aa_list.sort()
    names_aa = ['frac_aa_' + i for i in aa_list]
    col_names = ['Sequence', *names_aa, 'molecular_weight', 'aromaticity', 'instability_index', 'flexibility_mean', 'flexibility_std', 'flexibility_var',
                 'flexibility_max', 'flexibility_min', 'flexibility_median', 'isoelectric_point', 'secondary_structure_fraction_helix',
                 'secondary_structure_fraction_turn', 'secondary_structure_fraction_sheet', 'secondary_structure_fraction_disordered', 'gravy']
        
    seq_data.columns = col_names
    seq_data=seq_data.reset_index(drop=True)

    cleaned_data = pd.merge(cleaned_data, seq_data, on='Sequence')
    cleaned_data = cleaned_data.drop(columns=['Glycosylation', 'DNA binding'])

    return cleaned_data

def accession_expansion(data):

    for i in data['Accession']:
        if ';' in i: # find lists 
            corona_value = int(data.loc[data['Accession'] == i]['Corona'].to_numpy()) # see what their value is
            to_add_list = i.split(';') # break list

            for j in to_add_list: # go through each list member
                data = data.append({'Accession':j, 'Corona':corona_value}, ignore_index=True) # make a new column for each

            data = data.loc[data['Accession'] != i] # delete the original index (may not be needed will be auto deleted later anyway)
    return data


def confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    conf_int = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return conf_int

if __name__ == "__main__":
    print(type(ProteinAnalysis))
    print('found biopython')