import pandas as pd
import json
import urllib.request
import numpy as np

def uniprot_data_scraping(uniprot_list, data_to_pull='standard'):
    """
        Searches the uniprot database based on the supplied list of uniprot numbers and returns a Pandas DataFrame
        with the found data. Designed to suppress errors but print that it found one.

        :param accession_number: Uniprot accession number as a string
        :type accession_number: str
        
        :param data_to_pull: Which columns should be pulled from the UniProt Database. 
        This does not guarantee the data will populate. Rather this attempts to find the above data. Choices: all, standard, or list of columns.

        :return: DataFrame 
    """
    
    # Pull in lookup codes
    query_lookup_table = pd.read_excel('helper_data/query_lookup_table.xlsx', header=0)
    
    # Identify the columns that we want to search for
    if data_to_pull == 'standard': # standard columns
        data_to_pull = ['Entry', 'Protein names', 'Status', 'Protein families', 'Length',
       'Mass', 'Sequence', 'Binding site', 'Calcium binding', 'DNA binding',
       'Metal binding', 'Nucleotide binding', 'Site', 'Function [CC]',
       'Absorption', 'Active site',
       'Catalytic activity', 'Cofactor', 'EC number', 'Kinetics', 'Pathway',
       'pH dependence', 'Redox potential', 'Rhea Ids',
       'Temperature dependence', 'Interacts with', 'Subunit structure [CC]',
       'Induction', 'Tissue specificity', 'Gene ontology (biological process)',
       'Gene ontology (GO)', 'Gene ontology (molecular function)',
       'Gene ontology IDs', 'ChEBI', 'ChEBI (Catalytic activity)',
       'ChEBI (Cofactor)', 'ChEBI IDs', 'Intramembrane',
       'Subcellular location [CC]', 'Transmembrane', 'Topological domain',
       'Chain', 'Cross-link', 'Disulfide bond', 'Glycosylation',
       'Initiator methionine', 'Lipidation', 'Modified residue', 'Peptide',
       'Propeptide', 'Post-translational modification', 'Signal peptide',
       'Transit peptide', 'Beta strand', 'Helix', 'Turn', 'Coiled coil',
       'Compositional bias', 'Domain [CC]', 'Domain [FT]', 'Motif', 'Region',
       'Repeat', 'Zinc finger']
    
    elif data_to_pull == 'all': # all columns that we have keys for
        data_to_pull=query_lookup_table['Column names as displayed on website'].values
        
    elif type(data_to_pull) == list: # specific names that we can lookup
        data_to_pull = data_to_pull # do nothing
        
    else: # doesnt fit the specified values
        raise ValueError('Change data_to_pull to "any", "standard" or a list of column names')
    
    
    col_string = '' #create the search column string
    first_pass=True
    for col in data_to_pull: # connect together the columns that need to be separated by commas 
        if first_pass:
            col_string = query_lookup_table.loc[query_lookup_table['Column names as displayed on website'] == col,
                                 'Column names as displayed in URL'].values[0]
            first_pass = False
            
        else:
            try: col_string += ','+query_lookup_table.loc[query_lookup_table['Column names as displayed on website'] == col,
                                 'Column names as displayed in URL'].values[0]
                
            except: print(f"Could not find URL Name for {col}")
    
    col_string = col_string.replace(' ', '%20') #get rid of spaces
    
    # create prefix and suffix to web address 
    query_string_suffix= "&format=tab&columns="+col_string
    query_string_prefix='https://www.uniprot.org/uniprot/?query='
    
    first_pass=True
    for uniprot in uniprot_list: #go through list of all proteins
        if first_pass: # initialize the dataframe 
            try:
                total_data = pd.read_csv(query_string_prefix+uniprot+query_string_suffix, sep='\t',  thousands=',') # get data
                
                if total_data.shape[0] > 1: # Only take correct row, sometimes UniProt sends back more than 1 response
                    total_data = total_data.loc[total_data.Entry==uniprot]
                first_pass = False # dont come back here  
                
            except: 
                print(f"Could not complete URL request for {uniprot}")
            
            
        else:
            try: 
                revolve_data = pd.read_csv(query_string_prefix+uniprot+query_string_suffix, sep='\t',  thousands=',') # get data
                
                if revolve_data.shape[0] > 1: # Only take correct row, sometimes UniProt sends back more than 1 response
                    revolve_data = revolve_data.loc[revolve_data.Entry==uniprot]
                    print(revolve_data.shape)
                total_data = total_data.append(revolve_data, ignore_index=True) 
                
            except: 
                print(f"Could not complete URL request for {uniprot}")
            
            
    
    total_data=total_data.fillna(0) # fill nans with zeros
    
    return total_data
    
def protein_data_scraping_fasta(accession_number):
    """
        Searches the interpro database based on the uniprot accession number and returns a JSON with the found data.
        Designed to suppress errors but print that it found one.

        :param accession_number: Uniprot accession number as a string
        :type accession_number: str

        :return: FASTA data as a Sting
    """

    # url for API
    base_url = 'https://www.uniprot.org/uniprot/'
    accession_number = accession_number.replace(' ', '')
    current_url = base_url + accession_number + '.fasta'
    print(current_url)
    
    try:
        # Fetch data
        with urllib.request.urlopen(current_url) as url:
            data = url.read().decode()
            return data

    except json.JSONDecodeError:  # if data cant be found just pass it
        print(f'Data not found for {accession_number} (Decode error)')
        return None

def netsurfp_1point1_data_processing(unique_id_list, complete_netsurfp_df):
    first_pass=True
    for i in unique_id_list:

        # get correct values
        boolarray = complete_netsurfp_df['sequence name'].str.contains(i)
        filtered = complete_netsurfp_df[boolarray]
        
        # total amino acids
        total_aa = filtered.shape[0]
        
        #count exposed residues
        exposed = filtered[filtered['class assignment'] == 'E']
        exposed_count = exposed.shape[0]
        if total_aa == 0:
            print(i)
            frac_exposed=0
        else:
            frac_exposed = exposed_count / total_aa
        
        aa_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 
                'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0,
                'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
        
        if frac_exposed == 0: # bypass any that have had no exposed amino acids
            aa_exposed_frac = 0
            
        else:
            for j in exposed['amino acid']: # go through and count what needs to be added
                aa_dict[j] += 1

            aa_exposed_frac_total = {'fraction_total_exposed_'+key: value / total_aa for key, value in aa_dict.items()}
            aa_exposed_frac_exposed = {'fraction_exposed_exposed_'+key: value / exposed_count for key, value in aa_dict.items()}
        
        nonpolar_aa = ['A', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W']
        nonpolar_counts = 0
        for k in nonpolar_aa:
            nonpolar_counts += aa_dict[k]
        
        if nonpolar_counts == 0: 
            nonpolar_exposed_frac = 0
        else:
            nonpolar_exposed_frac_exposed = nonpolar_counts / exposed_count
            nonpolar_exposed_frac_total = nonpolar_counts / total_aa
            
        data_to_update ={'entry': i, 'fraction_exposed':np.around(frac_exposed, 3), 'fraction_buried': np.around(1-frac_exposed, 3),
                        'fraction_exposed_nonpolar_total':nonpolar_exposed_frac_total,
                        'fraction_exposed_nonpolar_exposed':nonpolar_exposed_frac_exposed,
                        "rsa_mean": np.around(filtered['relative surface accessibility'].mean(), 3),
                        "rsa_median": np.around(filtered['relative surface accessibility'].median(), 3),
                        "rsa_std": np.around(filtered['relative surface accessibility'].std(), 3),
                        "asa_sum": np.around(filtered['absolute surface accessibility'].sum(), 3),
                        **aa_exposed_frac_total, **aa_exposed_frac_exposed}
        
        
        if first_pass:
            netsurfp_processed_data = pd.DataFrame()
            netsurfp_processed_data = pd.DataFrame.from_dict(data_to_update, orient='index').transpose()
            first_pass = False
        
        else:
            netsurfp_processed_data = netsurfp_processed_data.append(pd.DataFrame.from_dict(data_to_update, orient='index').transpose(), ignore_index=True)


    return netsurfp_processed_data

def netsurfp_2_data_processing(unique_id_list, complete_netsurfp_df):
    first_pass=True
    for i in unique_id_list:
        # print('hi')
        # get correct values
        boolarray = complete_netsurfp_df['id'].str.contains(i)
        filtered = complete_netsurfp_df[boolarray]
        filtered['class assignment'] = np.where(filtered.rsa > 0.25, 'E', 'B')
        # total amino acids
        total_aa = filtered.shape[0]
        
        #count exposed residues
        exposed = filtered[filtered['class assignment'] == 'E']
        exposed_count = exposed.shape[0]
        if total_aa == 0:
            print('Error at ',i)

            frac_exposed=0
            frac_buried = 0

        else:
            frac_exposed = exposed_count / total_aa
            frac_buried = 1-frac_exposed
        
        aa_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 
                'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0,
                'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}

        if frac_exposed == 0: # bypass any that have had no exposed amino acids
            aa_exposed_frac = 0
            aa_exposed_frac_total = {'fraction_total_exposed_'+key: 0 for key, value in aa_dict.items()}
            aa_exposed_frac_exposed = {'fraction_exposed_exposed_'+key: 0 for key, value in aa_dict.items()}
            
        else:
            for j in exposed['seq']: # go through and count what needs to be added
                aa_dict[j] += 1

            aa_exposed_frac_total = {'fraction_total_exposed_'+key: value / total_aa for key, value in aa_dict.items()}
            aa_exposed_frac_exposed = {'fraction_exposed_exposed_'+key: value / exposed_count for key, value in aa_dict.items()}
        
        nonpolar_aa = ['A', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W']
        nonpolar_counts = 0
        for k in nonpolar_aa:
            nonpolar_counts += aa_dict[k]
        
        if nonpolar_counts == 0: 
            nonpolar_exposed_frac_exposed = 0
            nonpolar_exposed_frac_total = 0

        else:
            nonpolar_exposed_frac_exposed = nonpolar_counts / exposed_count
            nonpolar_exposed_frac_total = nonpolar_counts / total_aa
            
        data_to_update ={'entry': i, 'fraction_exposed':np.around(frac_exposed, 3), 'fraction_buried': np.around(frac_buried, 3),
                        'fraction_exposed_nonpolar_total':nonpolar_exposed_frac_total,
                        'fraction_exposed_nonpolar_exposed':nonpolar_exposed_frac_exposed,
                        "rsa_mean": np.around(filtered['rsa'].mean(), 3),
                        "rsa_median": np.around(filtered['rsa'].median(), 3),
                        "rsa_std": np.around(filtered['rsa'].std(), 3),
                        "asa_sum": np.around(filtered['asa'].sum(), 3),
                        **aa_exposed_frac_total, **aa_exposed_frac_exposed,
                        'nsp_secondary_structure_coil':np.around(np.sum(np.where(filtered.q3.str.contains('C'), True, False))/filtered.shape[0], 3),
                        'nsp_secondary_structure_sheet':np.around(np.sum(np.where(filtered.q3.str.contains('E'), True, False))/filtered.shape[0], 3),
                        'nsp_secondary_structure_helix':np.around(np.sum(np.where(filtered.q3.str.contains('H'), True, False))/filtered.shape[0], 3),
                        'nsp_disordered':np.around(np.sum(filtered.disorder.to_numpy() >= 0.5)/filtered.shape[0],3)}
        
        
        if first_pass:
            netsurfp_processed_data = pd.DataFrame()
            netsurfp_processed_data = pd.DataFrame.from_dict(data_to_update, orient='index').transpose()
            first_pass = False
        
        else:
            netsurfp_processed_data = netsurfp_processed_data.append(pd.DataFrame.from_dict(data_to_update, orient='index').transpose(), ignore_index=True)


    return netsurfp_processed_data

if __name__ == "__main__": 
    test_df = uniprot_data_scraping(['P61823'],data_to_pull='standard')
    print(test_df)