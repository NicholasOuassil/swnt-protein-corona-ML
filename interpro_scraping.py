import urllib.request
import json


def protein_data_scraping(accession_number):
    """
        Searches the interpro database based on the uniprot accession number and returns a JSON with the found data.
        Designed to suppress errors but print that it found one.

        :param accession_number: Uniprot accession number as a string
        :type accession_number: str

        :return: JSON / Dictionary
    """

    # url for API
    base_url = 'https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/'

    current_url = base_url + accession_number

    try:
        # Fetch data
        with urllib.request.urlopen(current_url) as url:
            data = json.loads(url.read().decode())  # unpack
            return data

    except json.JSONDecodeError:  # if data cant be found just pass it
        print(f'Data not found for {accession_number} (Decode error)')
        return None


def interpro_scraping(accession_number):
    """
        This function collects the InterPro numbers from the desired uniprot accession number. Input number as a string.

        :param accession_number: Uniprot accession number as a string
        :type accession_number: str

        :return: list of InterPro numbers
    """

    data = protein_data_scraping(accession_number)

    interpro_numbers = []  # data storage

    if data is not None:  # if there wasn't an error
        results = data['results']  # list of results
        count = data['count']  # find out how many intepro numbers there are

        for j in range(count):  # cycle through all found interpro numbers
            revolv_dict = results[j]  # unpacked dict (one per InterPro Number)
            interpro_numbers.append(revolv_dict['metadata']['accession'])  # add the found InterPro number

        return interpro_numbers

    else:  # if we couldn't find anything
        return interpro_numbers  # return a blank list


def interpro_scraping_pandas(input_data, column_of_interest, new_col_name):
    """
    Web scraping tool for translating Uniprot Accession Numbers to InterPro Numbers
    Expecting input data as a pandas dataframe and column of interest as the location string
    Returns the dataframe with a new column name
    """

    data_storage = []  # init list

    # Cycle through all input data
    for i in input_data[column_of_interest]:
        data_storage.append(interpro_scraping(i))  # add this urls data to the final list

    input_data[new_col_name] = data_storage  # add the new column
    return input_data




if __name__ == '__main__':
    print(interpro_scraping('Q9ZUA2'))  # test to see if it works
