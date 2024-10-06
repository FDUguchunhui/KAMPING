import re
from collections import defaultdict
import urllib.request as request
import numpy as np
import pandas as pd

def entry_id_conv_dict(root, unique=False):
    # This dictionary is the unique version
    # Every item is unique to reveal subgraphs
    entries = [(entry.get('id'), entry.get('name'), entry.get('type')) for entry in root.findall('entry')]
    entries = pd.DataFrame(entries, columns=['id', 'name', 'type'])
    # names are separated by space, make it as a list
    entries['name'] = entries['name'].str.split()

    # get unique names for each entry
    if unique:
        name = [[name + '-' + id for name in names ] for id, names in zip(entries['id'], entries['name'])]
    else:
        name = entries['name']

    conversion_dictionary = dict(zip(entries['id'], name))

    return conversion_dictionary



def names_dict(root, organism, conversion_dictionary):
    '''
    The `names_dict` function parses the names in the KEGG API XML file for the given root. It does the following:
        1. Iterates through the `relation` elements in the XML to collect `entry1` and `entry2` attributes.
        2. Combines these entries and maps them using the provided `conversion_dictionary`.
        3. For each unique entry, it fetches additional information from the KEGG API based on the type of entry (gene, compound, or pathway).
        4. Constructs a dictionary (`dd`) where keys are the entries and values are their corresponding names or descriptions fetched from the KEGG API.

This function helps in mapping KEGG entries to their human-readable names or descriptions.
    '''
    # d = conv_dict_unique(root)
    # d = conv_dict(root)
    e1 = []
    e2 = []
    for entry in root.findall('relation'):
        e1.append(entry.get('entry1'))
        e2.append(entry.get('entry2'))
    e = e1 + e2
    e_list = [conversion_dictionary[entry].split(' ') for entry in e]
    e_list1 = [l for sublist in e_list for l in sublist]
    e_conv = set(e_list1)
    dd = {}
    for n in e_conv:
        # Uses organism code since there are pathways, undefined, and others that will cause
        # an error if used here
        if n.startswith(organism):
            # Remove, if necessary, any terminal modifiers to avoid an error in api call
            n4url = re.sub(r'-[0-9]+', '', n)
            # Uses find to get gene info since other api tools give error
            url = 'https://rest.kegg.jp/find/genes/%s/'
            response = request.urlopen(url % n4url).read().decode('utf-8')
            split_response = response.split('\n')
            # Find only the query gene if given back several accessions that are similar
            s = filter(lambda x: x.startswith(n4url + '\t'), split_response)
            try:
                # Adds to dictionary the end entry, which is the written out name
                dd[n] = re.sub('^ ', '', list(s)[0].split(';')[1])
            except IndexError:
                # Some genes only have a name and no description
                dd[n] = split_response[0].split('\t')[1]
        # Only obtains compounds
        elif n.startswith('cpd:'):
            # Remove terminal modifiers, which are always added to compounds
            # unless a non-unique mixed pathway is chosen
            n4url = re.sub(r'-[0-9]+', '', n)
            # Uses find to get gene info since other api tools give error
            url = 'https://rest.kegg.jp/find/compound/%s'
            response = request.urlopen(url % n4url).read().decode('utf-8')
            subbed_response = re.sub(r'%s\t' % n4url, '', response)
            try:
                # Find only the query compound if given back several accessions that are similar
                split_response = re.sub('^ ', '', subbed_response.strip('\n').split(';')[1])
            except IndexError:
                # Some compounds only have one name
                split_response = subbed_response.strip('\n')
            # Adds to dictionary the end entry, which is the written out name
            dd[n] = split_response
        elif n.startswith('path:'):
            n4url1 = re.sub(r'-[0-9]+', '', n)
            n4url2 = re.sub(r'path:{}'.format(organism), '', n4url1)
            url = 'https://rest.kegg.jp/find/pathway/%s'
            response = request.urlopen(url % n4url2).read().decode('utf-8').strip('\n').split('\t')
            try:
                dd[n] = response[1]
            except IndexError:
                # One pathway has no metadata and gives an error if line not included
                dd[n] = np.nan
        else:
            dd[n] = np.nan
    return dd

def _parse_entries(root):
    '''
    Parses the entries in the KEGG API XML file for the given root.
    Returns the entry id, name, and type.
    '''
    entry_dict = defaultdict(list)
    for entries in root.findall('entry'):
        for key, items in entries.attrib.items():
            entry_dict[key].append(items)

    entry_id=[]
    entry_name=[]
    entry_type=[]
    for key, items in entry_dict.items():
        if key == 'id':
            for i in items:
                entry_id.append(i)
        if key == 'name':
            for i in items:
                entry_name.append(i)
        if key == 'type':
            for i in items:
                entry_type.append(i)

    return entry_id, entry_name, entry_type


def get_conversion_dictionary(species, target):
    '''
    Convert KEGG gene IDs to either NCBI gene IDs or UniProt IDs.
    '''
    if target == 'uniprot':
        url = 'http://rest.kegg.jp/conv/%s/uniprot'
    elif target == 'ncbi':
        url = 'http://rest.kegg.jp/conv/%s/ncbi-geneid'
    response = request.urlopen(url % species).read().decode('utf-8')
    response = response.rstrip().rsplit('\n')
    kegg = []
    uniprot = []
    for resp in response:
        uniprot.append(resp.rsplit()[0])
        kegg.append(resp.rsplit()[1])
    d = {}
    for key, value in zip(kegg, uniprot):
        if key not in d:
            d[key] = [value]
        else:
            d[key].append(value)
    return d
