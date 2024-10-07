import os
from pathlib import Path

import pandas as pd
import typer

from kamping.mol.utils import get_unique_compound_values, fetch_mol_file_string, get_smiles
from kamping.utils import read_all_tsv_files


def mol_scrapper(input_data:str, output_file: str):
    '''
    Function to scrap mol files from KEGG database
    '''

    # if wd is None, set it to current working directory
    if output_file is None:
        output_file = os.path.join(os.getcwd() + '/smile_files.tsv')
    else:
        # create the output directory if it does not exist
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # check if input_data is a folder
    if os.path.isdir(input_data):
        # get all files in the folder
        df = read_all_tsv_files(input_data)
        compounds = get_unique_compound_values(df)
    else:
        # get file name without extension
        df = pd.read_csv(input_data, sep='\t')
        compounds = get_unique_compound_values(df)

    # iterate over all compounds and fetch the mol file string
    if compounds is not None:
        # use typer to print the progress
        # open a file to save the smiles string
        with open(output_file, 'w') as f:
            # add header to the file
            f.write('id\tsmiles\n')
            for compound in compounds:
                # remove "cpd:" prefix
                compound = compound.replace('cpd:', '')
                typer.echo(typer.style(f'Now parsing: {compound}...', fg=typer.colors.GREEN, bold=False))
                mol_file_string = fetch_mol_file_string(compound)
                smiles_string = get_smiles(mol_file_string)
                # save the smiles string to a file as line
                f.write(f'{compound}\t{smiles_string}\n')

            # after writing all lines, close the file
            f.close()

        # print file process completion
        typer.echo(typer.style(f'Processing completed and output saved to {output_file}', fg=typer.colors.GREEN, bold=False))

if __name__ == '__main__':
    # for test
    mol_scrapper('data/mixed_output_PC', 'data/mol_files/smiles_file.tsv')