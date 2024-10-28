import logging

import pandas as pd
from pathlib import Path

class ProteinMetabliteParser:
    '''
    Create a parser to convert the output of the mixed command into a
    metabolite-protein interaction file.

    usage:
    parser = proteinMetabliteParser()
    parser.parse(file='mixed.tsv', wd='.')

    return: None
    '''
    def __init__(self, verbose: bool = False, keep_glycan=False, keep_PPI: bool = False):
        self.verbose = verbose
        self.keep_PPI = keep_PPI
        self.keep_glycan = keep_glycan


    def parse_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        # expand the relation in the DataFrame
        # create an empty dataframe
        new_df = pd.DataFrame()

        # if row number greater than 0:
        if (df.shape[0]) == 0:
            # throw an error
            return new_df



        for _, row in df.iterrows():
            # for PPrel type, if the name is compound (intermediated by metabolite), expand the relation
            if row['type'] == 'PPrel':
                if row['name'] == 'compound':
                    new_row_df = expand_relation_PPrel(row)
                    new_df = pd.concat([new_df, new_row_df], ignore_index=True)
                # if the name is not compound, keep the row when the keep_PPI is True
                elif self.keep_PPI:
                    new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)
            # for ECrel type, expand the relation
            elif row['type'] == 'ECrel' and row['name'] == 'compound':
                   # hsa05140.xml there is ECrel with activation... error?
                    #         <relation entry1="23" entry2="27" type="ECrel">
                    # <subtype name="compound" value="28"/>
                    # <subtype name="activation" value="--&gt;"/>
                    # <subtype name="indirect effect" value="..&gt;"/>
                new_row_df = expand_relation_ECrel(row)
                new_df = pd.concat([new_df, new_row_df])
            # for PCrel type, keep the row
            elif row['type'] == 'PCrel':
                new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)
            elif row['type'] == 'GErel':
                if self.keep_PPI:
                    new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

        if new_df.shape[0] == 0:
            raise ValueError(f'input dataframe does not contain data to parse!')

        # Create a new DataFrame from the list of new rows
        new_df.reset_index(inplace=True, drop=True)

        # remove rows with entry1 or entry2 start with gl
        if not self.keep_glycan:
            keep = ~ (new_df['entry1'].str.startswith('gl') | new_df['entry2'].str.startswith('gl'))
            new_df = new_df[keep]

        return new_df

    def parse_file(self, file:str, wd:str) -> None:
        # the input file should be the output of mixed command
        if not Path(file).exists():
            raise FileNotFoundError(f'File {file} not found!')

        # load the file into a DataFrame
        df = pd.read_csv(file, sep='\t')
        try:
            df = self.parse_dataframe(df)
        except ValueError as e:
            logging.warning(f'check {file}, {e.args[0]}')
            return
        # export the DataFrame to a TSV file
        df.to_csv(Path(wd) / f'{Path(file).stem}_mpi.tsv', sep='\t', index=False)


def expand_relation_ECrel(row: pd.Series):
    '''
    helper function to expand the relation for ECrel type
    '''
    new_row1 = row.copy()
    new_row1['entry2'] = row['value']
    new_row1['type'] = 'PCrel'
    new_row1['value'] = "custom"
    new_row1['name'] = "enzyme-enzyme expansion"
    new_row1['entry2_type'] = 'compound'

    new_row2 = row.copy()
    new_row2['entry2'] = row['entry2']
    new_row2['entry1'] = row['value']
    new_row2['type'] = 'PCrel'
    new_row2['value'] = 'custom'
    new_row2['name'] = 'enzyme-enzyme expansion'
    new_row2['entry1_type'] = 'compound'

    # combined two series into a DataFrame
    df = pd.concat([new_row1, new_row2], axis=1).transpose()
    return df

def expand_relation_PPrel(row: pd.Series):
    '''
    helper function to expand the relation for PCrel type
    '''
    new_row1 = row.copy()
    new_row1['entry2'] = row['value']
    new_row1['type'] = 'PCrel'
    new_row1['value'] = 'custom'
    new_row1['name'] = 'protein-protein expansion'
    new_row1['entry2_type'] = 'compound'

    new_row2 = row.copy()
    new_row2['entry2'] = row['entry2']
    new_row2['entry1'] = row['value']
    new_row2['type'] = 'PCrel'
    new_row2['value'] = 'custom'
    new_row2['name'] = 'protein-protein expansion'
    new_row2['entry1_type'] = 'compound'

    # combined two series into a DataFrame
    df = pd.concat([new_row1, new_row2], axis=1).transpose()

    return df

def remove_suffix(entry: str):
    '''
    helper function to remove the suffix from the entry
    '''
    return entry.split('-')[0]




def parse_to_mpi(input_data: str, wd: Path, keep_glycan=False, keep_PPI=False, verbose: bool = False):
    '''
    Converts a file or folder of mixed pathways tsv files to metabolite-protein interactions tsv files.
    The generalized metabolite-protein interactions includes the following types:
    - ECrel: Enzyme-Enzyme relations which is successive enzyme reaction connected by a compound.
    - PCrel: Protein-Compound relations which is a protein connected to a compound.
    - PPrel: Protein-Protein relations with a compound as intermediate.

    Parameters:
    input_data: str: The file or folder of mixed pathways tsv files.
    wd: Path: The directory to save the resulting files.
    keep_PPI: bool: Keep the protein-protein interactions (include PPrel without compound as intermidate and
    GErel in the resulting file.
    keep_glycan: bool: Keep relations with one of the entries starting with 'gl' in the resulting file. When used in successive
    metabolite embedding analysis, glycan relations should be removed since no MDL molfile is available .
    verbose: bool: Print verbose messages
    '''
    protein_metabolite_parser = ProteinMetabliteParser(keep_PPI=keep_PPI,
                                                       keep_glycan=keep_glycan, verbose=verbose)
    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.tsv'):
            try:
                protein_metabolite_parser.parse_file(file=file, wd=wd)
            except ValueError as e:
                logging.warning(e)
                continue
    else:
        protein_metabolite_parser.parse_file(file=input_data, wd=wd)