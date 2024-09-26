import pandas as pd
from pathlib import Path
import typer

class proteinMetabliteParser:
    '''
    Create a parser to convert the output of the mixed command into a
    metabolite-protein interaction file.

    usage:
    parser = proteinMetabliteParser()
    parser.parse(file='mixed.tsv', wd='.')

    return: None
    '''
    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def parse(self, file: str, wd: Path):
        # the input file should be the output of mixed command
        if not Path(file).exists():
            raise FileNotFoundError(f'File {file} not found!')

        # load the file into a DataFrame
        df = pd.read_csv(file, sep='\t')
        # remove rows with type "maplink"
        # which will leave "ECrel", "GErel", and "PPrel", and "PCrel"

        # if row number greater than 0:
        if (df.shape[0]) == 0:
            # throw an error
            raise ValueError(f'{file} does not contain data to parse!')

        # expand the relation in the DataFrame
        # create an empty dataframe
        new_df = pd.DataFrame()

        for _, row in df.iterrows():
            if row['name'] == 'compound':
                if row['type'] == 'PPrel':
                    # convert the row to a DataFrame
                    new_row_df = expand_relation_PPrel(row)
                    new_df = pd.concat([new_df, new_row_df])
                if row['type'] == 'ECrel':
                    new_row_df = expand_relation_ECrel(row)
                    new_df = pd.concat([new_df, new_row_df])
            elif row['type'] == 'PCrel':
                # remove when compound is glycan
                if row['entry1'].startswith('cpd:') or row['entry2'].startswith('cpd:'):
                    new_df = pd.concat([new_df, row.to_frame().T],  ignore_index=True)
                else:
                    continue
            else:
                new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

        if new_df.shape[0] == 0:
            raise ValueError(f'No data to parse in {file}!')

        # Create a new DataFrame from the list of new rows
        new_df.reset_index(inplace=True, drop=True)

        # remove the suffix from the entry
        new_df['entry1'] = new_df['entry1'].apply(lambda x: remove_suffix(x))
        new_df['entry2'] = new_df['entry2'].apply(lambda x: remove_suffix(x))


# export the DataFrame to a TSV file
        # get file name without extension
        file_name = Path(file).stem
        wd = Path(wd)
        new_df.to_csv(wd / f'{file_name}_mpi.tsv', sep='\t', index=False)

        return None

def expand_relation_ECrel(row: pd.Series):
    '''
    helper function to expand the relation for ECrel type
    '''

    new_row1 = row.copy()
    new_row1['entry2'] = row['value']
    new_row1['type'] = 'PCrel'
    new_row1['value'] = "custom"
    new_row1['name'] = "enzyme-enzyme expansion"

    new_row2 = row.copy()
    new_row2['entry1'] = row['entry2']
    new_row2['entry2'] = row['value']
    new_row2['type'] = 'PCrel'
    new_row2['value'] = 'custom'
    new_row2['name'] = 'enzyme-enzyme expansion'

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

    new_row2 = row.copy()
    new_row2['entry1'] = row['entry2']
    new_row2['entry2'] = row['value']
    new_row2['type'] = 'PCrel'
    new_row2['value'] = 'custom'
    new_row2['name'] = 'protein-protein expansion'

    # combined two series into a DataFrame
    df = pd.concat([new_row1, new_row2], axis=1).transpose()

    return df

def remove_suffix(entry: str):
    '''
    helper function to remove the suffix from the entry
    '''
    return entry.split('-')[0]


def parse_to_mpi(input_data: str, wd: Path,  verbose: bool = False):
    '''
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis.
    '''
    protein_metabolite_parser = proteinMetabliteParser()
    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.tsv'):
            try:
                protein_metabolite_parser.parse(file=file, wd=wd)
            except ValueError as e:
                typer.echo(typer.style(e, fg=typer.colors.RED, bold=True))
                continue
    else:
        protein_metabolite_parser.parse(file=input_data, wd=wd)

