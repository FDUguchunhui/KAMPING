import pandas as pd
from knext.utils import FileNotFound
from pathlib import Path

class MPIParser:

    def parse(self, file: str):
        # the input file should be the output of mixed command
        if not Path(file).exists():
            raise FileNotFound(f'File {file} not found!')

        df = pd.read_csv(file, sep='\t')

        # remove rows with type "maplink", "GErel", and "PPrel"
        df = df[~df['type'].isin(['maplink', 'GErel', 'PPrel'])]

        #


def expand_EC_relation(df: pd.DataFrame):

    # Create a list to store the new rows
    new_rows = []

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        # If the "name" field is "compound", create two new rows
        if row['name'] == 'compound':
            new_row1 = row.copy()
            new_row1['entry2'] = row['value']
            new_rows.append(new_row1)

            new_row2 = row.copy()
            new_row2['entry1'] = row['entry2']
            new_row2['entry2'] = row['value']
            new_rows.append(new_row2)
        else:
            new_rows.append(row)

    # Create a new DataFrame from the list of new rows
    new_df = pd.concat(new_rows, axis=1).transpose()
    new_df.reset_index(inplace=True, drop=True)
    return new_df