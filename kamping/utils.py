import os

import pandas as pd


def read_all_tsv_files(directory):
    '''
    Read all files in the directory as one dataframe.
    '''
    # list a files except those start with "."

    all_files = os.listdir(directory)
    df = pd.DataFrame()
    for file in all_files:
        if file.startswith('.'):
            continue
        temp_df = pd.read_csv(os.path.join(directory, file), sep='\t')
        # add a new column "source" to the DataFrame and set it to the file name
        temp_df['source'] = file.removesuffix('.tsv')
        df = pd.concat([df, temp_df])
    # reset the index of the DataFrame
    df.reset_index(drop=True, inplace=True)
    return df
