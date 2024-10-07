#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 20:49:11 2023

@author: everest_castaneda1
"""
import logging
from enum import Enum

logging.basicConfig(level=logging.INFO)
import sys
from typing import Union
from typing_extensions import  Annotated
import typer
from pathlib import Path

from kamping.parser.network import InteractionParser
from kamping.parser.call import kgml

app = typer.Typer()

class Type(str, Enum):
    gene_only = 'gene-only'
    MPI = 'MPI'
    original = 'original'

class Identifier_Conversion(str, Enum):
    uniprot = 'uniprot'
    ncbi = 'ncbi'
    none = None

@app.command()
def get_kgml(species: str = typer.Argument(..., help='the species to get kgml files'),
             out_dir: Union[str, None] = typer.Option(None, help='Directory to save results. '
                                                                 'If not provided, results will be saved to the current working directory.')
             ):
    """
    Acquires all KGML files for a given species. Use a KEGG species
    identification, usually a 3 to 4 letter organism code, as input. Handles
    directories first, hence make sure you have a proper species identifier or
    it may return empty, junk folders. The results flag is to save all results 
    to a directory. For more information about KEGG organism codes, visit: 
    https://www.genome.jp/kegg/catalog/org_list.html
    """

    kgml(species, out_dir)

#todo: check the default of boolean
@app.command()
def network(
            type: Annotated[Type, typer.Option(help='the type of network')],
            input_data: str = typer.Argument(..., help='Path to KGML file or folder of KGML files'),
            id_conversion: Annotated[Identifier_Conversion, typer.Option(help=' convert KEGG gene id to which identifier ')] = None,
            unique: Annotated[bool, typer.Option(help='Flag to return unique genes with terminal modifiers.')] = False,
            out_dir: Union[str, None] = typer.Option(None, help='Directory to save results. '
                                                              'If not provided, results will be saved in the current working directory.'),
            verbose: Annotated[bool, typer.Option(help='Flag to print progress.')] = False):
    """
    Converts a folder of KGML files or a single KGML file into a
    edgelist of genes that can be used in graph analysis. If -u/--unique flag
    is used genes are returned with terminal modifiers to enhance network
    visualization or analysis.


    """
    #todo: id_conversion dictionary should be created only once to improve performance
    if out_dir is None:
        out_dir = Path.cwd()
    else:
        out_dir = Path(out_dir)
        # create the output directory if it does not exist
        out_dir.mkdir(parents=True, exist_ok=True)

    if Path(input_data).is_dir():
        files = sorted(Path(input_data).glob('*.xml'))
        for file in files:
            try:
                logging.info(f'Parsing {file}...')
                gip = InteractionParser(type=type.value, input_data=file,
                                        id_conversion=id_conversion.value,
                                        unique=unique,
                                        verbose=verbose)
                df_out = gip.parse_file()
                df_out.to_csv(out_dir / f'{file.stem}.tsv', sep='\t', index=False)
            except Exception as e:
                typer.echo(typer.style(f'Error when parsing {file}: {e}', fg=typer.colors.RED, bold=True))
                continue
    else:
        logging.info(f'Parsing {input_data}...')
        gip = InteractionParser(type=type,
                                input_data=input_data,
                                id_conversion=id_conversion,
                                unique=unique,
                                verbose=verbose)
        df_out = gip.parse_file()
        df_out.to_csv(out_dir / f'{Path(input_data).stem}.tsv', sep='\t', index=False)


    # @cli.command()
    # @click.argument('file')
    # @click.option('-r', '--results', required = False)
    # def from_mixed_to_mpi(file: str, results: str = None):
    #     """
    #     Converts the mixed file to metabolite-protein interactions.
    #     """
    #     # remove the maplink, GErel, and PPrel
    #     df = MPIParser().parse(file)
    #     # export the DataFrame to a TSV file
    #     if results is not None:
    #         results = Path(results)
    #         if results.exists() == False:
    #             typer.echo(f'Directory {results} does not exist or is invalid. Please input a valid directory...')
    #             sys.exit()
    #         else:
    #             df.to_csv(results / 'mpi.tsv', sep='\t', index=False)
    #     else:
    #         wd = Path.cwd()
    #         results = wd / 'kgml_{}'.format(species)
    #         typer.echo(f'No output directory provided. All files will be saved to:\n{results}')
    #         results.mkdir(exist_ok = True)
    #         df.to_csv(results / 'mpi.tsv', sep='\t', index=False)



if __name__ == '__main__':
    network(input_data='data/kgml_hsa/hsa051d40.xml', type='MPI', id_conversion='uniprot', out_dir='data/converted')

