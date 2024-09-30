#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 20:49:11 2023

@author: everest_castaneda1
"""

import sys
from typing import Union

import typer
from pathlib import Path

from .kegg_parser.convert import genes_convert
from .kegg_parser.call import kgml
from .kegg_parser.pathway import parse_pathway

app = typer.Typer()

@app.command()
def get_kgml(species: str,
             results: Union[str, None] = typer.Option(None, help='Directory to save results. '
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

    kgml(species, results)

@app.command()
def pathway(input_data: str = typer.Argument(..., help='Path to KGML file or folder of KGML files'),
            mixed: bool = typer.Option(False, help='Flag to return mixed interactions.'),
            results: Union[str, None] = typer.Option(None, help='Directory to save results. '
                                                              'If not provided, results will be saved in the current working directory.'),
            unique: bool = typer.Option(False, help='Flag to return unique genes with terminal modifiers.'),
            verbose: bool = typer.Option(False, help='Flag to print progress.')
            ):
    """
    Converts a folder of KGML files or a single KGML file into a
    edgelist of genes that can be used in graph analysis. If -u/--unique flag
    is used genes are returned with terminal modifiers to enhance network
    visualization or analysis.
    """
    # work as a wrapper function with mixed=False call parse function parse the file(s)
    parse_pathway(input_data, results=results, mixed=mixed, unique=unique, verbose=False)


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
    app()
