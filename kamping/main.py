#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 20:49:11 2023

@author: everest_castaneda1
"""
import logging
from enum import Enum

import click

from kamping.parser.convert import Converter

logging.basicConfig(level=logging.INFO)
import sys
from typing import Union
from typing_extensions import  Annotated
import typer
from pathlib import Path

from kamping.parser.network import KeggGraph
from kamping.parser.utils import kgml


@click.group()
def cli():
    pass


@cli.command()
@click.argument('species', type=str)
@click.option('--out_dir', type=str, default=None, help='Directory to save results. If not provided, results will be saved in the current working directory.')
def get_kgml(species: str, out_dir: Union[str, None] = None):
    """
    Acquires all KGML files for a given species. Use a KEGG species
    identification, usually a 3 to 4 letter organism code, as input. Handles
    directories first, hence make sure you have a proper species identifier or
    it may return empty, junk folders. The results flag is to save all results 
    to a directory. For more information about KEGG organism codes, visit: 
    https://www.genome.jp/kegg/catalog/org_list.html
    """

    kgml(species, out_dir)

@cli.command()
@click.option('--type', type=click.Choice(['gene-only', 'mpi', 'original']), required=True, help='The type of network')
@click.argument('species', type=str)
@click.argument('input_data', type=str)
@click.option('--id_conversion', type=click.Choice(['ncbi', 'uniprot']), default=None, help='Convert KEGG gene id to which identifier')
@click.option('--unique', is_flag=True, help='Flag to return unique genes with terminal modifiers')
@click.option('--out_dir', type=str, default=None, help='Directory to save results. If not provided, results will be saved in the current working directory.')
@click.option('--verbose', is_flag=True, help='Flag to print progress')
def network(type: str, species: str, input_data: str, out_dir:str,
            id_conversion: Union[str, None] = None, unique: bool = False, verbose:bool = False):
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

    id_converter = None
    if id_conversion is None:
        id_converter = Converter(species, gene_target=id_conversion, unique=unique, verbose=verbose)

    if Path(input_data).is_dir():
        files = sorted(Path(input_data).glob('*.xml'))
        for file in files:
            try:
                logging.info(f'Parsing {file}...')
                interaction = KeggGraph(type=type, input_data=file,
                                        id_converter=id_converter,
                                        unique=unique,
                                        verbose=verbose)
                df_out = interaction.interaction
                df_out.to_csv(out_dir / f'{file.stem}.tsv', sep='\t', index=False)
            except Exception as e:
                typer.echo(typer.style(f'Error when parsing {file}: {e}', fg=typer.colors.RED, bold=True))
                continue
    else:
        logging.info(f'Parsing {input_data}...')
        interaction = KeggGraph(type=type,
                                input_data=input_data,
                                id_converter=id_converter,
                                unique=unique,
                                verbose=verbose)
        df_out = interaction.interaction
        df_out.to_csv(out_dir / f'{Path(input_data).stem}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    cli()