#!/usr/bin/env python3

import csv
import operator

import chopin_kmd_hmdb_api_client.client
from chopin_kmd_hmdb_api_client.api.default import (
    api_annotation_get,
    api_compound_find,
    api_taxonomy_get,
)

import click

__version__ = "1.0.0"


kmd_hmdb_client = chopin_kmd_hmdb_api_client.client.Client(
    "https://kmd-hmdb-rest-api.metabolomics-chopin.e-metabohub.fr",
    verify_ssl=False,
    timeout=500,
)

find_compound = (
    lambda *args, **kwargs:
        api_compound_find.sync(*args, **kwargs, client=kmd_hmdb_client)
)
get_taxonomy = (
    lambda *args, **kwargs:
        api_taxonomy_get.sync(*args, **kwargs, client=kmd_hmdb_client)
)
get_annotation = (
    lambda *args, **kwargs:
        api_annotation_get.sync(*args, **kwargs, client=kmd_hmdb_client)
)

positive_adducts = [
    "M+H",
    "M+2H",
    "M+H+NH4",
    "M+H+Na",
    "M+H+K",
    "M+ACN+2H",
    "M+2Na",
    "M+H-2H2O",
    "M+H-H2O",
    "M+NH4",
    "M+Na",
    "M+CH3OH+H",
    "M+K",
    "M+ACN+H",
    "M+2Na-H",
    "M+IsoProp+H",
    "M+ACN+Na",
    "M+2K+H",
    "M+DMSO+H",
    "M+2ACN+H",
    "2M+H",
    "2M+NH4",
    "2M+Na",
    "2M+K",
]

negative_adducts = [
    "M-H",
    "M-2H",
    "M-H2O-H",
    "M+Cl",
    "M+FA-H",
    "M+Hac-H",
    "M-H+HCOONa",
    "M+Br",
    "M+TFA-H",
    "2M-H",
    "2M+FA-H",
    "2M+Hac-H",
]

adduct_choices = positive_adducts + negative_adducts

taxonomy_column_choices = [
    "class",
    "kingdom",
    "molecular_framework",
    "sub_class",
    "super_class",
    "id",
]

annotation_column_choices = [
    "adduct",
    "kendricks_mass",
    "kendricks_mass_defect",
    "monisotopic_molecular_weight",
    "nominal_mass",
    "polarity",
    "annotation_id",
]

compound_column_choices = [

    "database",
    "metabolite_name",
    "chemical_formula",
    "hmdb_id",
    "inchikey",
    "compound_id",
] + annotation_column_choices


@click.group()
def cli():
    pass


@cli.command(help="")
@click.option(
    "--version",
    is_flag=True,
)
@click.option(
    "--mz-ratio",
    default=[303.05],
    show_default=True,
    multiple=True,
    help="Provide the mz-ratio."
)
@click.option(
    "--database",
    default=["farid"],
    show_default=True,
    multiple=True,
    help="Provide the database."
)
@click.option(
    "--mass-tolerance",
    default=10.5,
    show_default=True,
    help="Provide the mass-tolerance."
)
@click.option(
    "--adducts",
    default=["M+H"],
    type=click.Choice(adduct_choices),
    multiple=True,
    show_default=True,
    show_choices=False,
    help="Provide the adducts."
)
@click.option(
    "--columns",
    default=compound_column_choices[:],
    type=click.Choice(compound_column_choices),
    multiple=True,
    show_default=True,
    show_choices=False,
    help="Provide the outputed columns."
)
@click.option(
    "--output-path",
    help="Provide the output path."
)
def compound(*args, **kwargs):

    if kwargs.pop("version"):
        print(__version__)
        exit(0)

    adducts = kwargs.pop("adducts")
    polarity = get_polarity(adducts)

    other_kwargs, compound_kwargs = build_kwargs(
        adducts=adducts,
        polarity=polarity,
        **kwargs
    )
    columns = other_kwargs["columns"]
    result = find_compound(**compound_kwargs)
    result = explode_compounds(
        result,
        with_annotations=any(map(
          columns.__contains__,
          annotation_column_choices
        ))
    )
    check_columns_in_result(result, columns)
    output_csv_result(
        result,
        columns,
        other_kwargs.get("output_path")
    )


def explode_compounds(result, with_annotations):
    if with_annotations:
        return [{
            "database": cpd.database,
            "metabolite_name": cpd.metabolite_name,
            "chemical_formula": cpd.chemical_formula,
            "hmdb_id": cpd.hmdb_id,
            "inchikey": cpd.inchikey,
            "compound_id": cpd.id,
            "adduct": annotation.name,
            "kendricks_mass": annotation.kendricks_mass,
            "kendricks_mass_defect": annotation.kendricks_mass_defect,
            "monisotopic_molecular_weight":
                annotation.monisotopic_molecular_weight,
            "nominal_mass": annotation.nominal_mass,
            "polarity": annotation.polarity,
            "annotation_id": annotation.id,
            }
            for cpd in result
            for annotation in cpd.annotations
        ]
    else:
        return [{
            "database": cpd.database,
            "metabolite_name": cpd.metabolite_name,
            "chemical_formula": cpd.chemical_formula,
            "hmdb_id": cpd.hmdb_id,
            "inchikey": cpd.inchikey,
            "compound_id": cpd.id,
            }
            for cpd in result
        ]


@cli.command(help="")
@click.option(
    "--id",
    type=int,
    help="Provide the wanted annotation's id."
)
@click.option(
    "--columns",
    default=annotation_column_choices[:],
    type=click.Choice(annotation_column_choices),
    multiple=True,
    show_default=True,
    show_choices=False,
    help="Provide the outputed columns."
)
@click.option(
    "--output-path",
    help="Provide the output path."
)
def annotation(*args, **kwargs):
    result = get_annotation(id=kwargs.pop("id"))
    result = [result]
    columns = kwargs["columns"]
    check_columns_in_result(result, columns)
    output_csv_result(
        result,
        columns,
        kwargs.get("output_path")
    )


def get_polarity(adducts):
    if any(map(positive_adducts.__contains__, adducts)):
        return "positive"
    if any(map(negative_adducts.__contains__, adducts)):
        return "negative"
    # polarity = []
    # if any(map(positive_adducts.__contains__, adducts)):
    #     polarity.append("positive")
    # if any(map(negative_adducts.__contains__, adducts)):
    #     polarity.append("negative")


def build_kwargs(**kwargs):
    for original, replacement in (
        ("database", "database_list"),
        ("polarity", "polarity_list"),
    ):
        if original in kwargs:
            kwargs[replacement] = kwargs.pop(original)
    other_kwargs = {
        other_arg: kwargs.pop(other_arg)
        for other_arg in ("columns", "output_path", "with_annotations")
        if other_arg in kwargs
    }
    return other_kwargs, kwargs


def check_columns_in_result(result, columns):
    if not result:
        return
    if not isinstance(result[0], dict):
        result = [item.to_dict() for item in result]
    keys = result[0].keys()
    missing = [
        column for column in columns
        if column not in keys
    ]
    if missing:
        if len(missing) == 1:
            raise ValueError(
                f"Could not find the column {missing[0]} in the results."
            )
        else:
            raise ValueError(
                "Could not find any of the columns "
                + ','.join(missing)
                + " in the results."
            )


def output_csv_result(result, columns, output_path, **csv_parameters):
    if not output_path:
        raise ValueError("Missing output path. Cannot output CSV results.")
    with open(output_path, mode="w", newline='') as output_file:
        writer = csv.writer(output_file, **csv_parameters)
        write_result(result, columns, writer)


def write_result(result, columns, writer):
    getters = list(map(operator.itemgetter, columns))
    writer.writerow(columns)
    writer.writerows(
        (getter(compound) for getter in getters)
        for compound in result
    )


if __name__ == "__main__":
    cli()
