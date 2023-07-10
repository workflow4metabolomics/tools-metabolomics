#!/usr/bin/env python3

import csv
import itertools
import os

import click

import plotly.express
import plotly.graph_objects

__version__ = "1.0.0"


@click.group()
def cli():
    pass


@cli.command(help="")
@click.option(
    "--version",
    is_flag=True,
    default=False,
)
@click.option(
    "--input",
    default="./test.csv",
    help="Provide the mz-ratio."
)
@click.option(
    "--output",
    default="./test.html",
    help="Provide the database."
)
@click.option(
    "--x-column",
    default=["nominal_mass"],
    multiple=True,
    help="Provide the column names for the X axis.",
)
@click.option(
    "--y-column",
    default=["kendricks_mass_defect"],
    multiple=True,
    help="Provide the column names for the Y axis.",
)
@click.option(
    "--annotation-column",
    multiple=True,
    default=[
        "metabolite_name",
        "chemical_formula",
    ],
    help="Provide the columns name for the annotation."
)
def plot(*args, **kwargs):

    if kwargs.pop("version"):
        print(__version__)
        exit(0)

    input_path = kwargs.pop("input")
    data = read_input(input_path, kwargs)
    fig = build_fig(*data)
    build_html_plot(fig, kwargs.get("output"))


def read_input(path: str, kwargs: {}):
    if not os.path.exists(path):
        raise ValueError(f"The path '{path}' does not exist.")
    sep = detect_sep(path)
    with open(path) as csv_file:
        line_generator = csv.reader(csv_file, delimiter=sep)
        first_line = next(line_generator)
        all_lines = list(line_generator)
        hover_names = (
            "metabolite_name",
            "chemical_formula",
        )
        annotation_indexes = get_index_of(first_line, hover_names)
        (
            x_index,
            y_index,
            x_column,
            y_column,
        ) = get_indexes_names(
            first_line,
            list(kwargs.get("x_column")),
            list(kwargs.get("y_column")),
        )
        x_lists = [[] for i in range(len(x_index))]
        y_lists = [[] for i in range(len(y_index))]
        x_column = list(map(first_line.__getitem__, x_index))
        y_column = list(map(first_line.__getitem__, y_index))
        trace_names = [
            f"f({x_column[i]}) = {y_column[i]}"
            for i in range(len(x_index))
        ]
        hover_names = kwargs["annotation_column"]
        annotation_indexes = [
            get_index_of(first_line, column)[0]
            for column in hover_names
        ]
        hover_names = list(map(first_line.__getitem__, annotation_indexes))
        annotations = list()
        for line in all_lines:
            for i in range(len(x_index)):
                x_lists[i].append(float(line[x_index[i]]))
                y_lists[i].append(float(line[y_index[i]]))
            annotations.append("<br>".join(
                f"{hover_names[hover_index]}: {line[index]}"
                for hover_index, index in enumerate(annotation_indexes)
            ))
    return x_lists, y_lists, annotations, trace_names


def get_indexes_names(first_line, x_column, y_column):
    x_column, y_column = map(list, zip(*itertools.product(x_column, y_column)))
    x_index = get_index_of(first_line, x_column)
    y_index = get_index_of(first_line, y_column)
    for i in range(len(x_index))[::-1]:
        if x_index[i] == y_index[i]:
            del x_index[i], x_column[i], y_index[i], y_column[i],
    return (
        x_index,
        y_index,
        x_column,
        y_column,
    )


def get_index_of(first_line, column):
    if isinstance(column, (tuple, list)):
        return [get_index_of(first_line, x)[0] for x in list(column)]
    try:
        return [int(column) - 1]
    except ValueError:
        return [first_line.index(column)]


def build_fig(x_lists, y_lists, annotations, trace_names):
    fig = plotly.express.scatter()
    for i in range(len(x_lists)):
        fig.add_trace(
            plotly.graph_objects.Scatter(
                name=trace_names[i],
                x=x_lists[i],
                y=y_lists[i],
                hovertext=annotations,
                mode="markers",
            )
        )
    return fig


def detect_sep(tabular_file: str) -> str:
    with open(tabular_file, "r") as file:
        first_line = file.readline()
    if len(first_line.split(',')) > len(first_line.split('\t')):
        return ','
    return '\t'


def build_html_plot(fig, output: str):
    return plotly.offline.plot(
        fig,
        filename=output,
        auto_open=False,
    )


if __name__ == "__main__":
    cli()
