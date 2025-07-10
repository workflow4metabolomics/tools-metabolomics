# Convex Hull Plot Generator

This R script generates convex hull plots for metabolite intensity data across different batches. It allows customization of the output file name, output directory, and whether to display individual data points.

## Requirements

Ensure you have R installed along with the necessary packages:

```r
install.packages(c("ggplot2", "optparse"))
```

## Usage

Run the script using the command line:

```r
Rscript plot_convex_hull.R -f data.csv -m "Metabolite X" -o "my_plot.png" -d "Plots/" -p FALSE
```

## Command-Line Options:

| Option         | Description                                        | Default                      |
|---------------|----------------------------------------------------|------------------------------|
| `-f`, `--file`       | Path to the input CSV file                    | `data.csv`                   |
| `-m`, `--metabolite` | Name of the metabolite for the plot title     | `"Metabolite X"`             |
| `-o`, `--output`     | Name of the output image file                 | `"plot_convex_hull.png"`     |
| `-d`, `--directory`  | Directory where the plot will be saved        | `"Graphe/"`                  |
| `-p`, `--points`     | To display individual data points (TRUE/FALSE) | `TRUE`               |
