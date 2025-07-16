"""
Module for converting PhysioFit output into .mflux mtf file for influx_si
"""
import logging
from pathlib import Path
from typing import Union, List, Tuple

import pandas as pd

_logger = logging.getLogger("root.physiofit2mtf")

def physiofit2mtf(
        data: pd.DataFrame,
) -> List[Tuple[Union[str, pd.DataFrame]]]:
    """
    Generate dataframes with the .mflux structure (1 dataframe per experiment/mflux file)

    :param physiofit_res: PhysioFit results file (path or dataframe)
    :return: list containing the dfs
    """


    _logger.debug(f"PhysioFit data before indexing:\n{data}")
    data = data.loc[(~data["parameter name"].str.contains("_M0")) & (~data["parameter name"].str.contains("_0"))]
    _logger.debug(f"PhysioFit data after indexing:\n{data}")
    # build .mflux
    _logger.info("Building .mflux file...")
    mflux_file = pd.DataFrame(columns=["Id", "Comment", "Flux", "Value", "SD"])
    mflux_file["Id"] = data["experiments"]
    mflux_file["Flux"] = data["parameter name"]
    mflux_file["Value"] = data["optimal"]
    mflux_file["SD"] = data["sd"]
    mflux_dfs = [
        (exp_name, mflux_file[mflux_file["Id"] == exp_name].copy())
        for exp_name in sorted(mflux_file["Id"].unique())
    ]
    _logger.debug("List of PhysioFit experiments and associated dataframes:")
    for (exp, df) in mflux_dfs:
        _logger.debug(f"Experiment {exp:}\n{df}")
    return mflux_dfs

def normalize_data(
        physiofit_data: pd.DataFrame,
        norm_value: Union[str, float],
        ignored_columns: Union[List[str], None] = ["experiments", "parameter name"],
) -> pd.DataFrame:
    """
    Normalize extracellular fluxes by a specific value

    :param physiofit_data: PhysioFit results file 
    :param norm_value: value to normalize by
    :return: normalized dataframe
    """
    if ignored_columns is None:
        ignored_columns = []
    _logger.debug(f"Running normalization with value: {norm_value} of type {type(norm_value)}")
    if not isinstance(physiofit_data, pd.DataFrame):
        raise TypeError("Physiofit_data must be a pandas DataFrame")
    
    if isinstance(norm_value, str):
        if norm_value not in physiofit_data["parameter name"].values:
            raise ValueError(f"Normalization parameter '{norm_value}' not found in PhysioFit data.")
        # We run normalization on each experiment separately
        for exp in physiofit_data["experiments"].unique():
            _logger.debug(f"Normalizing experiment {exp} by {norm_value}")
            # Get data for the current experiment
            exp_data = physiofit_data[physiofit_data["experiments"] == exp] 
            # Get the normalization value for the current experiment
            target_norm_value = abs(exp_data.loc[exp_data["parameter name"] == norm_value, "optimal"].values[0])
            _logger.debug(f"Experiment data before normalization:\n{exp_data}")
            if target_norm_value == 0:
                raise ValueError(f"Normalization value for {norm_value} in experiment {exp} is zero, cannot normalize.")
            _logger.debug(f"Target normalization value for {norm_value} in experiment {exp}: {target_norm_value}")
            # Normalize all columns except 'experiments' and 'parameter name' because they are not numerical
            for col in physiofit_data.columns:
                if col not in ignored_columns:
                    physiofit_data.loc[physiofit_data["experiments"] == exp, col] = exp_data[col] / target_norm_value
            _logger.debug(f"Experiment data after normalization:\n{physiofit_data[physiofit_data['experiments'] == exp]}")
    else:
        # Assume norm_value is a float or int
        if norm_value == 0:
            raise ValueError("Normalization value cannot be zero.")
        for col in physiofit_data.columns:
            if col not in ignored_columns:
                physiofit_data[col] = physiofit_data[col] / norm_value
        _logger.debug(f"Data after normalization by constant {norm_value}:\n{physiofit_data}")

    return physiofit_data
