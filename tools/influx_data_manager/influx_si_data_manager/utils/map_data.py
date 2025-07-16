import pandas as pd
import logging

CONNECTOR_NAME = "influx_si_data_connector"

logger = logging.getLogger("root.influx_si_data_manager.utils.map_data")

def map_data(mapping_file, data, from_tool):
    """
    Replace cell values in dataset by following indications in mapping file
    :param mapping_file: Path to the mapping file
    :param data: dataset to modify
    :param from_tool: tool name to filter mapping entries
    :return: Mapped data
    """
    # Read the mapping file and filter for relevant entries
    mapping = pd.read_csv(mapping_file, sep="\t")
    filtered_mapping = mapping[(mapping["Connector"] == CONNECTOR_NAME) & 
                              (mapping["From_tool"] == from_tool)]
    
    # Replace values in each column
    for column in filtered_mapping["Column"].unique():
        print(f"Data to modify: {data}")
        logger.debug(f"Replacing values in column: {column}")
        logger.debug(f"Data to modify: {data}")
        replacements = filtered_mapping[filtered_mapping["Column"] == column]
        replace_dict = dict(zip(replacements["Replace"], replacements["By"]))
        logger.debug(f"Replacement dictionary: {replace_dict}")
        print(replace_dict)
        logger.debug(f"Column values before replacement: {data[column]}")
        data[column] = data[column].replace(replace_dict)
        logger.debug(f"Column values after replacement: {data[column]}")
    
    return data
            
            

if __name__ == "__main__":
    
    #isocor_data = pd.read_csv("../test_data/isocor_results.tabular", sep="\t")
    #print(isocor_data)
    physiofit_data = pd.read_csv(r"/home/llegregam/Documents/tools_w4m/packages/influx_data_manager/influx_si_data_manager/test_data/summary.csv", sep=",")
    print(physiofit_data)
    replaced = map_data(r"/home/llegregam/Documents/tools_w4m/packages/influx_data_manager/influx_si_data_manager/test_data/mapping.txt", data=physiofit_data, from_tool="physiofit")
    print(replaced)