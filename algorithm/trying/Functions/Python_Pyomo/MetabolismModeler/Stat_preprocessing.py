import pyomo.opt as po
import pyomo.environ as pe
import pandas as pd
import numpy as np
import openpyxl
import time
import json
import os
import re
from cobra import Model
from cobra.util import create_stoichiometric_matrix
from cobra.io import load_matlab_model, save_json_model, load_json_model
import matplotlib.pyplot as plt
from ast import literal_eval
from Network_analysis import visualize_network, identify_stoichiometric_matrix_of_sample_task_filtered_model
from typing import List, Tuple

def preprocess_and_save_sample_information_file(
    sample_information_file_location: str, expression_file_location: str
):

    sample_information_file = sample_information_file_location
    expression_file = expression_file_location

    if expression_file.endswith(".csv"):
        expression_df = pd.read_csv(expression_file, nrows=1)
    elif expression_file.endswith(".xlsx"):
        expression_df = pd.read_excel(expression_file, nrows=1)
    else:
        raise ValueError("Expression file is not in correct format")

    if sample_information_file.endswith(".csv"):
        sample_df = pd.read_csv(sample_information_file)
    elif sample_information_file.endswith(".xlsx"):
        sample_df = pd.read_excel(sample_information_file)
    else:
        raise ValueError("Sample information file is not in correct format")

    expression_df = expression_df.transpose()
    samples_in_expression = expression_df.index[1:]

    sample_df = sample_df[sample_df["ID"].isin(samples_in_expression)]
    samples_not_in_sample_df = set(samples_in_expression) - set(sample_df["ID"])
    indices_of_samples_not_in_sample_df = [
        expression_df.index.get_loc(x) + 1 for x in samples_not_in_sample_df
    ]
    print(
        f"Samples and indices not in sample_df: {list(zip(samples_not_in_sample_df, indices_of_samples_not_in_sample_df))}"
    )
    sample_df["Sample"] = sample_df["ID"].apply(
        lambda x: expression_df.index.get_loc(x)
    )
    sample_df = sample_df.reset_index(drop=True)
    sample_df = sample_df.drop(columns=["ID"])

    sample_information_file_location = os.path.join(
        os.path.dirname(sample_information_file), "Sample_information.xlsx"
    )
    sample_df.to_excel(sample_information_file_location)


def get_metabolites(x, model):
    try:
        metabolites = list(model.reactions.get_by_id(x["reaction_name"]).metabolites.keys())
        return metabolites
    except KeyError:
        return None


def add_metabolite_reaction_data_from_model(
        model: Model, df_ori: pd.DataFrame
) -> pd.DataFrame:
    df = pd.DataFrame(df_ori, copy=True)
    df["metabolites"] = df.apply(lambda x: get_metabolites(x, model), axis=1)

    print("added metabolites and reaction formulas to the dataframe")
    return df


def parse_json_and_create_dataframe(
    sample_information_file_location: str,
    json_file: str,
    task_number: int,
    cobra_model: Model,
    should_save: bool = True,
    method_name: str = "INIT_irrev",
) -> pd.DataFrame:
    # sample_information_file = "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\jelle\\main_statistics\\Sample_information.xlsx"
    print("reading sample information file")
    sample_information_file = sample_information_file_location
    sample_df = pd.read_excel(sample_information_file)

    ## assert all columns etc.
    # then prin
    print(
        f"loaded sample information file {sample_information_file}, all columns are present"
    )

    # json_file_location = "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\jelle\\test_runs\\INIT_new_consensus_model_task_250_small_example_sample_json\\INIT_new_consensus_model_task_250_small_example_sample_json_results.json"
    json_file_location = json_file
    with open(json_file_location, "r") as file:
        data = json.load(file)
    task_number = task_number

    ### add checks for the information/columns in the sample_information_file and json file and raise error if not found
    # assert
    list_ = ["sample", "task", "reaction_name", "flux", "expression"]
    for idx in range(len(list_)):
        print(list_[idx])

    for sample_number in data["sample_array"]:
        print(f"Processing sample {sample_number}")
        # current_sample = data["sample_array"][int(sample_number) - 1]
        whole_sample = data[method_name]["reactions_included_per_task"][
            str(sample_number)
        ][str(task_number)]
        reactions_in_sample_list = list(whole_sample.keys())
        for idx, reaction in enumerate(reactions_in_sample_list):
            list_for_each_reaction = list(
                data[method_name]["reactions_included_per_task"][str(sample_number)][
                    str(task_number)
                ][reaction]
            )
            # find any string with _f at the end and some part before it, then remove _f
            if re.search(r"_f", list_for_each_reaction[3]):
                list_for_each_reaction[3] = re.sub(r"_f", "", list_for_each_reaction[3])
            elif re.search(r"_b", list_for_each_reaction[3]):
                # TODO because we change flux here, and remove _b (also at _r) we need to change the sign of the flux
                # TODO but this also means that stoichiometric matrix of these values has to be changed so
                # TODO we should save the values per sample with negative flux so that when making filtered stoichiometric matrix
                # TODO we can change the values so that it works out correctly
                list_for_each_reaction[3] = re.sub(r"_b", "", list_for_each_reaction[3])
                list_for_each_reaction[1] = -1 * list_for_each_reaction[1]
            elif re.search(r"_r", list_for_each_reaction[3]):
                list_for_each_reaction[3] = re.sub(r"_r", "", list_for_each_reaction[3])
                list_for_each_reaction[1] = -1 * list_for_each_reaction[1]

            current_sample_information = sample_df.loc[
                sample_df["Sample"] == int(sample_number)
            ]
            if current_sample_information.empty:
                # create single row with nan values for all columns
                current_sample_information = pd.DataFrame(columns=sample_df.columns)
                current_sample_information = pd.concat(
                    [
                        current_sample_information,
                        pd.DataFrame(
                            [pd.Series([None] * len(sample_df.columns))],
                            columns=sample_df.columns,
                        ),
                    ]
                )

            filtered_current_dictionary = {
                "sample": [sample_number],
                "task": [task_number],
                "reaction_name": [list_for_each_reaction[3]],
                "flux": [list_for_each_reaction[1]],
                "expression": [list_for_each_reaction[2]],
                "expression_used_by_model": [list_for_each_reaction[6]],
                "healthy_status": [current_sample_information["etiology"].values[0]],
                "gender": [current_sample_information["gender"].values[0]],
                "race": [current_sample_information["race"].values[0]],
                "age": [current_sample_information["age"].values[0]],
                "weight": [current_sample_information["weight"].values[0]],
                "height": [current_sample_information["height"].values[0]],
                "diabetes": [current_sample_information["Diabetes"].values[0]],
                "hypertension": [current_sample_information["Hypertension"].values[0]],
                "LVEF": [current_sample_information["LVEF"].values[0]],
                "BMI": [current_sample_information["BMI"].values[0]],
                "lib_pool": [current_sample_information["LibPool"].values[0]],
                "reaction_formula": [list_for_each_reaction[5]],
            }
            ## add "total time" " total reactions"  " total expression carrying"  and " total expressionless reactions"

            if idx == 0 and sample_number == 1:
                df = pd.DataFrame.from_dict(filtered_current_dictionary)
            else:
                df_temp = pd.DataFrame(
                    filtered_current_dictionary,
                    columns=filtered_current_dictionary.keys(),
                )
                df = pd.concat([df, df_temp])
    df = df.reset_index(drop=True)
    df = df.astype('object')
    df = add_metabolite_reaction_data_from_model(cobra_model, df)

    if should_save:
        location_to_save = os.path.join(
            os.path.dirname(json_file_location), f"task_{task_number}_sample.xlsx"
        )
        try:
            df.to_excel(location_to_save)
        except PermissionError:
            print("Permission error")
        except Exception as e:
            print(e)

    return df


def get_reaction_formulas(x, model):
    try:
        reaction_formula = model.reactions.get_by_id(
            x["reaction_name"]
        ).build_reaction_string()
        return reaction_formula
    except KeyError:
        return None


## create simple function to plot total time


def create_graph_per_sample(
    df: pd.DataFrame, task_number: int, column_to_show: str
) -> plt.plot():

    # filter the df to only inlcude currently selected task

    # plot the column
    df[column_to_show].plot()

    # save plot as png
    return plt.plot()


### create sample_reaction tables (with fluxes, unqiuereactions and expression as input):
def create_sample_reaction_tables(
    df: pd.DataFrame, task_number: int, column_to_select="unique_reactions"
) -> pd.DataFrame:
    pass
    # special logic for unique reactions
    # otherwise just select the particular column


## create co uccrence table



def analyse_unique_reactions_for_sample_list(
    df_in: pd.DataFrame, task_number: int, samples_to_select: [int]
) -> dict:
    df = df_in.copy(deep=True)
    df = df[df["sample"].isin(samples_to_select)]
    df = df[df["task"] == task_number]
    sample_reaction_dict = (
        df.groupby("sample")["reaction_name"].apply(lambda x: set(x.unique())).to_dict()
    )
    all_reactions = set.intersection(*sample_reaction_dict.values())
    for sample, reactions in sample_reaction_dict.items():
        sample_reaction_dict[sample] = reactions.difference(all_reactions)
    return sample_reaction_dict


def analyse_unique_metabolites_for_sample_list(
    df_in: pd.DataFrame, task_number: int, samples_to_select: [int]
) -> dict:
    df = df_in.copy(deep=True)
    df = df[df["sample"].isin(samples_to_select)]
    df = df[df["task"] == task_number]
    pattern = re.compile(r'MAM\d{1,6}\[\w\]')
    def extract_metabolites(value):
        return pattern.findall(value)

    df['extracted_metabolites'] = df['metabolites'].apply(lambda x: extract_metabolites(x))

    sample_metabolite_dict = (
        df.groupby("sample")["extracted_metabolites"]
        .apply(lambda x: set(item for sublist in x for item in sublist))
        .to_dict()
    )

    all_metabolites = set.intersection(*sample_metabolite_dict.values())
    for sample, metabolites in sample_metabolite_dict.items():
        sample_metabolite_dict[sample] = metabolites.difference(all_metabolites)
    return sample_metabolite_dict


def analyse_core_reactions_common_to_all_samples(
    df_in: pd.DataFrame,
    task_number: int,
    samples_to_select: [int] = None,
) -> set:
    if samples_to_select is None:
        samples_to_select = list(df_in["sample"])
    df = df_in.copy(deep=True)
    df = df[df["sample"].isin(samples_to_select)]
    df = df[df["task"] == task_number]
    sample_reaction_dict = (
        df.groupby("sample")["reaction_name"].apply(lambda x: set(x.unique())).to_dict()
    )
    all_reactions = set.intersection(*sample_reaction_dict.values())
    return all_reactions

def analyse_core_metabolites_common_to_all_samples(
        df_in: pd.DataFrame,
        task_number: int,
        samples_to_select: [int] = None,
) -> set:
    if samples_to_select is None:
        samples_to_select = list(df_in["sample"])
    df = df_in.copy(deep=True)
    df = df[df["sample"].isin(samples_to_select)]
    df = df[df["task"] == task_number]
    pattern = re.compile(r'MAM\d{1,6}\[\w\]')
    def extract_metabolites(value):
        return pattern.findall(value)
    df['extracted_metabolites'] = df['metabolites'].apply(lambda x: extract_metabolites(x))
    sample_metabolite_dict = (
        df.groupby("sample")["extracted_metabolites"]
        .apply(lambda x: set(item for sublist in x for item in sublist))
        .to_dict()
    )
    all_reactions = set.intersection(*sample_metabolite_dict.values())
    return all_reactions


def identify_amount_of_temporary_exchange_reactions(
    df: pd.DataFrame, task_number: int
) -> list[str]:
    pass


def identify_expressionless_reactions(df: pd.DataFrame, task_number: int) -> list[str]:
    pass


if __name__ == "__main__":
    raw_expression_file = "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\jelle\\main_statistics\\geTMM_MAGNET_BC.csv"
    sample_information_file = "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\jelle\\main_statistics\\magnet_phenoData_complete.csv"
    json_file = "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\jelle\\test_runs\\ALLRUNS Regular fractional INIT local threshold\\INIT_new_consensus_model_task_1_all_samples_json\\INIT_new_consensus_model_task_1_all_samples_json_results.json"
    # model = load_matlab_model("E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\Base_files\\modelCellfie.mat\\")
    # save_json_model(model, "E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\Base_files\\modelCellfie.json", **{"allow_nan": True})
    cobra_model = load_json_model("E:\\Git\Metabolic_Task_Score\\Data\\Pyomo_python_files\\Base_files\\modelCellfie.json")
    task_number = 1
    # preprocess_and_save_sample_information_file(
    #     sample_information_file, raw_expression_file
    # )
    processed_sample_information_file = ("E:\\Git\Metabolic_Task_Score\\Data\\"
                                         "Pyomo_python_files\\jelle\\main_statistics\\"
                                         "Sample_information.xlsx")
    # df = parse_json_and_create_dataframe(
    #     processed_sample_information_file, json_file, task_number, cobra_model, True
    # )
    df = pd.read_excel(
        os.path.join(os.path.dirname(json_file), f"task_{task_number}_sample.xlsx")
    )
    df_in = df
    sample_number = 1
    task_number = 1
    # calculate_direct_path_using_least_connected_metabolites(model, df, 1, 1)
    # calculate_least_amount_of_edges_using_gurobi(model, df, 1, 1)

    (filtered_stoichiometric_matrix, reaction_index_dict,
            reaction_name_index_dict, metabolite_index_dict,
    metabolite_name_index_dict) = identify_stoichiometric_matrix_of_sample_task_filtered_model(
        cobra_model, df_in, task_number, sample_number
    )
    visualize_network(filtered_stoichiometric_matrix,
                      reaction_index_dict,
                      reaction_name_index_dict,
                        metabolite_index_dict,
                        metabolite_name_index_dict,
                      cobra_model)

    # unique_reactions = analyse_unique_reactions_for_sample_list(df, 1, [1, 10, 11, 30, 48, 316])
    # print("Unique reactions: ", unique_reactions)
    # unique_metabolites = analyse_unique_metabolites_for_sample_list(df, 1, [1, 10, 11, 30, 48, 316])
    # print("Unique metabolites: ", unique_metabolites)
    # common_reactions = analyse_core_reactions_common_to_all_samples(df, 1)
    # print("Common set of reactions: ", common_reactions)
    # common_metabolites = analyse_core_metabolites_common_to_all_samples(df, 1)
    # print("Common set of metabolites: ", common_metabolites)
    #
    # # convert metabolite or reaction to name from id


def create_solutions_of_other_samples_in_solution_matrix(
        solution_sets, #List[dict[str, List], # [{"sample_number: ["rxn1", "rxn2],}
        expression_matrix, # 322 * 13000
        cobra_model,
        expressionless_reaction_values_per_sample,
        objective_function_per_sample = "INIT_irrev",
) ->[List,List]:
    pass

    # solution_sets =
    #     {1: ""MAR04379", "MAR04358
    #          "MAR04363
    #      "MAR04368
    #          "MAR04391
    #      "MAR06412
    #          "MAR07747
    #      "MAR03829
    #          "MAR03957
    #      "MAR04139
    #          "MAR04145
    #      "MAR04152
    #          "MAR04209
    #      "MAR04410
    #          "MAR04589
    #      "MAR06413
    #          "MAR06414
    #      "MAR06914
    #          "MAR06916
    #      "MAR06918
    #          "MAR06921
    #      "MAR08746
    #          "MAR03825
    #      "MAR04855
    #          "MAR04898
    #      "MAR04926
    #          "MAR05043
    #      "MAR06328
    #          "MAR07638
    #      "MAR20069
    #          temporary_exchange_MAM01965[c]
    #      temporary_exchange_MAM02630[c]
    #          temporary_exchange_MAM02751[c]
    #      temporary_exchange_MAM01285[c]
    #          temporary_exchange_MAM02039[c]
    #      temporary_exchange_MAM01371[c]
    #          temporary_exchange_MAM02040[c]
    #      temporary_exchange_MAM01596[c]
    #          "MAR04365_r
    #      "MAR04373_r
    #          "MAR04375_r
    #      "MAR04381_r
    #          "MAR04002_r
    #      "MAR03827_r
    #          "MAR06409_r
    #      "MAR04141_r
    #          "MAR04408_r
    #      "MAR04458_r
    #          "MAR04652_r # >> find index of "MAR04652
    #      "MAR04852_r
    #          "MAR04922_r
    #      }
    #
    list_1 =["a", "b"]
    for value in list_1:
        print(value)
        # a
        # b
    for idx in range(len(list_1)):
        print(idx)
        # 0
        # 1
    for idx, value in enumerate(list_1):
        print(idx, value)
        # 0 a
        # 1 b
    # reactions = cobra.reactions
    # cobra_model.get_reaction_by_id("MAR04652")
    # indices = [idx for idx, rxn in enumerate(cobra_model.reactions) if rxn.id in solution_sets[1]]
    #expressions = expression_matrix[sample_number][indices]
    # E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\ALL RUNS Old runs\INIT_new_conensus_model_sample_1_to_501

    if objective_function_per_sample == "INIT_irrev":
        formula = min(rxn*1/gi)
    pass

    # make heatmap using plotly package

    # expression_matrix = cM.load_full_reversible_expression(main_data_folder, expression_type)
    # cobra_model = load_matlab_model(main_folder + "\\" + "modelCellfie.mat")

    # create a sample by sample matrix where for each row the values denote the objective if the column's solution
    # was calculated using the row's expression values

def create_jaccard_sample_x_sample_matrix():
    pass

def optimize_order_of_samples_in_jaccard_matrix():
    pass


def analyse_solutions_in_designated_folders():
    # return core reactions, unique reactions, unique metabolites, common reactions, common metabolites
    # return expressionless reactions, temporary exchange reactions
    # return amount of time to calculate each solution
    # also per group separately
    # option for comparing bottleneck with regular (lowest expression in both, alternatives between pairs)
    # comparing objectives for which it is necessary to recalculate the expressionless values \
    #               as now exchange reactions also get a different value. However task 10 shows some differences using the bottleneck!

    pass
def analyse_two_different_solutions():
    # return core reactions, unique reactions, unique metabolites, common reactions, common metabolites
    pass

    # test to check if solution of sample int he same sample gives the correct score
