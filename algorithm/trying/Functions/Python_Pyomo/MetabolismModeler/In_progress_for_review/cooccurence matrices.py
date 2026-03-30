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
from typing import List, Tuple
import glob

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

def select_sample_and_task_from_dataframe(
        df_in: pd.DataFrame, task_number: int, sample_number: int
) -> pd.DataFrame:
    df = df_in.copy(deep=True)
    df = df[df["sample"] == sample_number]
    df = df[df["task"] == task_number]
    return df

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




#[List,List]:
#) ->

#def co_oocurrence_matrix(
#    solution_sets: dict[int, List[str]],
#):

def create_solutions_of_other_samples_in_solution_matrix(
        solution_sets: dict[int, List[str]],
        expression_matrix: [List[float]],  # size 2d reaction * sample
        cobra_model: Model,  # contains reactions in same order as expression_matrix
        objective_function_per_sample: str,
        expressionless_reaction_value_per_sample: dict[int, float],
):
    union_of_reactions = []
    def unique_list_keep_order(lst):
        unique_list = []
        for item in lst:
            if item not in unique_list:
                unique_list.append(item)
        return unique_list

    for sample_number in solution_sets.keys():
        list_of_reactions = solution_sets[sample_number]
        for reaction_idx, reaction in enumerate(list_of_reactions):
            if re.search(r"_f", reaction):
                list_of_reactions[reaction_idx] = re.sub(r"_f", "", reaction)
            elif re.search(r"_b", list_of_reactions[reaction_idx]):
                list_of_reactions[reaction_idx] = re.sub(r"_b", "", list_of_reactions[reaction_idx])
            elif re.search(r"_r", list_of_reactions[reaction_idx]):
                list_of_reactions[reaction_idx] = re.sub(r"_r", "", list_of_reactions[reaction_idx])
        union_of_reactions.extend(list_of_reactions)
    union_of_reactions = unique_list_keep_order(union_of_reactions)
    samples_in_analysis = [sample for sample in solution_sets.keys()] #[3, 4, 5]
    # {3: [rxn1, rxn2], 4: [rxn1, rxn3], 5: [rxn2, rxn3]}
    solution_sets = {idx+1: list_ for idx, (sample, list_) in enumerate(solution_sets.items())}

    # find the expression value for the reaction
    reaction_indices = [idx for idx, reaction in enumerate(cobra_model.reactions) if reaction.id in union_of_reactions]
    union_of_reactions_without_temporary_exchange = [reaction for reaction in union_of_reactions if not reaction.startswith("temporary_exchange")]
    reaction_indices_tuple = [(reaction, idx) for reaction, idx in zip(union_of_reactions_without_temporary_exchange, reaction_indices)]
    # ([rxn1, 5000), (rxn2, 6000), (rxn3, 2304)]
    reaction_indices_tuple = sorted(reaction_indices_tuple, key=lambda x: x[1])
    sorted_reactions_dictionary = {val[0]: val[1] for val in reaction_indices_tuple}

    # sorting tuple based on index -- we would have indexes in chronological order + reactions in a tuple
    # sorted_reactions_dictionary = {tup[0]: idx for idx, tup in enumerate(reaction_indices_tuple)}

    matrix_with_values = np.zeros((len(solution_sets.keys()), len(solution_sets.keys())))
    for row_sample_number in range(len(solution_sets.keys())):
        # expression_of_sample = expression_matrix[:, row_sample_number]
        expression_of_sample = [value for sublist_ in expression_matrix for idx,
                                     value in enumerate(sublist_) if idx == samples_in_analysis[row_sample_number]-1]
        expressionless_reaction_value = expressionless_reaction_value_per_sample[row_sample_number+1]
        expression_of_sample = np.where(np.logical_or(expression_of_sample == -1, np.isnan(expression_of_sample)),
                                        expressionless_reaction_value,
                                        expression_of_sample)
        expression_of_sample = np.where(np.logical_or(expression_of_sample == -1, np.isnan(expression_of_sample)),
                                        expressionless_reaction_value,
                                        expression_of_sample)
        for col_sample_number in range(len(solution_sets.keys())):
            reactions_in_solution = solution_sets[col_sample_number+1]
            temporary_exchange_reactions = [reaction for reaction in reactions_in_solution if reaction.startswith("temporary_exchange")]
            reactions_in_solution = [reaction for reaction in reactions_in_solution if reaction in sorted_reactions_dictionary.keys()]

            expression_of_reactions = [expression_of_sample[sorted_reactions_dictionary[reaction]]
                                       for reaction in reactions_in_solution if reaction in sorted_reactions_dictionary.keys()]
            amount_to_add_to_sum = len(temporary_exchange_reactions) * (1/expressionless_reaction_value)
            objective_value = sum([1/expression for expression in expression_of_reactions]) + amount_to_add_to_sum
            print(f"Objective value for sample {row_sample_number} and sample {col_sample_number} is {objective_value}")
            matrix_with_values[col_sample_number, row_sample_number] = objective_value

    print(matrix_with_values)

    # TODO: create a heatmap using plotly package
    # TODO: co-occurence matrix
    # TODO: make heatmap of co-occurence matrix
    # create a dash app that does the following:
        # 1. allows user to upload filtered_excel solution files
        # 2. identify the task number, and sample number
        # 3. show the core reactions (common to all samples)
        # 4. and show the user the variable reactions per sample
        # 5. show the user a table with the expresssion values of the union of reactions of the
            # uploaded samples of each sample


    # sample 1 = [rxn1, rxn2, rxn3]
    # sample 2 = [rxn1, rxn3, rxn4]
    # common = [rxn1, rxn3]
    # variable_in_sample1 = [rxn2]
    # variable_in_sample2 = [rxn4]

    # union = [rxn1, rxn2, rxn3, rxn4]




    # make heatmap using plotly package

    # expression_matrix = cM.load_full_reversible_expression(main_data_folder, expression_type)
    # cobra_model = load_matlab_model(main_folder + "\\" + "modelCellfie.mat")

    # create a sample by sample matrix where for each row the values denote the objective if the column's solution
    # was calculated using the row's expression values

def create_jaccard_sample_x_sample_matrix():
    pass

def optimize_order_of_samples_in_jaccard_matrix():
    pass


    # test to check if solution of sample int he same sample gives the correct score
if __name__ == "__main__":
    path = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Base_files"
    expression_file_name = "full_expression_reversible_no_threshold.json"
    expression = json.load(open(os.path.join(path, expression_file_name)))
    cobra_model = load_json_model(os.path.join(path, "modelCellfie.json"))
    objective_string = "INIT_irrev"
    solution_sets = {3: [
        "MAR04388_f",
        "MAR04379",
        "MAR04358",
        "MAR04363_f",
        "MAR04368_f",
        "MAR04391_f",
        "MAR06412",
        "MAR07747",
        "MAR03957",
        "MAR04145",
        "MAR04152_f",
        "MAR04209",
        "MAR04410_f",
        "MAR04589_f",
        "MAR06413",
        "MAR06414_f",
        "MAR06914",
        "MAR06916",
        "MAR06918",
        "MAR06921",
        "MAR08746",
        "MAR04855_f",
        "MAR04898_f",
        "MAR05043",
        "MAR06328_f",
        "MAR07638",
        "MAR01485",
        "MAR20069",
        "temporary_exchange_MAM01965[c]",
        "temporary_exchange_MAM02630[c]",
        "temporary_exchange_MAM02751[c]",
        "temporary_exchange_MAM01285[c]",
        "temporary_exchange_MAM02039[c]",
        "temporary_exchange_MAM01371[c]",
        "temporary_exchange_MAM02040[c]",
        "temporary_exchange_MAM01596[c]",
        "MAR04365_r",
        "MAR04373_r",
        "MAR04375_r",
        "MAR04381_r",
        "MAR04280_r",
        "MAR04002_r",
        "MAR06409_r",
        "MAR04141_r",
        "MAR04408_r",
        "MAR04458_r",
        "MAR04652_r",
        "MAR04922_r",
    ], 2: ["MAR04922_r"]}
    expressionless_reaction_dict = {i+1: 1.1481 for i in range(len(solution_sets.keys()))}

    create_solutions_of_other_samples_in_solution_matrix(
        solution_sets,
        expression,
        cobra_model,
        objective_string,
        expressionless_reaction_dict
    )
