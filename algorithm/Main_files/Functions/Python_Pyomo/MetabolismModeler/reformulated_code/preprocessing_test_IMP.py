'''
This module contains the preprocessing functions that will take the following inputs:
    - A metabolic model starting with model_ (SBML compatible or cobrapy model), in json, matlab, or SBML format

    - A gene expression dataset (in the form of a json or excel file) with the following characteristics:
        - name of dataset starts with "data_"
        - Gene ID as rows of the first column
        - Expression values of each gene in from the second row and column onwards
        - Sample IDs as columns of the first row

    - A sample information dataset (in the form of a json or excel file) with the following characteristics:
        - name of dataset starts with "sample_info"
        - Sample ID as rows of the first column
        - Any other characteristics of the samples in the following columns (e.g. tissue, age, condition).

    - A tasklist (in the form of a json or excel file) with the following characteristics:
        - "ID" as the first column (starting from 1 and positive integers only)
        - "DESCRIPTION" as the second column with task names/descriptions
        - "IN" as the fourth column with the input metabolites in MAM0xxxx notation
        - "IN LB" as the fifth column with the lower bounds of the input metabolites
        - "IN UB" as the sixth column with the upper bounds of the input metabolites
        - "OUT" as the seventh column with the output metabolites in MAM0xxxx notation
        - "OUT LB" as the eighth column with the lower bounds of the output metabolites
        - "OUT UB" as the ninth column with the upper bounds of the output metabolites
        - Any columns starting with "CLASSIFIER_" are used to classify the tasks

The model should be in a directory of its name inside the "models" directory.
The gene expression and sample information should be in a directory of
        their name inside the "expression_datasets" directory.
The tasklist should be in a directory of its name inside the "tasklists" directory.

The preprocessing functions will then be called by running.py by providing a methods_dict with
__global_options from which the following files are created (if not already present):
    - A "combinations" directory based on the expression data, sample information, and tasklist that contains:
        - A closed base model (with the same name as the model) in the "models" directory
        - A irreversible version of the base model in the "models" directory
        - Expression data in the form of a feather file parsed to the base models:
            - The first column, second row onwards contains the reaction names of the model
            - The first row, second column onwards contains the sample IDs
            - The dataset has the same amount (+1) of rows as the model has reactions, with each value
                    being the expression value of the gene that codes for the reaction of the sample
                    in its corresponding column.
            - A reversible dataset
            - An irreversible dataset
            - Specific naming convention for the different methods of creating the expression data:
                - {transformation}__{thresholding}__{GPR_method}__{promiscuity_method}__{model_type}__name
                - transformation:   NT (no transformation),
                                    L2T (log2 transformation),
                                    L10T (log10 transformation),
                                    LNT (ln transformation),
                                    ZTL (z-score transformation local) only using data in sample,
                                    ZTG (z-score transformation global) using all data in expression

                - thresholding:     NT (no thresholding),
                                    GTPL_V (global thresholding percentile local), V is the percentile
                                            for the threshold taken from the local (this sample) data
                                    GTPG_V (global thresholding percentile global), V is the percentile
                                            for the threshold taken from the global data
                                    GTV_V (global thresholding value), V is the value for the threshold
                                    LTM (local thresholding mean),
                                    LTMMMP_V_V (local thresholding min max mean with percentile),
                                            V is the percentile low and high
                                    LTMMMV_V_V (local thresholding min max mean with value),
                                            V is the value low and high



                - GPR_method:
                                    min1Sum (taking the sum of OR values) min1 indicates to only take the min value as contributor
                                    min2Sum (taking the sum of OR values) min2 indicates that all values of the AND min rule are contributors
                                    min1Max (taking the maximum of OR values),
                                    min1Mean (taking the mean of OR values),
                                    min1Median (taking the median of OR values),
                                    trimmeanSum_V (taking the mean of AND values after trimming V%),
                                                and then taking the sum of OR values,
                                    trimmeanMax_V (taking the mean of AND values after trimming V%),
                                                and then taking the maximum of OR values
                                    trimmin1Sum_V (taking the min of AND values after trimming V%),
                                                and then taking the sum of OR values
                                    trimmin2Max_V (taking the min of AND values after trimming V%),
                                                and then taking the maximum of OR values

                - promiscuity_method:   Y (yes, a gene being the determining factor for multiple reactions
                                                is reduced by 1/n)
                                        N (no)

                - model_type:       rev (reversible),
                                    irr (irreversible)
            - A list  linking sample IDs to corresponding numbers as .json
            - A file with information on the added reactions and bounds within the model as .json
            - A metadata/log file with information on the creation of the combination files and whenever anything is ran
            - (Optional only for reversible) A rev_to_irrev file with information on the added reactions and bounds
                    within the model as .json for the irreversible model
            - (Optional only for irreversible) A list linking each reaction's index with its
                    corresponding irreversible reaction index as .json

Initially created on 2024-09-20 by Jelle Bonthuis (MaCSBio)
Updated on {date} by {author}:
    - {list of changes}
'''
version = "1.0.0"
# todo add metadata file that includes version and date of creation
# TODO add major gene per reaction contribution to CSV

# Imports
import os
import pandas as pd
from typing import List, Dict, Tuple, Union, Any
from copy import deepcopy
from logging import warning, error, info, debug
from cobra import Model, Reaction
from cobra.io import (save_json_model,
                      save_matlab_model,
                      write_sbml_model,
                      load_json_model,
                      load_matlab_model,
                      read_sbml_model,
                      )
import time
import cProfile  # TODO for testing only
from json import load, dump
from scipy.stats import zscore, trim_mean
import numpy as np
import re


'''
* testing see if code works > e.g. if folder empty then create files, does this do it with empty folder, and then tear-down (jargon)
* pytest 
* put files in separate folder
* sonarqube cyclomatic complexity and code coverage
* add Code2Flow to see flow of code or pycallgraph
* pythontutor for homework exercises studetns
'''

# Functions
class FileChecker:
    @staticmethod
    def read_model_file_in_folder(
            model_name: str = None,
            model_folder: str = None,
            model_file_name: str = None,
    ) -> Model:
        """
        Read the model file in the folder.

        :param model_name:
            The name of the model file.
        :param model_folder:
            The path to the folder containing the model file.
        :param model_file_name:
            The name of the model file.
        :return:
            The cobrapy model.
        """
        if model_file_name is not None:
            try:
                if model_file_name.endswith(".json"):
                    # return load_json_model(os.path.join(model_folder, model_file_name))
                    return load_json_model(model_file_name)
                elif model_file_name.endswith(".mat"):
                    return load_matlab_model(model_file_name)
                    # return load_matlab_model(os.path.join(model_folder, model_file_name))
                elif model_file_name.endswith(".xml"):
                    return read_sbml_model(model_file_name)
                    # return read_sbml_model(os.path.join(model_folder, model_file_name))
            except Exception:
                warning(f"An error occurred while reading the model file {model_file_name}.")

        if (model_folder is None and
                model_name is None):
            raise ValueError("Either model_name or model_folder should be provided.")
        # print("DEBUG: model_folder =", model_folder)
        # print("DEBUG: model_name =", model_name)
        # print("DEBUG: files in folder =", os.listdir(model_folder))

        if model_folder is None:
            model_folder = os.path.join("models", model_name)
        if model_name is None:
            model_files = [file for file in os.listdir(model_folder) if file.lower().startswith("model")]
            print(model_files)
        else:
            model_files = [file for file in os.listdir(model_folder) if model_name in file]
        if not model_files:
            raise FileNotFoundError(f"No model file found in the folder {model_folder} with the name {model_name}.")
        for model_file in model_files:
            try:
                if model_file.endswith(".json"):
                    return load_json_model(os.path.join(model_folder, model_file))
                elif model_file.endswith(".mat"):
                    return load_matlab_model(os.path.join(model_folder, model_file))
                elif model_file.endswith(".xml"):
                    return read_sbml_model(os.path.join(model_folder, model_file))
            except Exception as e:
                # pass
                print(f"DEBUG: Failed to load {model_file} with error: {e}")
        raise FileNotFoundError(f"No model file found in the folder {model_folder} with the name {model_name}.")

    @staticmethod
    def raise_warning_if_origin_files_not_present(
            file_or_folder_path: str,
    ):
        """
        Raise a warning if the file or folder path is not present.

        :param file_or_folder_path: str
            The path to the file or folder.
        """
        if not os.path.exists(file_or_folder_path):
            warning(f"The combination folder is present and populated, but the origin file or folder {file_or_folder_path} does not exist.\n"
                    f"While this is not problematic, it might indicate that the origin files were moved or deleted.")

    @staticmethod
    def raise_error_if_incorrect_folders(
            file_or_folder_path: str,
    ):
        """
        Raise an error if the file or folder path is not in the correct directory.

        :param file_or_folder_path: str
            The path to the file or folder.
        """
        if not os.path.exists(file_or_folder_path) and not os.path.isdir(file_or_folder_path):
            if os.path.exists(file_or_folder_path):
                raise FileNotFoundError(f"The file or folder {file_or_folder_path} does not exist.")
            elif not os.path.isdir(file_or_folder_path):
                raise NotADirectoryError(f"The path {file_or_folder_path} is not a directory.")

    @staticmethod
    def check_if_combination_files_exist(
            combination_folder: str,
            global_options_dict: dict,
    ) -> (dict, dict):
        check_if_all_correct_dict = {
            "expression_data": True,
            "model": True,
            "tasklist": True,
            "rev_to_irrev_dict": True,
        }
        selected_model = global_options_dict["__global_options"]["model_version"]
        selected_tasklist = global_options_dict["__global_options"]["tasklist"]
        expression_name_by_options = (f"{global_options_dict['__global_options']['transformation']}__"
                                      f"{global_options_dict['__global_options']['thresholding']}__"
                                      f"{global_options_dict['__global_options']['GPR_method']}__"
                                      f"{global_options_dict['__global_options']['promiscuity_method']}__"
                                      f"{global_options_dict['__global_options']['model_type']}__"
                                      f"reaction_expression_dataset.feather")
        expression_name_by_options = os.path.join(combination_folder, expression_name_by_options)
        if not os.path.exists(expression_name_by_options):
            check_if_all_correct_dict["expression_data"] = False
        irreverible_or_reversible = "irreversible" if global_options_dict["__global_options"]["model_type"] == "irr" else "reversible"
        model_name_by_options = (
            f"model_{selected_model}__{irreverible_or_reversible}"
            f".json"
        )
        model_name_by_options = os.path.join(combination_folder, model_name_by_options)
        if not os.path.exists(model_name_by_options):
            check_if_all_correct_dict["model"] = False
        tasklist_name_by_options = (
            f"task_structure_{selected_tasklist}.json"
        )
        tasklist_name_by_options = os.path.join(combination_folder, tasklist_name_by_options)
        rev_to_irrev_dict_name_by_options = (
            f"{selected_model}__"
            f"rev_to_irrev_dict.json"
        )
        if global_options_dict["__global_options"]["model_type"] == "irr":
            rev_to_irrev_dict_name_by_options = os.path.join(combination_folder, rev_to_irrev_dict_name_by_options)
            if not os.path.exists(rev_to_irrev_dict_name_by_options):
                check_if_all_correct_dict["rev_to_irrev_dict"] = False
        else:
            check_if_all_correct_dict.pop("rev_to_irrev_dict")

        names_dict_by_options = {
            "expression_data": expression_name_by_options,
            "model": model_name_by_options,
            "tasklist": tasklist_name_by_options,
        }
        if global_options_dict["__global_options"]["model_type"] == "irr":
            names_dict_by_options["rev_to_irrev_dict"] = rev_to_irrev_dict_name_by_options
        return (
            check_if_all_correct_dict,
            names_dict_by_options
        )

    @staticmethod
    def check_if_expression_data_is_correct(
            expression_data_file: str,
            alternative_reaction_notation: bool = False,
            already_loaded_expression_data: pd.DataFrame = None,
    ) -> (bool, (pd.DataFrame | None)):
        """
        Check if the expression data file is correct.
            - Is json readable
            - First column, second row onwards contains the reaction names of the model
            - First row, second column onwards contains the sample IDs
        :param expression_data_file: str
            The path to the expression data file.
        :param alternative_reaction_notation: bool
            Whether the reaction names are in a different notation than MARxxxx.
        :param already_loaded_expression_data: pd.DataFrame
            The already loaded expression data.
        :return: (is_correct: bool, expression_data: pd.DataFrame)
            is_correct is True if the expression data file is correct.
        """

        def check_if_string_or_integer(column):
            string_count = column.astype(str).str.contains("[a-zA-Z]").sum()
            integer_count = column.astype(str).str.contains("[0-9]").sum()
            return string_count >= 2 or integer_count >= 2

        is_correct = True
        if already_loaded_expression_data is None:
            try:
                if expression_data_file.endswith(".json"):
                    expression_data = pd.read_json(expression_data_file)
                elif expression_data_file.endswith(".xlsx"):
                    expression_data = pd.read_excel(expression_data_file)
                elif expression_data_file.endswith(".csv"):
                    expression_data = pd.read_csv(expression_data_file)
                elif expression_data_file.endswith(".feather"):
                    with open(expression_data_file, "rb") as f:
                        expression_data = pd.read_feather(f)
            except Exception as e:
                warning(f"An error occurred while reading the expression data file {expression_data_file}.\n"
                        f"Error message: {e}")
                return False, None
        else:
            expression_data = already_loaded_expression_data
        if alternative_reaction_notation:
            if not check_if_string_or_integer(expression_data.index):
                warning(f"The expression data file {expression_data_file} "
                        f"does not contain any likely reaction names.\n"
                        f"Consider characters or integers as reaction names of your model.\n"
                        f"This might indicate that extraction of reaction names from model failed.")
                is_correct = False
            raise NotImplementedError("Alternative reaction notation is not yet implemented.")
        elif not any(expression_data.index.str.contains("MAR")):
            warning(f"The expression data file {expression_data_file} "
                    f"does not contain any reaction names in MARxxxx notation.\n"
                    f"If this was intentional, please set alternative_reaction_notation to True."
                    f"in the global options.")
            is_correct = False
        if not check_if_string_or_integer(expression_data.columns):
            warning(f"The expression data file {expression_data_file} "
                    f"does not contain any likely sample IDs.\n"
                    f"Consider characters or integers as sample IDs.\n"
                    f"This might indicate that extraction of sample IDs from the expression data failed.")
            is_correct = False
        if not is_correct:
            info(f"Since one or more errors were found in the expression data file {expression_data_file}, "
                 f"the file will be renamed to {expression_data_file}_incorrect.json and "
                 f"a new expression data file will be created.")
            os.rename(expression_data_file, expression_data_file.split(".feather")[0]+ "_incorrect.feather")
            return False, None
        return True, expression_data

    @staticmethod
    def check_if_combination_model_is_correct(
            model_file_name: str,
            already_loaded_model: Model = None,
    ) -> (bool, (Model | None)):
        """
        Check if the model file is correct.
            - Is json readable
            - Is a cobrapy model
            - Is a closed model (all exchange reactions are closed)
            - All reactions have bounds
            - Is viable (LP minimization is possible)
        :param model_file_name: str
            The path to the model file.
        :param already_loaded_model: Model
            The already loaded model.
        :return: (is_correct: bool, cobra_model: Model)
            is_correct is True if the model file is correct.
        """
        if already_loaded_model is None:
            cobra_model = FileChecker.read_model_file_in_folder(None,
                                                                None,
                                                                model_file_name
                                                                )
        else:
            cobra_model = already_loaded_model
        is_correct = True
        # Check if the model is closed
        if not all(reaction.bounds == (0.0, 0.0) for reaction in cobra_model.boundary):
            warning(f"The model file {model_file_name} contains exchange reactions without bounds.")
            is_correct = False
        # Check if all reactions have bounds
        if not all(reaction.upper_bound is not None and reaction.lower_bound is not None for reaction in cobra_model.reactions):
            warning(f"The model file {model_file_name} contains reactions without bounds.")
            is_correct = False
        if any( rxn.lower_bound != 0 and rxn.upper_bound != 0 and ((rxn.upper_bound > 0 and rxn.lower_bound > 0) or (rxn.upper_bound < 0 and rxn.lower_bound < 0)) for rxn in cobra_model.reactions):
            warning(f"The model file {model_file_name} contains reactions with bounds that are forced open.")
            is_correct
        # Check if the model is viable
        try:
            cobra_model.slim_optimize()
        except Exception as e:
            warning(f"The model file {model_file_name} is not viable.\n"
                    f"Error message: {e}")
            is_correct = False
        if not is_correct:
            info(f"Since one or more errors were found in the model file {model_file_name}, "
                 f"the file will be renamed to {model_file_name}_incorrect.json and "
                 f"a new model file will be created.")
            os.rename(model_file_name, model_file_name.split(".json")[0] + "_incorrect.json")
            return False, None
        return True, cobra_model

    @staticmethod
    def check_if_combination_expression_and_model_are_congruent(
            cobra_model: Model,
            expression_data: pd.DataFrame,
    ) -> bool:
        """
        Check if the amount of reactions in the expression data is the same as in the model.

        :param cobra_model: Model
            The cobrapy model.
        :param expression_data: pd.DataFrame
            The expression data.
        """
        if (len(cobra_model.reactions) != len(expression_data)):
            warning(f"The amount of reactions in the model ({len(cobra_model.reactions)}) "
                    f"does not match the amount of reactions in the expression data ({len(expression_data) - 1})")
            return False
        return True

    @staticmethod
    def check_if_irrev_dict_is_correct(
            rev_to_irrev_dict_name_by_options: str,
            cobra_model: Model,
    ) -> (bool, (dict | None)):
        """
        Check if the irreversible dictionary is correct.
            - Is json readable
            - Contains integer keys
            - Has the same amount of keys as the amount of reactions in the model
        :param rev_to_irrev_dict_name_by_options: str
            The path to the irreversible dictionary.
        :param cobra_model: Model
            The cobrapy model.
        :return: (is_correct: bool, rev_to_irrev_dict: dict)
            is_correct is True if the irreversible dictionary is correct.
        """
        is_correct = True
        try:
            rev_to_irrev_dict = load(open(rev_to_irrev_dict_name_by_options))
        except Exception as e:
            warning(f"An error occurred while reading the irreversible dictionary {rev_to_irrev_dict_name_by_options}.\n"
                    f"Error message: {e}")
            return False, None
        if not all(isinstance(key[0], int) for key in rev_to_irrev_dict) or not all(isinstance(key[1], int) for key in rev_to_irrev_dict if len(key) > 1):
            warning(f"The irreversible list {rev_to_irrev_dict_name_by_options} "
                    f"does not contain integer")
            is_correct = False
        max_index = max([x[1] for x in rev_to_irrev_dict if len(x) >1])
        if len(rev_to_irrev_dict) != len(cobra_model.reactions) and max_index != len(cobra_model.reactions):
            warning(f"The irreversible list {rev_to_irrev_dict_name_by_options} "
                    f"does not contain the same amount of values as the amount of reactions in the model.")
            is_correct = False
        if not is_correct:
            info(f"Since one or more errors were found in the irreversible {rev_to_irrev_dict_name_by_options}, "
                 f"the file will be renamed to {rev_to_irrev_dict_name_by_options}_incorrect.json and "
                 f"a new irreversible dictionary will be created.")
            os.rename(rev_to_irrev_dict_name_by_options, rev_to_irrev_dict_name_by_options.split(".json")[0] + "_incorrect.json")
            return False, None
        return True, rev_to_irrev_dict

    @staticmethod
    def check_if_combination_taskstructure_is_correct(
            tasklist_file: str,
            already_loaded_tasklist: pd.DataFrame = None,
    ) -> (bool, (pd.DataFrame | None)):
        """
        Check if the created task_structure file is correct.
            - Is json readable
            - Is a dictionary with keys as sample IDs and the values as dicts containing:
            {"1.0": {"task_id": 1.0,
            "task_description": str,
            "task_should_fail_in_humans": int,
             "input_metabolites": ["O2[c]", "Pi[c]", "ADP[c]", "H+[c]"],
             "input_metabolite_lower_bounds": [6.0, 30.0, 30.0, 30.0],
              "input_metabolite_upper_bounds": [6.0, 32.0, 32.0, 32.0],
               "output_metabolites": ["H2O[c]", "CO2[c]"],
                "output_metabolite_lower_bounds": [36.0, 6.0],
                 "output_metabolite_upper_bounds": [38.0, 6.0]}

                - Task ID
                - Task description
                - Input metabolites
                - Input lower bounds
                - Input upper bounds
                - Output metabolites
                - Output lower bounds
                - Output upper bounds
        """
        is_correct = True
        task_structure = already_loaded_tasklist
        if already_loaded_tasklist is None:
            if tasklist_file.endswith(".json"):
                try:
                    task_structure = load(open(tasklist_file))
                except Exception as e:
                    warning(f"An error occurred while reading the task_structure file {tasklist_file}.\n"
                            f"Error message: {e}")
                    return False, None
            else:
                warning(f"Could not find task_structure file: {tasklist_file}.")
                return False, None
        if not isinstance(task_structure, dict):
            warning(f"The task_structure file {tasklist_file} is not a dictionary.")
            is_correct = False
        for sample_id, task in task_structure.items():
            if not isinstance(sample_id, int) and not isinstance(sample_id, str):
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a sample ID that is not a string.")
                is_correct = False
                break
            if not isinstance(task, dict):
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task that is not a dictionary.")
                is_correct = False
                break
            if "task_id" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without a task ID.")
                is_correct = False
                break
            if "task_description" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without a task description.")
                is_correct = False
                break
            if "input_metabolites" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without input metabolites.")
                is_correct = False
                break
            if "input_metabolite_lower_bounds" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without input metabolite lower bounds.")
                is_correct = False
                break
            if "input_metabolite_upper_bounds" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without input metabolite upper bounds.")
                is_correct = False
                break
            if "output_metabolites" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without output metabolites.")
                is_correct = False
                break
            if "output_metabolite_lower_bounds" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without output metabolite lower bounds.")
                is_correct = False
                break
            if "output_metabolite_upper_bounds" not in task:
                warning(f"The task_structure file {tasklist_file}, sample {sample_id} contains a task without output metabolite upper bounds.")
                is_correct = False
                break
        if not is_correct:
            info(f"Since one or more errors were found in the task_structure file {tasklist_file}, "
                 f"the file will be renamed to {tasklist_file}_incorrect.json and "
                 f"a new task_structure file will be created.")
            os.rename(tasklist_file, tasklist_file.split(".json")[0] + "_incorrect.json")
            return False, None
        return True, task_structure

    @staticmethod
    def check_if_tasklist_file_is_correct(
            tasklist_file_location: str,
            already_loaded_tasklist: pd.DataFrame = None,
    ):
        """
        Check if the tasklist file is correct.
            - Is json readable
            - "ID" as the first column (starting from 1 and positive integers only)
            - "DESCRIPTION" as the second column with task names/descriptions
            - "IN" as the fourth column with the input metabolites in MAM0xxxx notation
            - "IN LB" as the fifth column with the lower bounds of the input metabolites
            - "IN UB" as the sixth column with the upper bounds of the input metabolites
            - "OUT" as the seventh column with the output metabolites in MAM0xxxx notation
            - "OUT LB" as the eighth column with the lower bounds of the output metabolites
            - "OUT UB" as the ninth column with the upper bounds of the output metabolites
            - Any columns starting with "CLASSIFIER_" are used to classify the tasks

        :param tasklist_file: str
            The path to the tasklist file.
        :param already_loaded_tasklist: dict
            The already loaded tasklist.
        """
        is_correct = True
        tasklist = None
        if already_loaded_tasklist is None:
            for file in os.listdir(tasklist_file_location):
                if file.startswith("tasklist"):
                    if file.endswith(".json"):
                        tasklist = pd.read_json(os.path.join(tasklist_file_location, file))
                        break
                    elif file.endswith(".csv"):
                        tasklist = pd.read_csv(os.path.join(tasklist_file_location, file))
                        break
                    elif file.endswith(".xlsx"):
                        try:
                            tasklist = pd.read_excel(os.path.join(tasklist_file_location, file))
                        except Exception as e:
                            warning(f"An error occurred while reading the tasklist file {file}.\n"
                                    f"Error message: {e}")
                            return False, None
        else:
            tasklist = already_loaded_tasklist
        if tasklist is None:
            raise ValueError(f"The tasklist file in {tasklist_file_location} is not in a readable format.")

        if "ID" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'ID'.")
            is_correct = False
        if "DESCRIPTION" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'DESCRIPTION'.")
            is_correct = False
        if "IN" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'IN'.")
            is_correct = False
        if "IN LB" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'IN LB'.")
            is_correct = False
        if "IN UB" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'IN UB'.")
            is_correct = False
        if "OUT" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'OUT'.")
            is_correct = False
        if "OUT LB" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'OUT LB'.")
            is_correct = False
        if "OUT UB" not in tasklist.columns:
            warning(f"The tasklist file {tasklist_file_location} does not contain the column 'OUT UB'.")
            is_correct = False
        if not is_correct:
            raise("The tasklist file is not correct, please input a correct tasklist file.")

    @staticmethod
    def check_if_combination_files_are_correct(
            check_if_all_correct_dict: dict,
            names_dict: dict,
            global_options_dict: dict,
            inputs_dict: dict = None,
    ) -> (dict, dict):

        if inputs_dict is not None and all(value is not None for value in inputs_dict.values()):
            expression_data = inputs_dict["expression_data"]
            cobra_model = inputs_dict["model"]
            tasklist = inputs_dict["tasklist"]
            if global_options_dict["__global_options"]["model_type"] == "irr":
                rev_to_irrev_dict = inputs_dict["rev_to_irrev_dict"]
        else:
            inputs_dict = {}
            expression_data = None
            cobra_model = None
            tasklist = None
            if global_options_dict["__global_options"]["model_type"] == "irr":
                rev_to_irrev_dict = None

        expression_name_by_options = names_dict["expression_data"]
        model_name_by_options = names_dict["model"]
        tasklist_name_by_options = names_dict["tasklist"]

        (is_correct, expression_data) = FileChecker.check_if_expression_data_is_correct(
            expression_name_by_options,
            global_options_dict["__global_options"]["alternative_reaction_notation"],
            expression_data
        )
        if not is_correct:
            check_if_all_correct_dict["expression_data"] = False
        (is_correct, cobra_model) = FileChecker.check_if_combination_model_is_correct(model_name_by_options, cobra_model)
        if not is_correct:
            check_if_all_correct_dict["model"] = False
        (is_correct) = FileChecker.check_if_combination_expression_and_model_are_congruent(
            cobra_model,
            expression_data,
        )
        if not is_correct:
            check_if_all_correct_dict["expression_data"] = False
            check_if_all_correct_dict["model"] = False
        if global_options_dict["__global_options"]["model_type"] == "irr":
            rev_to_irrev_dict_name_by_options = names_dict["rev_to_irrev_dict"]
            (is_correct, rev_to_irrev_dict) = FileChecker.check_if_irrev_dict_is_correct(
                rev_to_irrev_dict_name_by_options,
                cobra_model,
            )
            inputs_dict["rev_to_irrev_dict"] = rev_to_irrev_dict
        if not is_correct:
            check_if_all_correct_dict["rev"] = False

        (is_correct, tasklist) = FileChecker.check_if_combination_taskstructure_is_correct(tasklist_name_by_options, tasklist)
        if not is_correct:
            check_if_all_correct_dict["tasklist"] = False

        inputs_dict["expression_data"] = expression_data
        inputs_dict["model"] = cobra_model
        inputs_dict["tasklist"] = tasklist

        return check_if_all_correct_dict, inputs_dict


class CombinationFileCreation:

    @staticmethod
    def create_combination_directory(
            main_data_folder: str,
            global_options_dict: dict,
    ) -> str:
        """
        Create the "combinations" directory and populate it with the necessary files for the preprocessing.
            - Creates a directory for the combination of the model, expression data, and tasklist.

        :param main_data_folder: str
            The main directory where the data is stored.
        :param global_options_dict: dict
            The dictionary with the global options for the preprocessing.
        """
        selected_model = global_options_dict["__global_options"]["model_version"]
        selected_expression_dataset = global_options_dict["__global_options"]["expression_dataset"]
        selected_tasklist = global_options_dict["__global_options"]["tasklist"]
        combinations_folder = os.path.join(main_data_folder, "combinations")
        os.mkdir(combinations_folder) if not os.path.exists(combinations_folder) else None
        combination_name = f"{selected_model}_{selected_expression_dataset}_{selected_tasklist}"
        combination_folder = os.path.join(combinations_folder, combination_name)
        os.mkdir(combination_folder) if not os.path.exists(combination_folder) else None
        return combination_folder

    @staticmethod
    def check_if_model_is_irreversible(
            cobra_model: Model,
    ) -> bool:
        """
        Check if the model is irreversible.

        :param cobra_model: Model
            The cobrapy model.
        """
        return all(reaction.lower_bound >= 0 for reaction in cobra_model.reactions if reaction.upper_bound > 0)

    @staticmethod
    def convert_reversible_model_to_irreversible(
            cobra_model: Model,
    ) -> (Model, dict):
        """
        Convert the reversible model to an irreversible model.

        :param cobra_model: Model
            The cobrapy model.
        :return: Model
            The irreversible cobrapy model.
        :return: dict
            The dictionary linking each reaction's index with its corresponding irreversible reaction index.
        """
        irreversible_model = cobra_model.copy()
        reactions_to_modify = [
            rxn for rxn in irreversible_model.reactions if rxn.lower_bound < 0 < rxn.upper_bound
        ]
        to_create_rev2irrev = [[idx + 1] for idx in range(len(cobra_model.reactions))]
        rxn_subs = [0] * len(reactions_to_modify)
        idx_2_counter = 0
        length_of_model = len(irreversible_model.reactions)
        start_time = time.time()
        for rxn in reactions_to_modify:
            idx = irreversible_model.reactions.index(rxn)
            new_reaction = Reaction(rxn.id + "_r")
            new_reaction.subsystem = rxn.subsystem
            new_reaction.gene_reaction_rule = rxn.gene_reaction_rule
            new_reaction.name = rxn.name
            new_reaction.bounds = (0, -rxn.lower_bound)
            new_reaction.add_metabolites(
                dict((met, -1 * coeff) for met, coeff in rxn.metabolites.items())
            )
            rxn.bounds = (0, rxn.upper_bound)
            rxn.id = rxn.id + "_f"
            rxn_subs[idx_2_counter] = new_reaction
            to_create_rev2irrev[idx].append(idx_2_counter + 1 + length_of_model)
            idx_2_counter += 1

        # for idx, rxn in enumerate(cobra_model.reactions):
        #     if rxn.lower_bound < 0 < rxn.upper_bound:
        #         new_reaction = Reaction(rxn.id + "_r")
        #         new_reaction.subsystem = rxn.subsystem
        #         new_reaction.gene_reaction_rule = rxn.gene_reaction_rule
        #         new_reaction.name = rxn.name
        #         new_reaction.bounds = (0, -rxn.lower_bound)
        #         new_reaction.add_metabolites(
        #             dict((met, -1 * coeff) for met, coeff in rxn.metabolites.items())
        #         )
        #         rxn.bounds = (0, rxn.upper_bound)
        #         rxn_subs[idx_2_counter] = new_reaction
        #         to_create_rev2irrev[idx].append(idx_2_counter + 1 + length_of_model)
        #         idx_2_counter += 1
        print(f"Time to create reactions: {time.time() - start_time}")
        start_time = time.time()
        irreversible_model.add_reactions(rxn_subs)
        print(f"Time to add reactions: {time.time() - start_time}")

        return irreversible_model, to_create_rev2irrev

    @staticmethod
    def create_rev2irrev_dict(
            cobra_model: Model,
    ) -> dict:
        """
        Create the irreversible dictionary.

        :param cobra_model: Model
            The cobrapy model.
        :return: dict
            The dictionary linking each reaction's index with its corresponding irreversible reaction index.
        """
        if CombinationFileCreation.check_if_model_is_irreversible(cobra_model):
            # find all reactions with _f and find the index with _r or _b and add to dictionary
            reaction_pair_indices = {idx: idx2 for idx in range(len(cobra_model.reactions)) for idx2, rxn in enumerate(cobra_model.reactions) if
                                     rxn.id == cobra_model.reactions[idx].id + "_r" or rxn.id == cobra_model.reactions[idx].id + "_b"}
            return reaction_pair_indices
        else:
            cobra_model, reaction_pair_indices = CombinationFileCreation.convert_reversible_model_to_irreversible(cobra_model)
            return reaction_pair_indices

    @staticmethod
    def create_closed_model(
            cobra_model: Model,
            irreversible: bool = False,
    ) -> (Model, dict):
        """
        Create a closed model from the cobrapy model
        and if necessary convert it to an irreversible model.

        :param cobra_model: Model
            The cobrapy model.
        :param global_options_dict: dict
            The dictionary with the global options for the preprocessing.

        :return: (Model, dict)
            The closed cobrapy model and the dictionary linking each reaction's
             index with its corresponding irreversible reaction index.
        """
        for reaction in cobra_model.exchanges:
            reaction.lower_bound = 0
            reaction.upper_bound = 0
        for reaction in cobra_model.boundary:
            reaction.lower_bound = 0
            reaction.upper_bound = 0
        for reaction in cobra_model.reactions:
            if reaction.lower_bound is None:
                reaction.lower_bound = -1000
            if reaction.upper_bound is None:
                reaction.upper_bound = 1000
            if reaction.lower_bound != 0 and reaction.upper_bound != 0 and (
                    (reaction.upper_bound > 0 and reaction.lower_bound > 0) or
                    (reaction.upper_bound < 0 and reaction.lower_bound < 0)
            ):
                if reaction.lower_bound < 0 and reaction.upper_bound < 0:
                    reaction.lower_bound = -1000
                    reaction.upper_bound = 0
                elif reaction.lower_bound > 0 and reaction.upper_bound > 0:
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
        if irreversible:
            cobra_model, rev_to_irrev_dict = CombinationFileCreation.convert_reversible_model_to_irreversible(cobra_model)
        else:
            rev_to_irrev_dict = None
        return cobra_model, rev_to_irrev_dict

    @staticmethod
    def set_threshold_and_apply_to_data(
    ):
        pass

    @staticmethod
    def tokenize_gpr(expression):
        # Tokenize by parentheses and logical operators
        tokens = re.findall(r'\(|\)|\band\b|\bor\b|[^\s()]+', expression)
        return tokens

    @staticmethod
    def parse_gpr(tokens):
        def parse_expression(tokens):
            stack = []
            operator = None
            while tokens:
                token = tokens.pop(0)

                if token == "(":
                    # Recursively parse the sub-expression inside parentheses
                    sub_expr = parse_expression(tokens)
                    stack.append(sub_expr)

                elif token == ")":
                    # End of the current expression; return what we have parsed so far
                    break

                elif token == "and" or token == "or":
                    # Set the current operator ("and"/"or")
                    operator = token

                else:
                    # It's a gene identifier
                    stack.append(token)

                if operator and len(stack) > 1:
                    # Apply the operator to the two most recent items in the stack
                    right = stack.pop()
                    left = stack.pop()
                    stack.append((operator, [left, right]))
                    operator = None

            return stack[0] if stack else None

        return parse_expression(tokens)

    @staticmethod
    def trim_and_min(values, trim_percentage):
        sorted_values = sorted(values)
        trim_count = int(len(sorted_values) * trim_percentage)
        if len(sorted_values) < 2:
            return min(sorted_values)
        elif trim_count == len(sorted_values):
            return min(sorted_values)
        elif trim_count == 0:
            return min(sorted_values)
        elif trim_count >= len(sorted_values) / 2:
            trimmed_values = sorted_values[trim_count:len(sorted_values)]
            return min(trimmed_values)
        trimmed_values = sorted_values[trim_count:len(sorted_values) - trim_count]
        return min(trimmed_values) if trimmed_values else None

    @staticmethod
    def trim_and_mean(values, trim_percentage):
        sorted_values = sorted(values)
        trim_count = int(len(sorted_values) * trim_percentage)
        if len(sorted_values) < 2:
            return np.mean(sorted_values)
        elif trim_count == len(sorted_values):
            return np.mean(sorted_values)
        elif trim_count == 0:
            return np.mean(sorted_values)
        elif trim_count >= len(sorted_values) / 2:
            trimmed_values = sorted_values[trim_count:len(sorted_values)]
            return np.mean(trimmed_values)
        trimmed_values = sorted_values[trim_count:len(sorted_values) - trim_count]
        return np.mean(trimmed_values) if trimmed_values else None

    @staticmethod
    def evaluate_contribution(
            rule: str,
            rule_value_1: (float | int | None),
            rule_value_2: (float | int | None),
            values: dict,
            gene_contribution_dict_input: dict,
            in_promiscuity: bool = False,
    ):
        '''
        Evaluate the contribution of the gene to the reaction,
        changes the gene_contribution_dict_input without returning anything

        :param operator: str
            The operator (and/or)
        :param rule: str
            The rule (min/trim)
        :param values: dict
            The dict containing operator and values
        :param gene_contribution_dict_input: dict
            The dictionary with the gene IDs and their contribution to all reactions in the sample.
        '''
        if in_promiscuity:
            # in the case of trimmean, it should resolve to mean with only those that contribute #(all vals that are not 0 contribution)
            # in the case of trimmin it will be a single value with min but that is also the only value (1 value is not 0 contribution)
            # int he case of trimmax it with be a single value with max but that is also the only value (1 value is not 0 contribution)
            # in the case of min choose the single value already chosen (1 is not 0 contribution)
            # in the case of mean choose the mean of all values that are not 0 contribution (all are not 0 contribution)
            # in the case of max choose the max of all values that are not 0 contribution (1 is not 0 contribution)
            # in the case of sum choose the sum of all values that are not 0 contribution (all are not 0 contribution)
            # in the case of median choose the median of all values that are not 0 contribution (all are not 0 contribution)
            if rule == "min1" or rule == "min2":
                min_value = min(values.values())
                for key, value in values.items():
                    if value is None:
                        continue
                    if value == min_value:
                        gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
                    else:
                        gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1) if key in gene_contribution_dict_input else 0
            elif rule == "mean" or rule == "median" or rule == "sum":
                for key, value in values.items():
                    if value is None:
                        continue
                    gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
            elif rule == "max":
                max_value = max(values.values())
                for key, value in values.items():
                    if value is None:
                        continue
                    if value == max_value:
                        gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
                    else:
                        gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1) if key in gene_contribution_dict_input else 0
            elif rule == "trimmean" or rule == "trimmin1" or rule == "trimmin2" or rule == "trimmax":
                values_without_none = {key: value for key, value in values.items() if value is not None}
                sorted_values = sorted(values_without_none.values())
                trim_count = int(len(sorted_values) * rule_value_1)
                if trim_count >= len(sorted_values) / 2:
                    trimmed_values = sorted_values[trim_count:len(sorted_values)]
                elif len(sorted_values) < 2 or trim_count == len(sorted_values) or trim_count == 0:
                    trimmed_values = sorted_values
                else:
                    trimmed_values = sorted_values[trim_count:len(sorted_values) - trim_count]
                if rule == "trimmean":
                    for key, value in values.items():
                        if value is None:
                            continue
                        if value in trimmed_values:
                            gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
                        else:
                            gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1) if key in gene_contribution_dict_input else 0
                elif rule == "trimmin1" or rule == "trimmin2":
                    min_value = min(trimmed_values)
                    for key, value in values.items():
                        if value is None:
                            continue
                        if value == min_value:
                            gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
                        else:
                            gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1) if key in gene_contribution_dict_input else 0
                elif rule == "trimmax":
                    max_value = max(trimmed_values)
                    for key, value in values.items():
                        if value is None:
                            continue
                        if value == max_value:
                            gene_contribution_dict_input[key] = (gene_contribution_dict_input[key]+1) if key in gene_contribution_dict_input else 1
                        else:
                            gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1) if key in gene_contribution_dict_input else 0

        elif not in_promiscuity:
            if rule == "min1":
                min_value = min(values.values())
                for key, value in values.items():
                    if value == min_value:
                        if key not in gene_contribution_dict_input:
                            gene_contribution_dict_input[key] = 1
                        else:
                            gene_contribution_dict_input[key] += 1
                    else:
                        if key not in gene_contribution_dict_input:
                            gene_contribution_dict_input[key] = 0
                        else:
                            gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1)

            elif rule == "min2":
                for key, value in values.items():
                    if key not in gene_contribution_dict_input:
                        gene_contribution_dict_input[key] = 1
                    else:
                        gene_contribution_dict_input[key] += 1

            elif rule == "mean" or rule == "median":
                val = 1 / len(values)
                for key, value in values.items():
                    if key not in gene_contribution_dict_input:
                        gene_contribution_dict_input[key] = val
                    else:
                        gene_contribution_dict_input[key] = val
            elif rule == "max":
                max_value = max(values.values())
                for key, value in values.items():
                    if value == max_value:
                        if key not in gene_contribution_dict_input:
                            gene_contribution_dict_input[key] = 1
                        else:
                            gene_contribution_dict_input[key] += 1
                    else:
                        if key not in gene_contribution_dict_input:
                            gene_contribution_dict_input[key] = 0
                        else:
                            gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1)

            elif rule == "trimmean" or rule == "trimmin1" or rule == "trimmin2" or rule == "trimmax":
                sorted_values = sorted(values.values())
                trim_count = int(len(sorted_values) * rule_value_1)
                if trim_count >= len(sorted_values) / 2:
                    trimmed_values = sorted_values[trim_count:len(sorted_values)]
                elif len(sorted_values) < 2 or trim_count == len(sorted_values) or trim_count == 0:
                    trimmed_values = sorted_values
                else:
                    trimmed_values = sorted_values[trim_count:len(sorted_values) - trim_count]
                if rule == "trimmean":
                    val = 1 / len(trimmed_values)
                    for key, value in values.items():
                        if value in trimmed_values:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = val
                            else:
                                gene_contribution_dict_input[key] += val
                        else:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = 0
                            else:
                                gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1)
                elif rule == "trimmin1":
                    min_value = min(trimmed_values)
                    for key, value in values.items():
                        if value == min_value:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = 1
                            else:
                                gene_contribution_dict_input[key] += 1
                        else:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = 0
                            else:
                                gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1)
                elif rule == "trimmin2":
                    for key, value in values.items():
                        if key not in gene_contribution_dict_input:
                            gene_contribution_dict_input[key] = 1
                        else:
                            gene_contribution_dict_input[key] += 1
                elif rule == "trimmax":
                    max_value = max(trimmed_values)
                    for key, value in values.items():
                        if value == max_value:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = 1
                            else:
                                gene_contribution_dict_input[key] += 1
                        else:
                            if key not in gene_contribution_dict_input:
                                gene_contribution_dict_input[key] = 0
                            else:
                                gene_contribution_dict_input[key] = max(0, gene_contribution_dict_input[key]-1)
            elif rule == "sum":
                total_sum = sum(values.values())
                for key, value in values.items():
                    contribution = value / total_sum
                    if key not in gene_contribution_dict_input:
                        gene_contribution_dict_input[key] = contribution
                    else:
                        gene_contribution_dict_input[key] = contribution

    @staticmethod
    def create_string_from_parsed_expression(parsed_expr: tuple | str) -> str:
        if isinstance(parsed_expr, tuple):
            operator, operands = parsed_expr
            return f"({CombinationFileCreation.create_string_from_parsed_expression(operands[0])} {operator} {CombinationFileCreation.create_string_from_parsed_expression(operands[1])})"
        else:
            return parsed_expr

    @staticmethod
    def evaluate_expression(
            parsed_expr: tuple,
            gene_values: dict,
            gene_contribution_dict: dict,
            gene_contribution_dict_in_reaction_for_promiscuity: dict,
            or_rule: str,
            and_rule: str,
            rule_value: (float | int | None),
            rule_value_2: (float | int | None),
            for_promiscuity: bool = False,
    ) -> (float, list):
        '''
        Warning this function has side (intentional) effects on the gene_contribution_dict
        Evaluate the parsed expression and apply the given rules and gene values

        The GPR method argument.
               min1Sum (taking the sum of OR values) min1 indicates to only take the min value as contributor
               min2Sum (taking the sum of OR values) min2 indicates that all values of the AND min rule are contributors
               min1Max (taking the maximum of OR values),
               min1Mean (taking the mean of OR values),
               min1Median (taking the median of OR values),
               trimmeanSum_V (taking the mean of AND values after trimming V%),
                            and then taking the sum of OR values,
               trimmeanMax_V (taking the mean of AND values after trimming V%),
                            and then taking the maximum of OR values
               trimminSum_V (taking the min of AND values after trimming V%),
                            and then taking the sum of OR values
               trimminMax_V (taking the min of AND values after trimming V%),
                            and then taking the maximum of OR values
               min2Sum (taking the sum of OR values),
        :param parsed_expr: tuple
            The parsed expression.
        :param gene_values: dict
            The dictionary with the gene IDs and their value in the sample.
        :param gene_contribution_dict: dict
            The dictionary with the gene IDs and their contribution to all reactions in the sample.
        :param or_rule: str
            The OR rule.
        :param and_rule: str
            The AND rule.
        :param rule_value: float
            The rule value.
        :param rule_value_2: float
            The second rule value.
        :return value: float
            The value of the expression after applying the GPR rules
        '''
        # TODO add test for all options (iterate through and check for a single sample if these work)
        if isinstance(parsed_expr, tuple):
            operator, operands = parsed_expr
            values_dict = {
                CombinationFileCreation.create_string_from_parsed_expression(op): CombinationFileCreation.evaluate_expression(
                    op, gene_values, gene_contribution_dict,
                    gene_contribution_dict_in_reaction_for_promiscuity,
                    or_rule, and_rule, rule_value,
                    rule_value_2,
                    for_promiscuity
                ) for op in operands
            }
            values_dict = {key: value for key, value in values_dict.items() if value is not None}
            if operator == "or":
                rule = or_rule
            elif operator == "and":
                rule = and_rule
            CombinationFileCreation.evaluate_contribution(
                rule, rule_value,
                rule_value_2, values_dict,
                gene_contribution_dict,
            )
            if for_promiscuity:
                CombinationFileCreation.evaluate_contribution(
                    rule, rule_value,
                    rule_value_2, values_dict,
                    gene_contribution_dict_in_reaction_for_promiscuity,
                    for_promiscuity,
                )

            if operator == "or":
                if or_rule == "sum":
                    if not values_dict.values():
                        return None
                    values = sum(values_dict.values())
                    return values
                elif or_rule == "max":
                    values = max(values_dict.values())
                    return values
                elif or_rule == "mean":
                    values = sum(values_dict.values()) / len(values_dict)
                    return values
                elif or_rule == "median":
                    values = np.median(list(values_dict.values()))
                    return values
            elif operator == "and":

                if and_rule == "min1" or and_rule == "min2":
                    if not values_dict.values():
                        return None
                    values = min(values_dict.values())
                    return values
                elif and_rule == "trimmean":
                    values = CombinationFileCreation.trim_and_mean(list(values_dict.values()), rule_value)
                    return values
                elif and_rule == "trimmin":
                    values = CombinationFileCreation.trim_and_min(list(values_dict.values()), rule_value)
                    return values
        if gene_values.get(parsed_expr, None) is not None:
            if parsed_expr not in gene_contribution_dict:
                gene_contribution_dict[parsed_expr] = 1
            else:
                gene_contribution_dict[parsed_expr] = 1
            if for_promiscuity:
                if parsed_expr not in gene_contribution_dict_in_reaction_for_promiscuity:
                    gene_contribution_dict_in_reaction_for_promiscuity[parsed_expr] = 1
                else:
                    gene_contribution_dict_in_reaction_for_promiscuity[parsed_expr] += 1
        return gene_values.get(parsed_expr, None)

    @staticmethod
    def evaluate_expression_in_promiscuity(
            parsed_expr: tuple,
            adjusted_gene_values: dict,
            gene_contribution_dict_in_reaction_input: dict,
            GPR_OR_argument: str,
            GPR_AND_argument: str,
            GPR_value: (float | int | None),
            GPR_value_2: (float | int | None),
    ):
        """
        Evaluate the parsed expression and apply the given rules and gene values

        :param parsed_expr: tuple
            The parsed expression.
        :param adjusted_gene_values: dict
            The dictionary with the gene IDs and their adjusted value in the sample.
        :param gene_contribution_dict_in_reaction: dict
            The dictionary with the gene IDs and their contribution to the reaction in the sample.
        :param GPR_OR_argument: str
            The OR rule.
        :param GPR_AND_argument: str
            The AND rule.
        :param GPR_value: float
            The rule value.
        :param GPR_value_2: float
            The second rule value.
        :return value: float
            The value of the expression after applying the GPR rules
        """
        if isinstance(parsed_expr, tuple):
            operator, operands = parsed_expr
            for op in operands:
                if CombinationFileCreation.create_string_from_parsed_expression(op) in gene_contribution_dict_in_reaction_input.keys():
                    if gene_contribution_dict_in_reaction_input[CombinationFileCreation.create_string_from_parsed_expression(op)] > 0:
                        pass
            values_dict = {
                CombinationFileCreation.create_string_from_parsed_expression(op): CombinationFileCreation.evaluate_expression_in_promiscuity(
                    op, adjusted_gene_values, gene_contribution_dict_in_reaction_input,
                    GPR_OR_argument, GPR_AND_argument, GPR_value, GPR_value_2
                ) for op in operands
                if CombinationFileCreation.create_string_from_parsed_expression(op)
                   in gene_contribution_dict_in_reaction_input.keys()
                   and gene_contribution_dict_in_reaction_input[CombinationFileCreation.create_string_from_parsed_expression(op)] > 0
            }
            values_dict = {key: value for key, value in values_dict.items() if value is not None}
            if len(values_dict) == 0:
                return None
            if operator == "or":
                if GPR_OR_argument == "sum":
                    values = sum(values_dict.values())
                    return values
                elif GPR_OR_argument == "max":
                    values = max(values_dict.values())
                    return values
                elif GPR_OR_argument == "mean":
                    values = sum(values_dict.values()) / len(values_dict)
                    return values
                elif GPR_OR_argument == "median":
                    values = np.median(list(values_dict.values()))
                    return values
            elif operator == "and":
                if GPR_AND_argument == "min1" or GPR_AND_argument == "min2":
                    values = min(values_dict.values())
                    return values
                elif GPR_AND_argument == "trimmean":
                    values = sum(values_dict.values()) / len(values_dict)
                    return values
                elif GPR_AND_argument == "trimmin":
                    values = min(values_dict.values())
                    return values
        return adjusted_gene_values.get(parsed_expr, None)

    @staticmethod
    def parse_gene_contribution_dict(
            gene_contribution_dict: dict,
    ) -> dict:
        """
        Parse the gene contribution dictionary by multiplying genes' contributions

        :param gene_contribution_dict: dict
            The dictionary with operands containing and their contribution to all reactions in the sample.
        """
        # split "(", "or", "and", ")"
        genes_in_reaction = []
        for key, value in gene_contribution_dict.items():
            potential_genes = key.replace("(", "").replace(")", "").replace("or", "").replace("and", "").split()
            if len(potential_genes) == 1:
                genes_in_reaction.append(potential_genes[0]) if potential_genes and potential_genes not in genes_in_reaction else None

        parsed_gene_contribution_dict = {}
        for gene in genes_in_reaction:
            parsed_gene_contribution_dict[gene] = 1
            for key, value in gene_contribution_dict.items():
                if gene in key:
                    parsed_gene_contribution_dict[gene] *= value

        return parsed_gene_contribution_dict

    @staticmethod
    def allocate_gene_values_to_reactions(
            cobra_model: Model,
            gene_score: pd.DataFrame,
            expression_data: pd.DataFrame,
            gene_ids: list,
            sample_number: int,
            promiscuity_method: str,
            GPR_OR_argument: str,
            GPR_AND_argument: str,
            GPR_value: (float | int | None),
            GPR_value_2: (float | int | None),
    ) -> (np.ndarray, np.ndarray,
          list[list[dict]], list[list[dict]]):

        gene_contribution_matrix = np.zeros((len(gene_ids), sample_number))
        reaction_expression_matrix = np.zeros((len(cobra_model.reactions), sample_number))
        per_reaction_gene_contribution_tuple_matrix = [[{} for _ in range(sample_number)] for _ in range(len(cobra_model.reactions))]
        per_reaction_adjusted_gene_contribution_tuple_matrix = [[{} for _ in range(sample_number)] for _ in range(len(cobra_model.reactions))]
        promiscuity_correction = True if promiscuity_method == "Y" else False

        for idx in range(expression_data.shape[1]):
            genes_in_reaction_dict_for_promiscuity = {}
            per_reaction_gene_contribution_in_sample = {reaction.id: {} for reaction in cobra_model.reactions}
            gene_contribution_dict_in_sample = {gene: 0 for gene in gene_ids}
            adjusted_expression_values_in_sample = {reaction: 0 for reaction in cobra_model.reactions}
            expression_in_sample = gene_score.iloc[:, idx]
            gene_values = dict(zip(gene_ids, expression_in_sample))
            reaction_token_rule_dict = {
                reaction.id: CombinationFileCreation.tokenize_gpr(reaction.gene_reaction_rule)
                for reaction in cobra_model.reactions
            }

            for idx2 in range(len(cobra_model.reactions)):
                rxn = cobra_model.reactions[idx2]
                tokens = reaction_token_rule_dict[rxn.id]
                if len(tokens) == 0:
                    result = -1
                else:
                    gene_contribution_dict_in_reaction = {}
                    gene_contribution_dict_in_reaction_for_promiscuity = {}
                    parsed_expr = CombinationFileCreation.parse_gpr(tokens)
                    result = CombinationFileCreation.evaluate_expression(
                        parsed_expr, gene_values,
                        gene_contribution_dict_in_reaction,
                        gene_contribution_dict_in_reaction_for_promiscuity,
                        GPR_OR_argument, GPR_AND_argument,
                        GPR_value, GPR_value_2,
                        promiscuity_correction
                    )
                    string_parsed_expr = str(parsed_expr[0]) + str(parsed_expr[1])
                    gene_contribution_dict_in_reaction = {key: min(1, value) for key, value in gene_contribution_dict_in_reaction.items()}

                    genes_in_reaction_dict_for_promiscuity[string_parsed_expr] = gene_contribution_dict_in_reaction_for_promiscuity
                    parsed_gene_contribution_dict_in_reaction = CombinationFileCreation.parse_gene_contribution_dict(
                        gene_contribution_dict_in_reaction
                    )
                    current_reaction_dict = {}

                    for key, value in parsed_gene_contribution_dict_in_reaction.items():
                        gene_contribution_dict_in_sample[key] += value
                        current_reaction_dict[key] = (key, value, gene_values[key])
                    per_reaction_gene_contribution_in_sample[rxn.id] = current_reaction_dict
                if promiscuity_method == "N":
                    adjusted_expression_values_in_sample[cobra_model.reactions[idx2]] = result
            gene_contribution_matrix[:, idx] = [gene_contribution_dict_in_sample[gene] for gene in gene_ids]
            for i, reaction in enumerate(cobra_model.reactions):
                per_reaction_gene_contribution_tuple_matrix[i][idx] = per_reaction_gene_contribution_in_sample[reaction.id]

            if promiscuity_method == "N":
                reaction_expression_matrix[:, idx] = [adjusted_expression_values_in_sample[reaction] for reaction in cobra_model.reactions]
                for i, reaction in enumerate(cobra_model.reactions):
                    per_reaction_adjusted_gene_contribution_tuple_matrix[i][idx] = per_reaction_gene_contribution_in_sample[reaction.id]
            elif promiscuity_method == "Y":
                reaction_token_rule_dict = {
                    reaction.id: CombinationFileCreation.tokenize_gpr(reaction.gene_reaction_rule)
                    for reaction in cobra_model.reactions
                }
                gene_contribution_dict_in_sample = {key: max(1, value) for key, value in gene_contribution_dict_in_sample.items()}
                adjusted_gene_values_original = {
                    gene: gene_values[gene] / gene_contribution_dict_in_sample[gene]
                    for gene in gene_values.keys()
                }
                for idx2 in range(len(cobra_model.reactions)):
                    tokens = reaction_token_rule_dict[cobra_model.reactions[idx2].id]
                    if len(tokens) == 0:
                        result = -1
                    else:
                        parsed_expr = CombinationFileCreation.parse_gpr(tokens)
                        string_parsed_expr = str(parsed_expr[0]) + str(parsed_expr[1])
                        promiscuity_gene_contribution_dict_in_reaction = genes_in_reaction_dict_for_promiscuity[string_parsed_expr]
                        # in the case of trimmean, it should resolve to mean with only those that contribute #(all vals that are not 0 contribution)
                        # in the case of trimmin it will be a single value with min but that is also the only value (1 value is not 0 contribution)
                        # int he case of trimmax it with be a single value with max but that is also the only value (1 value is not 0 contribution)
                        # in the case of min choose the single value already chosen (1 is not 0 contribution)
                        # in the case of mean choose the mean of all values that are not 0 contribution (all are not 0 contribution)
                        # in the case of max choose the max of all values that are not 0 contribution (1 is not 0 contribution)
                        # in the case of sum choose the sum of all values that are not 0 contribution (all are not 0 contribution)
                        # in the case of median choose the median of all values that are not 0 contribution (all are not 0 contribution)
                        rxn1 = cobra_model.reactions[idx2]
                        if idx2 == 12923:
                            pass
                        result = CombinationFileCreation.evaluate_expression_in_promiscuity(
                            parsed_expr, adjusted_gene_values_original,
                            promiscuity_gene_contribution_dict_in_reaction,
                            GPR_OR_argument, GPR_AND_argument,
                            GPR_value, GPR_value_2,
                        )
                    adjusted_expression_values_in_sample[cobra_model.reactions[idx2]] = result
                    dict_ = per_reaction_gene_contribution_tuple_matrix[idx2][idx]
                    new_reaction_dict = {}
                    for key, value in dict_.items():
                        new_reaction_dict[key] = (value[0], value[1], adjusted_gene_values_original[key])
                        # per_reaction_adjusted_gene_contribution_tuple_matrix[idx2][idx][key] = (key, value, adjusted_gene_values_original[key])
                    per_reaction_adjusted_gene_contribution_tuple_matrix[idx2][idx] = new_reaction_dict

                reaction_expression_matrix[:, idx] = [adjusted_expression_values_in_sample[reaction] for reaction in cobra_model.reactions]

            print(f"sample {idx} done")
        return gene_contribution_matrix, reaction_expression_matrix, per_reaction_gene_contribution_tuple_matrix, per_reaction_adjusted_gene_contribution_tuple_matrix

    @staticmethod
    def create_GPR_expression_data(
            expression_data_location: str,
            cobra_model: Model,
            global_options_dict: dict,
            selected_model: str,
    ) -> (np.ndarray, np.ndarray, list,
          list[list[dict]], list[list[dict]]):
        """
        Create the GPR expression data for the model.

        :param expression_data_location: str
            The path to the expression data.
        :param cobra_model: Model
            The cobrapy model.
        :param global_options_dict: dict
            The dictionary with the global options for the preprocessing.
        """

        def zscore_columnwise(data):
            return data.apply(zscore, axis=0)

        def zscore_dataframe(data):
            overall_mean = data.values.mean()  # Mean of all values
            overall_std = data.values.std()  # Std dev of all values
            return (data - overall_mean) / overall_std

        model_is_irreversible = CombinationFileCreation.check_if_model_is_irreversible(cobra_model)
        if model_is_irreversible:
            model_location = os.path.join(main_data_folder, "models", selected_model)
            FileChecker.raise_error_if_incorrect_folders(model_location)
            cobra_model = FileChecker.read_model_file_in_folder(None, model_location, None)
            cobra_model, rev_to_irrev_dict = CombinationFileCreation.create_closed_model(cobra_model, True)

        expression_files = [file for file in os.listdir(expression_data_location) if file.lower().startswith("data")]
        for file in expression_files:
            file = os.path.join(expression_data_location, file)
            if file.endswith(".json"):
                expression_data = pd.read_json(file)
                break
            elif file.endswith(".xlsx"):
                expression_data = pd.read_excel(file)
                break
            elif file.endswith(".csv"):
                expression_data = pd.read_csv(file)
                break
            else:
                raise ValueError(f"The expression data file {file} is not in a readable format.")

        transformation = global_options_dict["__global_options"]["transformation"]
        thresholding = global_options_dict["__global_options"]["thresholding"]
        thresholding_method = thresholding.split("_")[0]
        thresholding_value = thresholding.split("_")[1] if len(thresholding.split("_")) > 1 else None
        thresholding_value_2 = thresholding.split("_")[2] if len(thresholding.split("_")) > 2 else None
        if thresholding_method in ["GPTL", "GPTG, GTV"]:
            local_global_no = "global"
        elif thresholding_method in ["LTM", "LTMMP", "LTMMV"]:
            local_global_no = "local"
        elif thresholding_method in ["NT"]:
            local_global_no = "none"
        else:
            raise ValueError(f"The thresholding method {thresholding_method} is not recognized.")

        if thresholding_method in ["GPTL", "GPTG", "LTMMP"]:
            percentile_value_mean = "percentile"
        elif thresholding_method in ["GTV", "LTMMV"]:
            percentile_value_mean = "value"
        elif thresholding_method in ["LTM"]:
            percentile_value_mean = "mean"

        if thresholding_method in ["GPTL"]:
            second_local_global_no = "local"
        elif thresholding_method in ["GPTG"]:
            second_local_global_no = "global"

        gene_ids = expression_data.iloc[1:, 0].tolist()
        sample_ids = expression_data.columns[1:].tolist()

        # TODO add extra check + specific error raising if gene ids or sample ids turn out to be incorrect here

        genes_not_in_model = [gene for gene in gene_ids if gene not in cobra_model.genes]
        expression_data = expression_data[~expression_data.iloc[:, 0].isin(genes_not_in_model)]
        gene_ids = expression_data.iloc[0:, 0].tolist()
        expression_data = expression_data.iloc[0:, 1:]

        if len(sample_ids) > 1:
            linear_data = expression_data.values.flatten()
        else:
            linear_data = expression_data.values

        if local_global_no == "global":
            if percentile_value_mean == "percentile":
                if second_local_global_no == "local":
                    l_global = np.percentile(np.log10(expression_data), thresholding_value, axis=0)
                    thresholds = (10 ** l_global, None)
                elif second_local_global_no == "global":
                    l_global = np.percentile(np.log10(linear_data), thresholding_value)
                    thresholds = (10 ** l_global, None)
            elif percentile_value_mean == "value":
                thresholds = (thresholding_value, None)
            else:
                raise ValueError(f"The thresholding method {thresholding_method} is not recognized.")
        elif local_global_no == "local":
            if thresholding_method == "LTM":
                thresholds = (expression_data.mean(), None)
            elif thresholding_method == "LTMMP":
                l_high = np.percentile(np.log10(linear_data), thresholding_value_2)
                data_ths_high = 10 ** l_high
                l_low = np.percentile(np.log10(linear_data), thresholding_value)
                data_ths_low = 10 ** l_low
                thresholds = (data_ths_low, data_ths_high)
            elif thresholding_method == "LTMMV":
                thresholds = (thresholding_value, thresholding_value_2)
            else:
                raise ValueError(f"The thresholding method {thresholding_method} is not recognized.")
        elif local_global_no == "none":
            thresholds = (1, None)
        else:
            raise ValueError(f"The thresholding method {thresholding_method} is not recognized.")

        if local_global_no == "global" or "none":
            if transformation == "NT":
                gene_score = expression_data / thresholds[0]
            elif transformation == "LNT":
                gene_score = np.log(1 + (expression_data / thresholds[0]))
            elif transformation == "L2T":
                gene_score = np.log2(1 + (expression_data / thresholds[0]))
            elif transformation == "L10T":
                gene_score = np.log10(1 + (expression_data / thresholds[0]))
            elif transformation == "ZTL":
                z_score_expression = zscore_columnwise(expression_data)
                minimum_gene_score = np.min(z_score_expression)
                z_score_expression = z_score_expression + abs(minimum_gene_score) + 0.01
                gene_score = z_score_expression / thresholds[0]
            elif transformation == "ZTG":
                z_score_expression = zscore_dataframe(expression_data)
                minimum_gene_score = np.min(z_score_expression)
                z_score_expression = z_score_expression + abs(minimum_gene_score) + 0.01
                gene_score = z_score_expression / thresholds[0]
        elif local_global_no == "local":
            threshold = [0] * expression_data.shape[0]
            if thresholds[1] is not None:
                for idx in range((expression_data.shape)[0]):
                    expression_value = expression_data.iloc[idx, :]  # Get expression values for gene i
                    mean_expression = expression_value.mean()  # Calculate mean of expression values
                    if mean_expression >= thresholding_value_2:
                        threshold[idx] = thresholding_value_2
                    else:
                        threshold[idx] = max(mean_expression, thresholding_value)
            else:
                threshold = thresholds[0]
            gene_score = np.zeros((expression_data.shape[0], expression_data.shape[1]))
            for i in range(expression_data.shape[1]):  # Loop over samples
                gene_score[:, i] = np.log1p(expression_data.iloc[:, i].values / threshold)

        sample_number = expression_data.shape[1]
        GPR_method = global_options_dict["__global_options"]["GPR_method"]
        GPR_OR_argument = re.findall('[A-Z][^A-Z]*', GPR_method)
        GPR_OR_argument = GPR_OR_argument[0].replace("_", "")
        GPR_AND_argument = GPR_method.split(GPR_OR_argument)[0]
        GPR_value = GPR_method.split("_")[1] if len(GPR_method.split("_")) > 1 else None
        GPR_value_2 = GPR_method.split("_")[2] if len(GPR_method.split("_")) > 2 else None
        GPR_OR_argument = GPR_OR_argument.lower()
        promiscuity_method = global_options_dict["__global_options"]["promiscuity_method"]

        (gene_contribution_matrix, reaction_expression_matrix,
         per_reaction_gene_contribution_tuple_matrix, per_reaction_adjusted_gene_contribution_tuple_matrix ) = CombinationFileCreation.allocate_gene_values_to_reactions(
            cobra_model, gene_score, expression_data,
            gene_ids, sample_number, promiscuity_method,
            GPR_OR_argument, GPR_AND_argument, GPR_value,
            GPR_value_2
        )

        return (reaction_expression_matrix, gene_contribution_matrix, gene_ids,
                per_reaction_gene_contribution_tuple_matrix,
                per_reaction_adjusted_gene_contribution_tuple_matrix)

    @staticmethod
    def create_and_save_dataframes_from_matrices(
            expression_data_location: str,
            global_options_dict: dict,
            genes_ids: list,
            gene_contribution_matrix: np.ndarray,
            reaction_expression_matrix: np.ndarray,
            cobra_model: Model,
            is_irreversible: bool,
    ):
        """
        Create and save the dataframes from the matrices.

        :param gene_contribution_matrix: np.ndarray
            The gene contribution matrix.
        :param reaction_expression_matrix: np.ndarray
            The reaction expression matrix.
        :param cobra_model: Model
            The cobrapy model.
        """

        rev_or_irrev = "irr" if is_irreversible else "rev"
        base_name = (f"{global_options_dict['__global_options']['transformation']}__"
                     f"{global_options_dict['__global_options']['thresholding']}__"
                     f"{global_options_dict['__global_options']['GPR_method']}__"
                     f"{global_options_dict['__global_options']['promiscuity_method']}__"
                     f"{rev_or_irrev}__"
                     "{}_dataset.feather")
        gene_contribution_file_name = base_name.format("gene_contribution")
        reaction_expression_file_name = base_name.format("reaction_expression")

        sample_ids = [f"sample_{i+1}" for i in range(reaction_expression_matrix.shape[1])]
        if gene_contribution_matrix is not None:
            gene_contribution_df = pd.DataFrame(gene_contribution_matrix, index=genes_ids, columns=sample_ids)
            gene_contribution_location = os.path.join(expression_data_location, gene_contribution_file_name)
            gene_contribution_df.to_feather(gene_contribution_location)
            if global_options_dict["__global_options"].get("create_csv_files", False):
                gene_contribution_df.to_csv(gene_contribution_location.replace(".feather", ".csv"))
        else:
            gene_contribution_df = None
        reaction_expression_df = pd.DataFrame(reaction_expression_matrix, index=[reaction.id for reaction in cobra_model.reactions], columns=sample_ids)
        # replace all nan values with -1
        reaction_expression_df.fillna(-1, inplace=True)

        reaction_expression_location = os.path.join(expression_data_location, reaction_expression_file_name)
        reaction_expression_df.to_feather(reaction_expression_location)
        if global_options_dict["__global_options"].get("create_csv_files", False):
            reaction_expression_df.to_csv(reaction_expression_location.replace(".feather", ".csv"))
        samples_location = os.path.join(expression_data_location, "expression_samples")
        os.mkdir(samples_location) if not os.path.exists(samples_location) else None
        base_sample_name = (f"{global_options_dict['__global_options']['transformation']}__"
                            f"{global_options_dict['__global_options']['thresholding']}__"
                            f"{global_options_dict['__global_options']['GPR_method']}__"
                            f"{global_options_dict['__global_options']['promiscuity_method']}__"
                            f"{rev_or_irrev}__"
                            "{}")
        gene_contribution_sample_file_name = base_sample_name.format("gene_contribution_sample_{:03}.csv")
        reaction_expression_sample_file_name = base_sample_name.format("reaction_expression_sample_{:03}.csv")
        for sample in range(reaction_expression_matrix.shape[1]):
            reaction_expression_df_sample = pd.DataFrame(reaction_expression_matrix[:, sample], index=[reaction.id for reaction in cobra_model.reactions],
                                                         columns=["reaction_expression"])
            reaction_expression_df_sample.fillna(-1, inplace=True)
            if gene_contribution_matrix is not None:
                gene_contribution_df_sample = pd.DataFrame(gene_contribution_matrix[:, sample], index=genes_ids, columns=["gene_contribution"])
                gene_contribution_df_location = os.path.join(samples_location, gene_contribution_sample_file_name.format(sample + 1))
                gene_contribution_df_sample.to_csv(gene_contribution_df_location)
            reaction_expression_df_location = os.path.join(samples_location, reaction_expression_sample_file_name.format(sample + 1))
            reaction_expression_df_sample.to_csv(reaction_expression_df_location)

        return gene_contribution_df, reaction_expression_df

    @staticmethod
    def create_and_save_per_reaction_matrices(
            combination_folder: str,
            global_options_dict: dict,
            per_reaction_gene_contribution_tuple_matrix: list[list[dict]],
            per_reaction_adjusted_gene_contribution_tuple_matrix: list[list[dict]],
            cobra_model: Model,
    ):
        base_name = (f"{global_options_dict['__global_options']['transformation']}__"
                     f"{global_options_dict['__global_options']['thresholding']}__"
                     f"{global_options_dict['__global_options']['GPR_method']}__"
                     f"{global_options_dict['__global_options']['promiscuity_method']}__"
                     f"rev__"
                     "{}.parquet")

        per_reaction_gene_contribution_file_name = base_name.format(
            "per_reaction_gene_contribution"
        )
        per_reaction_adjusted_gene_contribution_file_name = base_name.format(
            "per_reaction_adjusted_gene_contribution"
        )
        simple_per_reaction_gene_contribution_file_name = base_name.format(
            "simple_per_reaction_gene_contribution"
        )
        simple_per_reaction_adjusted_gene_contribution_file_name = base_name.format(
            "simple_per_reaction_adjusted_gene_contribution"
        )

        per_reaction_gene_contribution_df = pd.DataFrame(
            per_reaction_gene_contribution_tuple_matrix,
            index=[(reaction.id, reaction.gene_reaction_rule)
                   for reaction in cobra_model.reactions],
            columns=[f"sample_{i+1}"
                     for i in range(len(per_reaction_gene_contribution_tuple_matrix[0]))]
        )
        per_reaction_adjusted_gene_contribution_df = pd.DataFrame(
            per_reaction_adjusted_gene_contribution_tuple_matrix,
            index=[(reaction.id, reaction.gene_reaction_rule)
                   for reaction in cobra_model.reactions],
            columns=[f"sample_{i+1}"
                     for i in range(len(per_reaction_adjusted_gene_contribution_tuple_matrix[0]))]
        )
        simple_per_reaction_gene_contribution = deepcopy(
            per_reaction_gene_contribution_tuple_matrix
        )
        simple_per_reaction_adjusted_gene_contribution = deepcopy(
            per_reaction_adjusted_gene_contribution_tuple_matrix
        )

        for idx in range(len(per_reaction_gene_contribution_tuple_matrix)):
            for idx2 in range(len(per_reaction_gene_contribution_tuple_matrix[idx])):
                simple_per_reaction_gene_contribution[idx][idx2] = {
                    key: values for key, values in per_reaction_gene_contribution_tuple_matrix[idx][idx2].items()
                    if values[1] > 0
                }
                simple_per_reaction_adjusted_gene_contribution[idx][idx2] = {
                    key: values for key, values in per_reaction_adjusted_gene_contribution_tuple_matrix[idx][idx2].items()
                    if values[1] > 0
                }
        simple_per_reaction_gene_contribution_df = pd.DataFrame(
            simple_per_reaction_gene_contribution,
            index=[(reaction.id, reaction.gene_reaction_rule)
                   for reaction in cobra_model.reactions],
            columns=[f"sample_{i+1}"
                     for i in range(len(simple_per_reaction_gene_contribution[0]))]
        )
        simple_per_reaction_adjusted_gene_contribution_df = pd.DataFrame(
            simple_per_reaction_adjusted_gene_contribution,
            index=[(reaction.id, reaction.gene_reaction_rule)
                   for reaction in cobra_model.reactions],
            columns=[f"sample_{i+1}"
                     for i in range(len(simple_per_reaction_adjusted_gene_contribution[0]))]
        )

        per_reaction_gene_contribution_location = os.path.join(combination_folder, per_reaction_gene_contribution_file_name)
        per_reaction_adjusted_gene_contribution_location = os.path.join(combination_folder, per_reaction_adjusted_gene_contribution_file_name)
        simple_per_reaction_gene_contribution_location = os.path.join(combination_folder, simple_per_reaction_gene_contribution_file_name)
        simple_per_reaction_adjusted_gene_contribution_location = os.path.join(combination_folder, simple_per_reaction_adjusted_gene_contribution_file_name)

        def process_and_save_df(df, file_location):
            df = df.map(lambda x: str(x) if isinstance(x, dict) else x)
            df = df.map(lambda x: str(x) if isinstance(x, tuple) else x)
            df = df.map(lambda x: str(x) if isinstance(x, dict) else x)
            df = df.map(lambda x: str(x) if isinstance(x, tuple) else x)
            df = df.astype(str)  # Ensure all columns are strings
            df.to_parquet(file_location)

        process_and_save_df(per_reaction_gene_contribution_df, per_reaction_gene_contribution_location)
        process_and_save_df(per_reaction_adjusted_gene_contribution_df, per_reaction_adjusted_gene_contribution_location)
        process_and_save_df(simple_per_reaction_gene_contribution_df, simple_per_reaction_gene_contribution_location)
        process_and_save_df(simple_per_reaction_adjusted_gene_contribution_df, simple_per_reaction_adjusted_gene_contribution_location)


        # do the same for csv
        if global_options_dict["__global_options"].get("create_csv_files", False):
            per_reaction_gene_contribution_df.to_csv(
                per_reaction_gene_contribution_location.replace(
                    ".parquet", ".csv")
            )
            per_reaction_adjusted_gene_contribution_df.to_csv(
                per_reaction_adjusted_gene_contribution_location.replace(
                    ".parquet", ".csv")
            )
            simple_per_reaction_gene_contribution_df.to_csv(
                simple_per_reaction_gene_contribution_location.replace(
                    ".parquet", ".csv")
            )
            simple_per_reaction_adjusted_gene_contribution_df.to_csv(
                simple_per_reaction_adjusted_gene_contribution_location.replace(
                    ".parquet", ".csv")
            )

    @staticmethod
    def transform_reversible_expression_to_irreversible(
            reaction_expression_matrix: np.ndarray,
            rev_to_irrev_list: list,
    ) -> np.ndarray:
        """
        Transform the reversible expression matrix to irreversible.

        :param reaction_expression_matrix: np.ndarray
            The reaction expression matrix.
        """
        size_ = max(x[1] for x in rev_to_irrev_list if len(x) > 1)
        irreversible_expression_matrix = np.zeros((size_, reaction_expression_matrix.shape[1]))

        for idx, reaction in enumerate(rev_to_irrev_list):
            irreversible_expression_matrix[idx, :] = reaction_expression_matrix[idx, :]
            if len(rev_to_irrev_list[idx]) > 1:
                irreversible_expression_matrix[rev_to_irrev_list[idx][1]-1, :] = reaction_expression_matrix[rev_to_irrev_list[idx][0]-1, :]
        return irreversible_expression_matrix

    @staticmethod
    def create_taskstructure_from_tasklist(
            tasklist_location: str,
            cobra_model: Model,
    ) -> dict:
        """
        Create the task structure from the tasklist.

        :param tasklist_location: str
            The path to the tasklist.
        """
        value_if_no_bounds = "unbounded"
        tasklist = None
        for file in os.listdir(tasklist_location):
            if file.startswith("tasklist"):
                if file.endswith(".xlsx"):
                    tasklist = pd.read_excel(os.path.join(tasklist_location, file))
                    break
                elif file.endswith(".csv"):
                    tasklist = pd.read_csv(os.path.join(tasklist_location, file))
                    break

        if tasklist is None:
            raise ValueError(f"The tasklist file in {tasklist_location} is not in a readable format.")

        individual_task_structure = {
            "task_id": int,
            "task_description": None,
            "task_should_fail_in_humans": None,
            "input_metabolites": [],
            "input_metabolite_lower_bounds": [],
            "input_metabolite_upper_bounds": [],
            "output_metabolites": [],
            "output_metabolite_lower_bounds": [],
            "output_metabolite_upper_bounds": [],
        }
        task_structure = {}
        for idx in range(len(tasklist)):
            task_id = tasklist.iloc[idx, 0]
            if not pd.isna(task_id):
                task_id = int(task_id)
                # from row of the current task scan all rows that are na for task id or we reach the end of the file
                current_task_structure = deepcopy(individual_task_structure)
                current_task_structure["task_id"] = task_id
                for idx2 in range(idx, len(tasklist)):
                    if idx2 == idx or pd.isna(tasklist.iloc[idx2, 0]):
                        if not pd.isna(tasklist.iloc[idx2, 1]):
                            current_task_structure["task_description"] = tasklist.iloc[idx2, 1]
                        if not pd.isna(tasklist.iloc[idx2, 2]):
                            current_task_structure["task_should_fail_in_humans"] = int(tasklist.iloc[idx2, 2])
                        if not pd.isna(tasklist.iloc[idx2, 3]):
                            current_task_structure["input_metabolites"].append(tasklist.iloc[idx2, 3])
                            if not pd.isna(tasklist.iloc[idx2, 4]):
                                current_task_structure["input_metabolite_lower_bounds"].append(int(tasklist.iloc[idx2, 4]))
                            else:
                                current_task_structure["input_metabolite_lower_bounds"].append(value_if_no_bounds)
                            if not pd.isna(tasklist.iloc[idx2, 5]):
                                current_task_structure["input_metabolite_upper_bounds"].append(int(tasklist.iloc[idx2, 5]))
                            else:
                                current_task_structure["input_metabolite_upper_bounds"].append(value_if_no_bounds)
                        if not pd.isna(tasklist.iloc[idx2, 6]):
                            current_task_structure["output_metabolites"].append(tasklist.iloc[idx2, 6])
                            if not pd.isna(tasklist.iloc[idx2, 7]):
                                current_task_structure["output_metabolite_lower_bounds"].append(int(tasklist.iloc[idx2, 7]))
                            else:
                                current_task_structure["output_metabolite_lower_bounds"].append(value_if_no_bounds)
                            if not pd.isna(tasklist.iloc[idx2, 8]):
                                current_task_structure["output_metabolite_upper_bounds"].append(int(tasklist.iloc[idx2, 8]))
                            else:
                                current_task_structure["output_metabolite_upper_bounds"].append(value_if_no_bounds)
                    elif idx2 != idx and not pd.isna(tasklist.iloc[idx2, 0]):
                        break
                if not all([x.startswith("MAM") for x in current_task_structure["input_metabolites"]]):
                    for idx3, met in enumerate(current_task_structure["input_metabolites"]):
                        if not met.startswith("MAM"):
                            replaced = False
                            met_without_compartment = "[".join(met.split("[")[:-1]).lower()
                            compartment = "".join(met.split("[")[-1]).lower().split("]")[0]
                            for modelmet in cobra_model.metabolites:
                                if met_without_compartment == modelmet.name.lower() and compartment == modelmet.compartment.lower():
                                    replaced = True
                                    current_task_structure["input_metabolites"][idx3] = modelmet.id
                                    break
                            if not replaced:
                                raise ValueError(f"The input metabolite {met} of task {task_id} is not present in the model.")
                if not all([x.startswith("MAM") for x in current_task_structure["output_metabolites"]]):
                    for idx3, met in enumerate(current_task_structure["output_metabolites"]):
                        if not met.startswith("MAM"):
                            replaced = False
                            met_without_compartment = "[".join(met.split("[")[:-1]).lower()
                            compartment = "".join(met.split("[")[-1]).lower().split("]")[0]
                            for modelmet in cobra_model.metabolites:
                                if met_without_compartment == modelmet.name.lower() and compartment == modelmet.compartment.lower():
                                    replaced = True
                                    current_task_structure["output_metabolites"][idx3] = modelmet.id
                                    break
                            if not replaced:
                                # # "MAM02053[c]"
                                # print(met)
                                # print(tempmet.name)
                                tempmet = None
                                # for modelmet in cobra_model.metabolites:
                                #     if modelmet.id == "MAM00198[l]":
                                #         tempmet = modelmet
                                #         print(tempmet.name)
                                #         print(met_without_compartment)
                                #         break
                                raise ValueError(f"The output metabolite {met} of task {task_id} is not present in the model.")

                task_structure[task_id] = current_task_structure

        return task_structure

    @staticmethod
    def save_task_structure(
            task_structure: dict,
            task_name: str,
            combination_folder: str,
    ):
        """
        Save the task structure to a file.

        :param task_structure: dict
            The task structure.
        :param combination_folder: str
            The path to the combination folder.
        """
        task_structure_location = os.path.join(combination_folder, f"task_structure_{task_name}.json")
        with open(task_structure_location, "w") as file:
            dump(task_structure, file)

    @staticmethod
    def populate_combination_files(
            main_data_folder: str,
            global_options_dict: dict,
            inputs_dict: dict,
            names_dict_by_options: dict,
            combination_folder: str,

    ) -> dict:
        """
        Populate the combination files with the necessary files for the preprocessing.
            - Creates a closed base model in the "models" directory.

        :param main_data_folder: str
            The main directory where the data is stored.
        :param global_options_dict: dict
            The dictionary with the global options for the preprocessing.
        :param inputs_dict: dict
            The dictionary with the inputs for the preprocessing.
        :param names_dict_by_options: dict
            The dictionary with the names of the files based on the global options.
        """
        selected_model = global_options_dict["__global_options"]["model_version"]
        selected_expression_dataset = global_options_dict["__global_options"]["expression_dataset"]
        selected_tasklist = global_options_dict["__global_options"]["tasklist"]

        if inputs_dict["model"] is None:
            model_location = os.path.join(main_data_folder, "models", selected_model)
            FileChecker.raise_error_if_incorrect_folders(model_location)
            cobra_model = FileChecker.read_model_file_in_folder(None, model_location, None)
            cobra_model_reversible, rev_to_irrev_dict = CombinationFileCreation.create_closed_model(
                cobra_model,
                False
            )
            # cobra_model = FileChecker.read_model_file_in_folder(None, model_location, None)
            cobra_model_irreversible, rev_to_irrev_dict = CombinationFileCreation.create_closed_model(
                cobra_model,
                True
            )
            inputs_dict["rev_to_irrev_dict"] = rev_to_irrev_dict
            rev_to_irrev_name = f"{global_options_dict['__global_options']['model_version']}__rev_to_irrev_dict.json"
            model_name_by_options_rev = (
                f"model_{global_options_dict['__global_options']['model_version']}__reversible.json"
            )
            model_name_by_options_irrev = (
                f"model_{global_options_dict['__global_options']['model_version']}__irreversible.json"
            )
            save_json_model(cobra_model_reversible, os.path.join(combination_folder, model_name_by_options_rev))
            save_json_model(cobra_model_irreversible, os.path.join(combination_folder, model_name_by_options_irrev))
            with open(os.path.join(combination_folder, rev_to_irrev_name), "w") as file:
                dump(rev_to_irrev_dict, file)
            if global_options_dict['__global_options']['model_type'] == "rev":
                inputs_dict["model"] = cobra_model_reversible
                cobra_model = cobra_model_reversible
            elif global_options_dict['__global_options']['model_type'] == "irr":
                inputs_dict["model"] = cobra_model_irreversible
                cobra_model = cobra_model_irreversible
        else:
            cobra_model = inputs_dict["model"]

        if inputs_dict["rev_to_irrev_dict"] is None:
            cobra_model, rev_to_irrev_dict = CombinationFileCreation.create_closed_model(cobra_model, global_options_dict)
            inputs_dict["model"] = cobra_model
            inputs_dict["rev_to_irrev_dict"] = rev_to_irrev_dict
            rev_to_irrev_name = f"{global_options_dict['__global_options']['model_version']}__rev_to_irrev_dict.json"
            save_json_model(rev_to_irrev_dict, os.path.join(combination_folder, rev_to_irrev_name))
        else:
            rev_to_irrev_dict = inputs_dict["rev_to_irrev_dict"]

        if inputs_dict["expression_data"] is None:
            expression_data_location = os.path.join(
                main_data_folder,
                "expression_datasets",
                selected_expression_dataset
            )
            FileChecker.raise_error_if_incorrect_folders(expression_data_location)
            (
                GPR_expression_data,
                GPR_contribution_data,
                genes_ids,
                per_reaction_gene_contribution_tuple_matrix,
                per_reaction_adjusted_gene_contribution_tuple_matrix
            ) = CombinationFileCreation.create_GPR_expression_data(
                    expression_data_location,
                    cobra_model_reversible,
                    global_options_dict,
                    selected_model,
            )
            (
                GPR_expression_data_df,
                GPR_contribution_data_df
            ) = CombinationFileCreation.create_and_save_dataframes_from_matrices(
                    combination_folder,
                    global_options_dict,
                    genes_ids,
                    GPR_contribution_data,
                    GPR_expression_data,
                    cobra_model_reversible,
                    False
            )
            CombinationFileCreation.create_and_save_per_reaction_matrices(
                combination_folder,
                global_options_dict,
                per_reaction_gene_contribution_tuple_matrix,
                per_reaction_adjusted_gene_contribution_tuple_matrix,
                cobra_model_reversible
            )
            (
                irreversible_GPR_expression_data
            ) = CombinationFileCreation.transform_reversible_expression_to_irreversible(
                    GPR_expression_data,
                    rev_to_irrev_dict
            )
            (
                irreversible_GPR_expression_data_df,
                GPR_contribution_data
            ) = CombinationFileCreation.create_and_save_dataframes_from_matrices(
                    combination_folder,
                    global_options_dict,
                    genes_ids,
                    None,
                    irreversible_GPR_expression_data,
                    cobra_model_irreversible,
                    True
            )
            if global_options_dict['__global_options']['model_type'] == "rev":
                inputs_dict["expression_data"] = GPR_expression_data_df
            elif global_options_dict['__global_options']['model_type'] == "irr":
                inputs_dict["expression_data"] = irreversible_GPR_expression_data_df

        if inputs_dict["tasklist"] is None:
            tasklist_location = os.path.join(main_data_folder, "tasklists", selected_tasklist)
            FileChecker.raise_error_if_incorrect_folders(tasklist_location)
            FileChecker.check_if_tasklist_file_is_correct(tasklist_location)
            tasklist = CombinationFileCreation.create_taskstructure_from_tasklist(tasklist_location, cobra_model)
            CombinationFileCreation.save_task_structure(tasklist, selected_tasklist, combination_folder)
            inputs_dict["tasklist"] = tasklist

        else:
            tasklist = inputs_dict["tasklist"]

        return inputs_dict

    @staticmethod
    def master_create_combination_directory(
            main_data_folder: str,
            global_options_dict: dict,
    ):
        """
        Create the "combinations" directory and populate it with the necessary files for the preprocessing.
            - Creates a closed base model in the "models" directory.

        :param main_data_folder: str
            The main directory where the data is stored.
        :param global_options_dict: dict
            The dictionary with the global options for the preprocessing.
        """
        combinations_folder = os.path.join(main_data_folder, "combinations")
        os.mkdir(combinations_folder) if not os.path.exists(combinations_folder) else None
        combination_folder = CombinationFileCreation.create_combination_directory(main_data_folder, global_options_dict)

        # check if files exist
        (check_if_all_correct_dict, names_dict_by_options) = FileChecker.check_if_combination_files_exist(
            combination_folder, global_options_dict
        )
        # check if files correct
        if (global_options_dict['__global_options']["check_if_input_files_correct"] and
                all(value is True for value in check_if_all_correct_dict.values())):
            (check_if_all_correct_dict, inputs_dict) = FileChecker.check_if_combination_files_are_correct(
                check_if_all_correct_dict, names_dict_by_options, global_options_dict
            )
        else:
            inputs_dict = {
                "expression_data": None,
                "model": None,
                "tasklist": None,
            }
        # if all files exist, check them for correctness, if not create them and check
        # if after creation they are still incorrect, do not continue
        if all(value is True for value in check_if_all_correct_dict.values()):
            (check_if_all_correct_dict, inputs_dict) = FileChecker.check_if_combination_files_are_correct(
                check_if_all_correct_dict, names_dict_by_options,
                global_options_dict, inputs_dict
            )
            if all(value is True for value in check_if_all_correct_dict.values()):
                return inputs_dict
            else:
                CombinationFileCreation.populate_combination_files(
                    main_data_folder, global_options_dict, inputs_dict, names_dict_by_options, combination_folder
                )
                (check_if_all_correct_dict, inputs_dict) = FileChecker.check_if_combination_files_are_correct(
                    check_if_all_correct_dict, names_dict_by_options,
                    global_options_dict, inputs_dict
                )
                if all(value is True for value in check_if_all_correct_dict.values()):
                    return inputs_dict
                else:
                    raise FileNotFoundError(f"An error occurred while creating the combination files."
                                            f"Or the created files did not pass the test.")

        elif not all(value is True for value in check_if_all_correct_dict.values()):
            CombinationFileCreation.populate_combination_files(
                main_data_folder, global_options_dict, inputs_dict, names_dict_by_options, combination_folder
            )
            check_if_all_correct_dict = {key: True for key in check_if_all_correct_dict.keys()}
            (check_if_all_correct_dict, inputs_dict) = FileChecker.check_if_combination_files_are_correct(
                check_if_all_correct_dict, names_dict_by_options,
                global_options_dict, inputs_dict
            )
            if all(value is True for value in check_if_all_correct_dict.values()):
                return inputs_dict
            else:
                raise FileNotFoundError(f"An error occurred while creating the combination files. "
                                        f"Or the created files did not pass the test.")

class IndividualFileUpdating:
    @staticmethod
    def update_expression_data(
    ):
        pass

    @staticmethod
    def update_model(
    ):
        pass

    @staticmethod
    def update_taskstructure(
    ):
        pass

    @staticmethod
    def update_rev_to_irrev(

    ):
        pass


if __name__ == "__main__":
    start_task =2
    end_task = 2

    methods_dict_for_improved = {
        "__global_options": {
            "epsilon": 1e-4,
            "expressionless_reaction_value": 1.14,
            "create_log_file": True,
            "tasks": list(range(int(start_task), int(end_task + 1))),
            "samples": [1,2,3,4],  # list(range(1, 1010)),
            "unbounded_lower_bound_value": 1,
            "unbounded_upper_bound_value": 1000,
            "solver_settings": {
                "TimeLimit": 2200,
                # "IntFeasTol": 1e-7,
                # "FeasibilityTol": 1e-7,
                # "MIPFocus": 2,
                # "Heuristics": 0.01,
                # "Cuts": 3,
                # "Presolve": 2,
                # "ScaleFlag": 3,
            },
            "transformation": "NT",
            "thresholding": "NT",
            "GPR_method": "min2Sum",
            "promiscuity_method": "N",
            "check_if_input_files_correct": True,
            "model_type": "irr",
            "tasklist": "paracetamol_final_final",
            "model_version": "updated_paracetamol_specific",
            "expression_dataset": "liver_kidney",
            "alternative_reaction_notation": False,
            "result_output_folder": r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\results",
            "name_of_run": "paracetamol_final_final_hepato",
            "create_csv_files": True,
            "overwrite_run_files": True,
            # if True, will not try to create new file names but instead overwrite existing files
            "skip_existing_files": False,  # if True, will skip the task and sample if the files already exist
            "version_of_code": version,
            "ignore_sample_warnings": True,
            "minimize_filtered_model": True,
        },
    }
    main_data_folder = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files"
    inputs_dict = CombinationFileCreation.master_create_combination_directory(
        main_data_folder,
        methods_dict_for_improved
    )
    print('done with combinations folder')

    # matlab_file = r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.mat"
    #
    #
    # # model = load_matlab_model(matlab_file)
    # # save_json_model(model, r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.json")
    # # write_sbml_model(model, r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.xml")
    # # save_matlab_model(model, r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.mat")
    #
    # # test performance of reading the three types of model
    # def performance_test_read():
    #     for i in range(10):
    #         start = time.time()
    #         load_matlab_model(matlab_file)
    #         print(f"Matlab: {time.time() - start}")
    #         start = time.time()
    #         load_json_model(r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.json")
    #         print(f"Json: {time.time() - start}")
    #         start = time.time()
    #         read_sbml_model(r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.xml")
    #         print(f"SBML: {time.time() - start}")
    #
    #
    # # cProfile.run("performance_test_read()")
    # cProfile.run("load_json_model(r'E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\consensus_2024\model_Human-GEM_COBRA version 17.json')")
