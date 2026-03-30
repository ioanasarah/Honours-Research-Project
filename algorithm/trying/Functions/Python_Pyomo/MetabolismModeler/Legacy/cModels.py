import glob
import os
import re
import logging
import warnings
import scipy.io
import cobra
import numpy as np
import json
import time
import datetime
import copy
import pyomo.environ as pe
import pyomo.opt as po
from .sResults import get_active_reactions_from_solved_model, save_solver_results


__all__ = [
    "load_irreversible_model_and_expression_from_mat_files_legacy",
    "load_reversible_model_and_expression_from_mat_files_legacy",
    "run_test_suite_different_model_formulations_legacy",
    "create_basic_pyomo_model_SV_constraints",
]


def create_basic_pyomo_model_SV_constraints(
    stoichiometric_matrix: np.ndarray,
    cobra_model: cobra.Model,
    using_cobra: bool = True,
) -> (pe.ConcreteModel, dict):
    model = pe.ConcreteModel()
    num_metabolites = len(stoichiometric_matrix)  # Example number of metabolites
    num_reactions = len(stoichiometric_matrix[0])  # Example number of reactions

    model.Metabolites = pe.RangeSet(1, num_metabolites)
    model.Rxns = pe.RangeSet(1, num_reactions)
    stoich_dict = {
        met: {
            rxn: stoichiometric_matrix[met, rxn]
            for rxn in range(num_reactions)
            if stoichiometric_matrix[met, rxn] != 0
        }
        for met in range(num_metabolites)
    }
    model.V_flux = pe.Var(model.Rxns, domain=pe.Reals)
    for rxn in model.Rxns:
        model.V_flux[rxn].setlb(cobra_model.reactions[rxn - 1].lower_bound)
        model.V_flux[rxn].setub(cobra_model.reactions[rxn - 1].upper_bound)

    def sv_zero_rule(model, met):
        return (
            sum(
                stoich_dict[met - 1][rxn] * model.V_flux[rxn + 1]
                for rxn in stoich_dict[met - 1]
            )
            == 0
        )

    model.SV_constraint = pe.Constraint(model.Metabolites, rule=sv_zero_rule)

    return model, stoich_dict


def load_irreversible_model_and_expression_from_mat_files_legacy(
    task_number: int, main_folder: str = "", directory_path: str = ""
) -> (cobra.Model, list, list):
    # if main_folder is "" then use current directory and search for task_number if directory_path is ""
    # if diretor_path is not "" then use that directory
    # if main_folder != "": then search there for task number
    filenames_to_get = [
        "expValueMeanAdjustedTaskIrrevKapprox.mat",
        "newIrrevModelKapprox.mat",
        "rev2irrevIrrevKapprox.mat",
    ]
    model = None
    expression = None
    rev2irrev = None
    if directory_path == "":
        if main_folder == "":
            folder_to_search = os.getcwd()
        else:
            folder_to_search = main_folder
        all_dirs = [dirpath for dirpath, dirnames, files in os.walk(folder_to_search)]
        pattern = rf"Task0*{task_number}S"
        for dir_ in all_dirs:
            if re.search(pattern, dir_):
                directory_path = dir_
                break
    print("Loading files from: ", directory_path)
    for file_to_get in filenames_to_get:
        search_pattern = os.path.join(directory_path, "**", file_to_get)
        files = glob.glob(search_pattern, recursive=True)
        if files:
            for file in files:
                if file_to_get == "newIrrevModelKapprox.mat":
                    with warnings.catch_warnings():
                        original_logging_level = logging.getLogger().getEffectiveLevel()
                        logging.getLogger().setLevel(logging.ERROR)
                        warnings.filterwarnings(
                            "ignore", "No defined compartments in model"
                        )
                        model = cobra.io.load_matlab_model(file)
                        logging.getLogger().setLevel(original_logging_level)

                elif file_to_get == "expValueMeanAdjustedTaskIrrevKapprox.mat":
                    expression_cells = scipy.io.loadmat(file)
                    expression = [
                        value[0]
                        for value in expression_cells[
                            "expValueMeanAdjustedTaskIrrevKapprox"
                        ]
                    ]

                elif file_to_get == "rev2irrevIrrevKapprox.mat":
                    rev2irrev_cells = scipy.io.loadmat(file)
                    rev2irrev = [
                        value[0][0]
                        for value in rev2irrev_cells["rev2irrevIrrevKapprox"]
                    ]
    if model is None or expression is None or rev2irrev is None:
        if model is None:
            print("Model was not found")
        if expression is None:
            print("Expression was not found")
        if rev2irrev is None:
            print("rev2irrev was not found")
        raise ValueError("One of the necessary files was not found")
    else:
        print("Model, expression and rev2irrev loaded successfully")
    return model, expression, rev2irrev


if __name__ == "__main__":
    placeholder = 5
    # main_folder_1 = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing"
    # work_model, expression, rev2irrev = load_irreversible_model_and_expression_from_mat_files(1, main_folder_1)
    # # print(model)
    # # print(expression)
    # # print(rev2irrev)


def load_reversible_model_and_expression_from_mat_files_legacy(
    task_number: int, main_folder: str = "", directory_path: str = ""
) -> (cobra.Model, list):
    # if main_folder is "" then use current directory and search for task_number if directory_path is ""
    # if diretor_path is not "" then use that directory
    # if main_folder != "": then search there for task number
    filenames_to_get = [
        "expValueMeanAdjustedTaskReversibleKApproximation.mat",
        "newReversibleModelForKApproxmiation.mat",
    ]
    model = None
    expression = None
    if directory_path == "":
        if main_folder == "":
            folder_to_search = os.getcwd()
        else:
            folder_to_search = main_folder
        all_dirs = [dirpath for dirpath, dirnames, files in os.walk(folder_to_search)]
        pattern = rf"Task0*{task_number}S"
        for dir_ in all_dirs:
            if re.search(pattern, dir_):
                directory_path = dir_
                break
    print("Loading files from: ", directory_path)
    for file_to_get in filenames_to_get:
        search_pattern = os.path.join(directory_path, "**", file_to_get)
        files = glob.glob(search_pattern, recursive=True)
        if files:
            for file in files:
                if file_to_get == "newReversibleModelForKApproxmiation.mat":
                    with warnings.catch_warnings():
                        original_logging_level = logging.getLogger().getEffectiveLevel()
                        logging.getLogger().setLevel(logging.ERROR)
                        warnings.filterwarnings(
                            "ignore", "No defined compartments in model"
                        )
                        model = cobra.io.load_matlab_model(file)
                        logging.getLogger().setLevel(original_logging_level)

                elif (
                    file_to_get
                    == "expValueMeanAdjustedTaskReversibleKApproximation.mat"
                ):
                    expression_cells = scipy.io.loadmat(file)
                    expression = [
                        value[0]
                        for value in expression_cells[
                            "expValueMeanAdjustedTaskReversibleKApproximation"
                        ]
                    ]

    if model is None or expression is None:
        if model is None:
            print("Model was not found")
        if expression is None:
            print("Expression was not found")
        raise ValueError("One of the necessary files was not found")
    else:
        print("Model, expression loaded successfully")
    return model, expression


def run_test_suite_different_model_formulations_legacy(
    name_of_run: str,
    test_array: list,
    sample_array: list,
    methods_dict: dict,
    main_data_folder: str,
    main_output_folder_original: str,
) -> None:
    ### Legacy function works with models created in the matlab code, use the non legacy versions for new models
    ### Main differences are the manner by which expression json files are used. Legacy functions load the matfiles
    ### and extract the expression values from there. New functions use base expression json files directly
    ### to save space, using base expression and creating task specific expression is done in the loading of the models
    ### specifically by calculating the differences of amount of reactions and the rev2irrev files

    # methods_dict contains version_names as keys and functions as values
    # test_array contains the test numbers to run
    # users create functions which can be added to the dict and ran
    solver = po.SolverFactory("gurobi")
    results_dict = {}
    results_dict["run_name"] = name_of_run
    results_dict["test_array"] = test_array
    results_dict["sample_array"] = sample_array
    results_dict["run_start_time_date"] = str(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    main_output_folder = os.path.join(main_output_folder_original, name_of_run)
    os.makedirs(main_output_folder, exist_ok=True)

    for version_name, tup in methods_dict.items():
        results_dict[version_name] = {}
        for sample_number in sample_array:
            results_dict[version_name][sample_number] = {
                "model_creation_time_per_task": {},
                "solving_time_per_task": {},
                "reactions_included_per_task": {},
                "additional_settings_per_task": {},
            }

    for task_number in test_array:
        cobra_model_reversible, expression_reversible_ori = (
            load_reversible_model_and_expression_from_mat_files_legacy(
                task_number, main_data_folder
            )
        )
        (
            cobra_model_irreversible,
            expression_irreversible_ori,
            rev2irrev_irreversible,
        ) = load_irreversible_model_and_expression_from_mat_files_legacy(
            task_number, main_data_folder
        )
        for sample_number in sample_array:
            if (
                np.size(expression_reversible_ori[0], 1) > sample_number
                and np.size(expression_irreversible_ori[0], 1) > sample_number
            ):
                expression_reversible = [
                    expression[sample_number]
                    for expression in expression_reversible_ori
                ]
                expression_irreversible = [
                    expression[sample_number]
                    for expression in expression_irreversible_ori
                ]
            else:
                expression_reversible = copy.deepcopy(expression_reversible_ori)
                expression_irreversible = copy.deepcopy(expression_irreversible_ori)

            for version_name, tup in methods_dict.items():
                results_dict[version_name]["additional_settings_per_task"][
                    task_number
                ] = {}
                model = None
                chosen_expression = None
                chosen_expression_old = None
                chosen_model = None

                type_of_cobra_model = tup[0]
                function_to_run = tup[1]
                settings_kwargs_dict = tup[2]

                if type_of_cobra_model == "irreversible":
                    chosen_model = cobra_model_irreversible
                    chosen_expression_old = expression_irreversible
                    time_start = time.time()
                    model, chosen_expression = function_to_run(
                        chosen_model,
                        chosen_expression_old,
                        rev2irrev_irreversible,
                        expression_reversible,
                        **settings_kwargs_dict,
                    )
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = (time.time() - time_start)
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = chosen_expression[
                        expression_irreversible.index(-1)
                    ]
                elif type_of_cobra_model == "reversible":
                    chosen_model = cobra_model_reversible
                    chosen_expression_old = expression_reversible
                    time_start = time.time()
                    model, chosen_expression = function_to_run(
                        cobra_model_reversible,
                        expression_reversible,
                        **settings_kwargs_dict,
                    )
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = (time.time() - time_start)
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = chosen_expression[
                        expression_reversible.index(-1)
                    ]
                else:
                    warnings.warn("Invalid type of cobra model")
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = -999
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = -999
                if model is None:
                    warnings.warn(f"Model for task {task_number} was not created")
                    results_dict[version_name]["solving_time_per_task"][sample_number][
                        task_number
                    ] = -999
                    continue
                else:
                    instance = model.create_instance()
                    result = solver.solve(instance, tee=True)
                    results_dict[version_name]["solving_time_per_task"][sample_number][
                        task_number
                    ] = result.Solver.Time
                    if "epsilon" in settings_kwargs_dict:
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["epsilon"] = settings_kwargs_dict["epsilon"]

                    if result.Solver.Termination_condition == "optimal":
                        current_reactions, (
                            amount_of_reactions,
                            amount_of_expressionless_reactions,
                            amount_of_expression_carrying_reactions,
                        ) = get_active_reactions_from_solved_model(
                            instance,
                            chosen_expression_old,
                            chosen_expression,
                            chosen_model,
                            **settings_kwargs_dict,
                        )
                        results_dict[version_name]["reactions_included_per_task"][
                            sample_number
                        ][task_number] = current_reactions
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["amount_of_reactions"] = amount_of_reactions
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number][
                            "amount_of_expressionless_reactions"
                        ] = amount_of_expressionless_reactions
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number][
                            "amount_of_expression_carrying_reactions"
                        ] = amount_of_expression_carrying_reactions

                        save_solver_results(
                            instance,
                            result,
                            chosen_model,
                            f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}",
                            main_output_folder,
                            chosen_expression,
                            chosen_expression_old,
                        )
                    else:
                        warnings.warn(
                            f"Task {task_number} for {version_name} did not converge"
                        )
                        results_dict[version_name]["reactions_included_per_task"][
                            task_number
                        ][sample_number] = {}
                        results_dict[version_name]["additional_settings_per_task"][
                            task_number
                        ][sample_number]["amount_of_reactions"] = -999
                        results_dict[version_name]["additional_settings_per_task"][
                            task_number
                        ][sample_number]["amount_of_expressionless_reactions"] = -999
                        results_dict[version_name]["additional_settings_per_task"][
                            task_number
                        ][sample_number][
                            "amount_of_expression_carrying_reactions"
                        ] = -999

    results_dict["run_end_time_date"] = str(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    print(results_dict)
    print("")
    location_to_save = os.path.join(main_output_folder, f"{name_of_run}_results")
    # to json
    counter = 1
    filename_without_extension = os.path.splitext(location_to_save)[0]
    saved_file = False
    try:
        filename = f"{filename_without_extension}.json"
        if not os.path.exists(filename):
            with open(filename, "w") as f:
                json.dump(results_dict, f)
            saved_file = True
    except Exception as e:
        print(e)
    while not saved_file:
        try:
            filename_with_counter = f"{filename_without_extension}_Try{counter}.json"
            if not os.path.exists(filename_with_counter):
                with open(filename_with_counter, "w") as f:
                    json.dump(results_dict, f)
                break
            else:
                counter += 1
        except Exception as e:
            print(e)
            break
