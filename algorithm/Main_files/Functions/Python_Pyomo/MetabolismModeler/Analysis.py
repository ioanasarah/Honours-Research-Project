from typing import List, Dict, Tuple
from datetime import datetime
from warnings import warn
from json import dump
import numpy as np
import pandas as pd
from copy import deepcopy
from time import time
from os import path, makedirs, listdir
from pyomo.opt import SolverFactory
from . import cModels as cM
from . import sResults as sR
from cobra.io import save_json_model, load_json_model
from cobra import Model
from re import findall
import gurobi_logtools as glt_gurb



def create_empty_dictionary_for_saving_results(
    name_of_run: str,
    test_array: list,
    sample_array: list,
    methods_dict: dict,
) -> dict:
    results_dict = {}
    results_dict["run_name"] = name_of_run
    results_dict["test_array"] = test_array
    results_dict["sample_array"] = sample_array
    global_options_dict = methods_dict.get("__global_options", {})
    results_dict["expression_type"] = global_options_dict.get("expression_type", "local_threshold")
    results_dict["expression_dataset"] = global_options_dict.get("expression_dataset", "DCM")
    results_dict["model_version"] = global_options_dict.get("model_version", "consensus_2024")

    results_dict["run_start_time_date"] = str(
        datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    results_dict["version_name"] = methods_dict.keys()[1] if len(methods_dict) > 1 else methods_dict.keys()[0]
    for sample_number in sample_array:
        results_dict["model_creation_time_per_task"][
            sample_number
        ] = {}
        results_dict["solving_time_per_task"][sample_number] = {}
        results_dict["work_units_per_task"][sample_number] = {}
        results_dict["reactions_included_per_task"][
            sample_number
        ] = {}
        results_dict["additional_settings_per_task"][
            sample_number
        ] = {}
        for task_number in test_array:
            results_dict["model_creation_time_per_task"][
                sample_number
            ][task_number] = {}
            results_dict["solving_time_per_task"][sample_number][
                task_number
            ] = {}
            results_dict["work_units_per_task"][sample_number][
                task_number
            ] = {}
            results_dict["reactions_included_per_task"][
                sample_number
            ][task_number] = {}
            results_dict["additional_settings_per_task"][
                sample_number
            ][task_number] = {}
    return results_dict

def load_expression_set(
    main_data_folder: str,
    model_version: str,
    model_type: str,
    expression_type: str,
) -> Tuple:
    pass

def load_main_model(
    main_data_folder: str,
    model_version: str,
    model_type: str,
) -> Model:
    pass

def check_if_method_inputs_are_correct(methods_dict: dict) -> None:
    pass

def run_single_method_optimized(
        name_of_run: str,
        test_array: list,
        sample_array: list,
        methods_dict: dict,
        main_data_folder: str,
        main_output_folder_original: str,
) -> None:

    results_dict = create_empty_dictionary_for_saving_results(
        name_of_run, test_array, sample_array, methods_dict
    )

    main_output_folder = path.join(main_output_folder_original, name_of_run)
    makedirs(main_output_folder, exist_ok=True)

    global_options_dict = methods_dict.get("__global_options", {})
    expression_type = global_options_dict.get("expression_type", "local_threshold")
    model_version = global_options_dict.get("model_version", "consensus_2024")
    model_type = (methods_dict.keys()[1] if len(methods_dict) > 1 else methods_dict.keys()[0])[0]
    settings_dict = (methods_dict.keys()[1] if len(methods_dict) > 1 else methods_dict.keys()[0])[2]
    expression_set = load_expression_set(main_data_folder, model_version, model_type, expression_type)
    cobra_model = load_main_model(main_data_folder, model_version, model_type)
    check_if_method_inputs_are_correct(methods_dict)



def run_test_suite_different_model_formulations(
    name_of_run: str,
    test_array: list,
    sample_array: list,
    methods_dict: dict,
    main_data_folder: str,
    main_output_folder_original: str,
) -> None:
    ### This is NOT a legacy function, and thus calculates the specific expression matrix from the rev2irrev file
    ### Legacy function works with models created in the matlab code, use the non legacy versions for new models
    ### Main differences are the manner by which expression json files are used. Legacy functions load the matfiles
    ### and extract the expression values from there. New functions use base expression json files directly
    ### to save space, using base expression and creating task specific expression is done in the loading of the models
    ### specifically by calculating the differences of amount of reactions and the rev2irrev files

    # methods_dict contains version_names as keys and functions as values
    # test_array contains the test numbers to run
    # users create functions which can be added to the dict and ran
    results_dict = {}
    results_dict["run_name"] = name_of_run
    results_dict["test_array"] = test_array
    results_dict["sample_array"] = sample_array
    results_dict["run_start_time_date"] = str(
        datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    main_output_folder = path.join(main_output_folder_original, name_of_run)
    makedirs(main_output_folder, exist_ok=True)
    global_options_dict = methods_dict.get("__global_options", {})
    expression_type = global_options_dict.get("expression_type", "local_threshold")
    if expression_type == "global_no_threshold":
        for version_name, tup in methods_dict.items():
            if version_name.startswith("__global_options"):
                continue
            if "expressionless_reaction_value" in tup[2]:
                old_expressionless_reaction_value = tup[2]["expressionless_reaction_value"]
                tup[2]["expressionless_reaction_value"] = old_expressionless_reaction_value * 3.827 #TODO floating constant is to increase
                # the expressionless_reaction_value when the values are not log transformed, due to the steps that the value went through during and after log
                # transform, which now is not the case, this approximation is going to be somewhat limited. A indepth investigation of the distribuion of these values
                # will be important!

    full_reversible_expression = cM.load_full_reversible_expression(main_data_folder, expression_type)

    for version_name, tup in methods_dict.items():
        if version_name.startswith("__global_options"):
            continue

        results_dict[version_name] = {
            "model_creation_time_per_task": {},
            "solving_time_per_task": {},
            "work_units_per_task": {},
            "reactions_included_per_task": {},
            "additional_settings_per_task": {},
        }
        for sample_number in sample_array:
            results_dict[version_name]["model_creation_time_per_task"][
                sample_number
            ] = {}
            results_dict[version_name]["solving_time_per_task"][sample_number] = {}
            results_dict[version_name]["work_units_per_task"][sample_number] = {}
            results_dict[version_name]["reactions_included_per_task"][
                sample_number
            ] = {}
            results_dict[version_name]["additional_settings_per_task"][
                sample_number
            ] = {}

    for task_number in test_array:
        pass_task, amount_of_reactions_added = (
            cM.check_if_pass_task_and_get_amount_of_reactions_added(
                task_number, main_data_folder
            )
        )
        if not pass_task:
            continue
        cobra_model_irreversible, cobra_model_reversible, rev2irrev = (
            cM.load_irrev_and_rev_model_from_task_specific_mat_files(
                task_number, main_data_folder
            )
        )


        expression_irreversible_all_samples, expression_reversible_all_samples = (
            cM.create_task_specific_expression_from_full_reversible(
                full_reversible_expression, rev2irrev, amount_of_reactions_added
            )
        )

        cobra_model_irreversible = sR.adjust_backward_reaction_names(
            cobra_model_irreversible
        )
        print(f"Running task {task_number}\n\n")

        for sample_number in sample_array:
            try:
                intermediate_json_name = (
                    f"Task_{task_number:03}_Sample{sample_number:03}_{expression_type}.json"
                )
                with open(
                    path.join(main_output_folder, intermediate_json_name), "w"
                ) as f:
                    dump(results_dict, f)
            except Exception as e:
                print(e)
                print("Could not save intermediate json, already exists")
            print(f"Running task {task_number} with sample {sample_number}\n\n")
            if type(expression_reversible_all_samples[0]) == list:
                if (
                    len(expression_reversible_all_samples[0]) >= sample_number
                    and len(expression_irreversible_all_samples[0]) >= sample_number
                ):
                    expression_reversible = [
                        expression[sample_number - 1]
                        for expression in expression_reversible_all_samples
                    ]
                    expression_irreversible = [
                        expression[sample_number - 1]
                        for expression in expression_irreversible_all_samples
                    ]
                else:
                    expression_reversible = deepcopy(expression_reversible_all_samples)
                    expression_irreversible = deepcopy(
                        expression_irreversible_all_samples
                    )
            elif type(expression_reversible_all_samples[0]) == np.ndarray:
                if (
                    np.size(expression_reversible_all_samples[0]) >= sample_number
                    and np.size(expression_irreversible_all_samples[0]) >= sample_number
                ):
                    expression_reversible = [
                        expression[sample_number - 1]
                        for expression in expression_reversible_all_samples
                    ]
                    expression_irreversible = [
                        expression[sample_number - 1]
                        for expression in expression_irreversible_all_samples
                    ]
                else:
                    expression_reversible = deepcopy(expression_reversible_all_samples)
                    expression_irreversible = deepcopy(
                        expression_irreversible_all_samples
                    )

            for version_name, tup in methods_dict.items():
                if version_name.startswith("__global_options"):
                    continue
                print(
                    f"Running task {expression_type}_{task_number} for {version_name} with sample {sample_number}"
                )
                solver = SolverFactory("gurobi")
                # results_dict[version_name]["model_creation_time_per_task"][sample_number][task_number] = {}
                # results_dict[version_name]["solving_time_per_task"][sample_number][task_number] = {}
                # results_dict[version_name]["reactions_included_per_task"][sample_number][task_number] = {}
                results_dict[version_name]["additional_settings_per_task"][
                    sample_number
                ][task_number] = {}
                model = None
                chosen_expression = None
                chosen_expression_old = None
                chosen_model = None

                type_of_cobra_model = tup[0]
                function_to_run = tup[1]
                settings_kwargs_dict = tup[2]
                if settings_kwargs_dict.get("create_log_file", False):
                    log_file_iterator = 1
                    log_file_name = f"{main_output_folder}\\{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_gurobi.log"
                    while True:
                        if path.exists(log_file_name):
                            log_file_name = f"{main_output_folder}\\{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_gurobi_Try{log_file_iterator}.log"
                            log_file_iterator += 1
                        else:
                            break
                    solver.options["LogFile"] = log_file_name

                if "solver_settings" in settings_kwargs_dict:
                    for key, value in settings_kwargs_dict["solver_settings"].items():
                        solver.options[key] = value
                settings_kwargs_dict["solver"] = solver
                if type_of_cobra_model == "irreversible":
                    chosen_model = cobra_model_irreversible
                    chosen_expression_old = expression_irreversible
                    time_start = time()
                    model, chosen_expression = function_to_run(
                        chosen_model,
                        chosen_expression_old,
                        rev2irrev,
                        expression_reversible,
                        **settings_kwargs_dict,
                    )
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = (time() - time_start)
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = chosen_expression[
                        expression_irreversible.index(-1)
                    ]
                elif type_of_cobra_model == "reversible":
                    chosen_model = cobra_model_reversible
                    chosen_expression_old = expression_reversible
                    time_start = time()
                    model, chosen_expression = function_to_run(
                        cobra_model_reversible,
                        expression_reversible,
                        **settings_kwargs_dict,
                    )
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = (time() - time_start)
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = chosen_expression[
                        expression_reversible.index(-1)
                    ]
                else:
                    warn("Invalid type of cobra model")
                    results_dict[version_name]["model_creation_time_per_task"][
                        sample_number
                    ][task_number] = -999
                    results_dict[version_name]["additional_settings_per_task"][
                        sample_number
                    ][task_number]["expressionless_reaction_value"] = -999
                if model is None:
                    warn(f"Model for task {task_number} was not created")
                    results_dict[version_name]["solving_time_per_task"][sample_number][
                        task_number
                    ] = -999
                    results_dict[version_name]["work_units_per_task"][sample_number][
                        task_number
                    ] = -999
                    continue
                else:
                    instance = model.create_instance()
                    result = solver.solve(instance, tee=True)
                    results_dict[version_name]["solving_time_per_task"][sample_number][
                        task_number
                    ] = result.Solver.Time
                    solver_output = solver._log
                    pattern = r'\((\d+\.\d+) work units\)'
                    matches = findall(pattern, solver_output)
                    second_work_units = 999
                    if matches is not None and len(matches) >= 2:
                        second_work_units = matches[1] # TODO this is not working as intended as in bottleneck we find 999 sometimes
                    results_dict[version_name]["work_units_per_task"][sample_number][
                        task_number
                    ] = second_work_units
                    if "epsilon" in settings_kwargs_dict:
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["epsilon"] = settings_kwargs_dict["epsilon"]

                    if result.Solver.Termination_condition == "optimal":
                        (
                            current_reactions,
                            (
                                amount_of_reactions,
                                amount_of_expressionless_reactions,
                                amount_of_expression_carrying_reactions,
                            ),
                            new_model,
                        ) = sR.get_active_reactions_from_solved_model(
                            instance,
                            chosen_expression_old,
                            chosen_expression,
                            chosen_model,
                            **settings_kwargs_dict,
                        )
                        file_name_ = f"model_json_{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}.json"
                        location_to_save_filtered = path.join(
                            main_output_folder, file_name_
                        )
                        file_iterator = 1
                        while True:
                            if path.exists(location_to_save_filtered):
                                file_name_ = f"model_json_{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_Try{file_iterator}.json"
                                location_to_save_filtered = path.join(
                                    main_output_folder, file_name_
                                )
                                file_iterator += 1
                            else:
                                break

                        save_json_model(new_model, location_to_save_filtered)
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
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["expression_type"] = expression_type
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["identifier"] = f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}"
                        epsilon = settings_kwargs_dict.get("epsilon", 1e-3)
                        sR.save_solver_results(
                            instance,
                            result,
                            solver,
                            chosen_model,
                            f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}",
                            main_output_folder,
                            chosen_expression,
                            chosen_expression_old,
                            epsilon,
                        )
                    else:
                        warn(f"Task {task_number} for {version_name} did not converge")
                        results_dict[version_name]["reactions_included_per_task"][
                            sample_number
                        ][task_number] = {}
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["amount_of_reactions"] = -999
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["amount_of_expressionless_reactions"] = -999
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["amount_of_expression_carrying_reactions"] = -999
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["expression_type"] = expression_type
                        results_dict[version_name]["additional_settings_per_task"][
                            sample_number
                        ][task_number]["identifier"] = f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}"

    results_dict["run_end_time_date"] = str(
        datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    print(results_dict)
    print("")
    location_to_save = path.join(main_output_folder, f"{name_of_run}_results_{expression_type}")
    # to json
    counter = 1
    filename_without_extension = path.splitext(location_to_save)[0]
    saved_file = False
    try:
        filename = f"{filename_without_extension}.json"
        if not path.exists(filename):
            with open(filename, "w") as f:
                dump(results_dict, f)
            saved_file = True
    except Exception as e:
        print(e)
    while not saved_file:
        try:
            filename_with_counter = f"{filename_without_extension}_Try{counter}.json"
            if not path.exists(filename_with_counter):
                with open(filename_with_counter, "w") as f:
                    dump(results_dict, f)
                break
            else:
                counter += 1
        except Exception as e:
            print(e)
            break

    # main_folder = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\Base_files"
    # model = cobra.io.load_matlab_model(main_folder + "\\" + "model_loaded.mat")
    # expression, significance, task_structure = load_json_data(main_folder)
    # model, expression = prepare_data_for_fastcc(model, expression, significance)
    # model = close_exchange_reactions(model)
    # irrev_model, rev2irrev = convert_to_irreversible(model)
    # amount_of_reactions_difference = len(irrev_model.reactions) - len(model.reactions)
    # for idx in range(amount_of_reactions_difference):
    #     expression.append(["Fill"] * len(expression[0]))
    # for idx in range(len(rev2irrev)):
    #     if len(rev2irrev[idx]) > 1:
    #         expression[rev2irrev[idx][1]] = expression[rev2irrev[idx][0]]
    #
    # cobra.io.save_matlab_model(irrev_model, main_folder + "\\" + "irrev_model.mat")
    # json.dump(expression, open(main_folder + "\\" + "expression_irreversible.json", 'w'))
    # json.dump(rev2irrev, open(main_folder + "\\" + "rev2irrev.json", 'w'))
    #

    # load full reversible expression
    # select task,
    # check if pass task and load amount of reactions added
    # load task_specific_reversible model
    # load task_specific_irreversible model
    # load task_specific_rev2irrev
    # create task_specific reversible expression: add amount of reactions added to the end of the expression
    # take task specific rev2irrev and check difference in length with reversible (if rev2irrev contains more reactions,
    #       then some added temporary reactions are also reversible)
    # create task_specific irreversible expression

def basic_loading_and_preparing_models(
        task_number: int,
        sample_number: int,
        methods_dict: dict,  # should not contain any of the gurobi options that will be checked
        main_data_folder: str,
) -> Tuple:
    global_options_dict = methods_dict.get("__global_options", {})
    expression_type = global_options_dict.get("expression_type", "local_threshold")
    if expression_type == "global_no_threshold":
        for version_name, tup in methods_dict.items():
            if version_name.startswith("__global_options"):
                continue
            if "expressionless_reaction_value" in tup[2]:
                old_expressionless_reaction_value = tup[2]["expressionless_reaction_value"]
                tup[2]["expressionless_reaction_value"] = old_expressionless_reaction_value * 3.827
    full_reversible_expression = cM.load_full_reversible_expression(main_data_folder, expression_type)
    pass_task, amount_of_reactions_added = (
        cM.check_if_pass_task_and_get_amount_of_reactions_added(
            task_number, main_data_folder
        )
    )
    if not pass_task:
        print("This task does pass the pass_task check, aborted tests")
        return (False, False, False, False, False, False)
    cobra_model_irreversible, cobra_model_reversible, rev2irrev = (
        cM.load_irrev_and_rev_model_from_task_specific_mat_files(
            task_number, main_data_folder
        )
    )
    expression_irreversible_all_samples, expression_reversible_all_samples = (
        cM.create_task_specific_expression_from_full_reversible(
            full_reversible_expression, rev2irrev, amount_of_reactions_added
        )
    )
    cobra_model_irreversible = sR.adjust_backward_reaction_names(
        cobra_model_irreversible
    )
    if type(expression_reversible_all_samples[0]) == list:
        if (
                len(expression_reversible_all_samples[0]) >= sample_number
                and len(expression_irreversible_all_samples[0]) >= sample_number
        ):
            expression_reversible = [
                expression[sample_number - 1]
                for expression in expression_reversible_all_samples
            ]
            expression_irreversible = [
                expression[sample_number - 1]
                for expression in expression_irreversible_all_samples
            ]
        else:
            expression_reversible = deepcopy(expression_reversible_all_samples)
            expression_irreversible = deepcopy(
                expression_irreversible_all_samples
            )
    elif type(expression_reversible_all_samples[0]) == np.ndarray:
        if (
                np.size(expression_reversible_all_samples[0]) >= sample_number
                and np.size(expression_irreversible_all_samples[0]) >= sample_number
        ):
            expression_reversible = [
                expression[sample_number - 1]
                for expression in expression_reversible_all_samples
            ]
            expression_irreversible = [
                expression[sample_number - 1]
                for expression in expression_irreversible_all_samples
            ]
        else:
            expression_reversible = deepcopy(expression_reversible_all_samples)
            expression_irreversible = deepcopy(
                expression_irreversible_all_samples
            )
    return(
        cobra_model_irreversible,
        cobra_model_reversible,
        rev2irrev,
        expression_reversible,
        expression_irreversible,
        expression_type,
    )

def perform_single_run(
    model,
    solver: SolverFactory,
    settings_kwargs_dict: dict,
    time_limit: int = 0,
    percentage_longer_to_stop: float = 0.1,
    best_time_with_version: float = 999999999,
    log_file_name: str = ""
):
    # perform a single run with regular settings to get the best time
    if settings_kwargs_dict.get("create_log_file", False):
        log_file_iterator = 1
        log_file_name_to_use = f'{log_file_name}.log'
        while True:
            if path.exists(log_file_name):
                log_file_name_to_use = f"{log_file_name}_Try{log_file_iterator}.log"
                log_file_iterator += 1
            else:
                break
        solver.options["LogFile"] = log_file_name_to_use
    if time_limit > 0:
        solver.options["TimeLimit"] = time_limit
    else:
        solver.options.pop("TimeLimit", None)
    if "solver_settings" in settings_kwargs_dict:
        for key, value in settings_kwargs_dict["solver_settings"].items():
            solver.options[key] = value

    instance = model.create_instance()
    try:
        result = solver.solve(instance, tee=True)
        if result.Solver.Termination_condition == "optimal":
            print("Optimal solution found")
            print(result.Solver.Time, best_time_with_version)
            if result.Solver.Time < best_time_with_version:
                print("New best time found")
                best_time_with_version = result.Solver.Time * (1 + percentage_longer_to_stop)
        current_run_time = result.Solver.Time
    except:
        current_run_time = 999999999


    return (best_time_with_version, current_run_time, solver)

def run_with_several_seeds(
    run_for_comparison: bool,
    amount_of_seeds_to_try: int, # default
    version_name: str,
    task_number: int,
    sample_number: int,
    model,
    expression_type: str,
    solver: SolverFactory,
    settings_kwargs_dict: dict,
    main_output_folder: str,
    time_limit: int = 0,
    percentage_longer_to_stop: float = 0.1,
):
    seed_range = range(1, amount_of_seeds_to_try + 1)
    runtimes = []
    for seed_value in seed_range:
        solver.options["Seed"] = seed_value
        if run_for_comparison:
            log_file_name = (f"{main_output_folder}\\Comparison_{version_name}_"
                             f"Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_gurobi_seed{seed_value}")
        else:
            log_file_name = (f"{main_output_folder}\\{version_name}_"
                             f"Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_gurobi_seed{seed_value}")
            # we should use the worst of a set of seeds as the best time, in order for other runs with different seeds to not be aborted prematurely
        best_time, current_run_time, solver = perform_single_run(
            model,
            solver,
            settings_kwargs_dict,
            time_limit,
            percentage_longer_to_stop,
            log_file_name=log_file_name,
        )
        runtimes.append(current_run_time)
    return runtimes

def calculate_trimmed_mean_and_std(
    runtimes: List[float],
    percentage_to_trim: float = 0.1,
) -> Tuple[float, float]:
    runtimes.sort()
    amount_to_trim = int(len(runtimes) * percentage_to_trim)
    trimmed_runtimes = runtimes[amount_to_trim:-amount_to_trim]
    mean = np.mean(trimmed_runtimes)
    std = np.std(trimmed_runtimes)
    return (mean, std)

def calculate_descriptive_statistics_and_save_as_txt(
        main_output_folder: str,
        version_name: str,
        task_number: int,
        sample_number: int,
        expression_type: str,
        runtimes: List[float],
        percentage_to_trim: float = 0.1,
):
    trimmed_mean, trimmed_sd = calculate_trimmed_mean_and_std(runtimes, percentage_to_trim)
    mean = np.mean(runtimes)
    std = np.std(runtimes)
    sem = np.std(runtimes) / np.sqrt(len(runtimes))
    min_time = min(runtimes)
    max_time = max(runtimes)
    median_time = np.median(runtimes)
    list_with_stats = [mean, std, sem, min_time, max_time, median_time, trimmed_mean, trimmed_sd]
    [print(f"{name}: {value}") for name, value in zip(["mean", "std", "sem", "min", "max", "median", "trimmed_mean", "trimmed_sd"], list_with_stats)]
    df = pd.DataFrame(
        {
            "mean": [mean],
            "std": [std],
            "sem": [sem],
            "min": [min_time],
            "max": [max_time],
            "median": [median_time],
            "trimmed_mean": [trimmed_mean],
            "trimmed_sd": [trimmed_sd],
        }
    )
    df.to_csv(
        path.join(
            main_output_folder,
            f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_descriptive_statistics.csv",
        ),
        sep="\t",
    )
    with open(
        path.join(
            main_output_folder,
            f"{version_name}_Task_{task_number:03}_Sample{sample_number:03}_{expression_type}_descriptive_statistics.txt",
        ),
        "a",
    ) as f:
        f.write(f"\n"
                f"runtimes: {runtimes}\n"
                f"mean: {mean}\n"
                f"std: {std}\n"
                f"sem: {sem}\n"
                f"min: {min_time}\n"
                f"max: {max_time}\n"
                f"median: {median_time}\n"
                f"trimmed_mean: {trimmed_mean}\n"
                f"trimmed_sd: {trimmed_sd}\n ")


def test_different_options_for_best_solving_or_best_incumbent(
    name_of_run: str,
    task_number: int,
    sample_number: int,
    methods_dict: dict, # should not contain any of the gurobi options that will be checked
    to_iterate_and_check_options: dict, # should contain the gurobi options that will be checked {option: [list of values to check]}
    main_data_folder: str,
    main_output_folder_original: str,
    time_limit: int, # in the case of finding the best/fastest incumbent, this should be set to some non-zero value
    percentage_longer_to_stop: float = 0.1,
) -> None:
    main_output_folder = path.join(main_output_folder_original, name_of_run)
    makedirs(main_output_folder, exist_ok=True)
    (cobra_model_irreversible, cobra_model_reversible, rev2irrev, expression_reversible,
     expression_irreversible, expression_type) = basic_loading_and_preparing_models(
        task_number, sample_number, methods_dict, main_data_folder
    )
    if cobra_model_reversible is False:
        print("This task does pass the pass_task check, aborted tests")
        return
    permutation_of_all_options = []
    for option, values in to_iterate_and_check_options.items():
        if len(permutation_of_all_options) == 0:
            for value in values:
                permutation_of_all_options.append({option: value})
        else:
            new_permutation = []
            for value in values:
                for old_permutation in permutation_of_all_options:
                    new_permutation.append({**old_permutation, option: value})
            permutation_of_all_options = new_permutation

    for version_name, tup in methods_dict.items():
        model = None
        if version_name.startswith("__global_options"):
            continue
        type_of_cobra_model = tup[0]
        function_to_run = tup[1]
        settings_kwargs_dict = tup[2]
        if type_of_cobra_model == "irreversible":
            chosen_model = cobra_model_irreversible
            chosen_expression_old = expression_irreversible
            model, chosen_expression = function_to_run(
                chosen_model,
                chosen_expression_old,
                rev2irrev,
                expression_reversible,
                **settings_kwargs_dict,
            )
        solver = SolverFactory("gurobi")
        runtimes_for_comparison = run_with_several_seeds(
            True,
            10,
            version_name,
            task_number,
            sample_number,
            model,
            expression_type,
            solver,
            settings_kwargs_dict,
            main_output_folder,
            time_limit,
        )
        solver.options = {}
        worst_seed_best_time_with_version = max(runtimes_for_comparison) * (1 + percentage_longer_to_stop)
        calculate_descriptive_statistics_and_save_as_txt(
            main_output_folder,
            version_name,
            task_number,
            sample_number,
            expression_type,
            runtimes_for_comparison,
        )
        for permutation in permutation_of_all_options:
            if worst_seed_best_time_with_version > 0:
                solver.options["TimeLimit"] = worst_seed_best_time_with_version
            else:
                solver.options.pop("TimeLimit", None)
            if "solver_settings" in settings_kwargs_dict:
                for key, value in settings_kwargs_dict["solver_settings"].items():
                    solver.options[key] = value
            for key, value in permutation.items():
                solver.options[key] = value
            runtimes = run_with_several_seeds(
                False,
                5,
                version_name,
                task_number,
                sample_number,
                model,
                expression_type,
                solver,
                settings_kwargs_dict,
                main_output_folder,
                worst_seed_best_time_with_version,
            )
            if (max(runtimes) * (1 + percentage_longer_to_stop)) < worst_seed_best_time_with_version:
                worst_seed_best_time_with_version = max(runtimes) * (1 + percentage_longer_to_stop)
                print(f"New best time found: {worst_seed_best_time_with_version} for {version_name} with {permutation}")

def load_timelines_and_create_distilled_df(
        data_folder: str,
        log_file_to_compare_to_location: str,
        time_limit: int = None,
        task_number: int = None,
        sample_number: int = None,
        expression_type: str = None,
        name_of_run: str = None,

):
    if name_of_run is not None:
        name_of_run_data_folder = name_of_run
    else:
        name_of_run_from_data_folder = path.split(data_folder)[1]
    # name_of_run = f"linear_best_incumbent_Time{time_limit}_T{task_number:03}_S{sample_number:03}_E{expression_type}"
    names_dict = {
        "time_limit": time_limit,
        "task_number": task_number,
        "sample_number": sample_number,
        "expression_type": expression_type,
    }
    searches_list = [
        r"Time(\d+)_",
        r"T(\d+)_",
        r"S(\d+)_",
        r"E(\w+)",
    ]
    for search in searches_list:
        try:
            names_dict[search] = int(findall(search, name_of_run_from_data_folder)[0])
        except:
            pass
    time_limit = names_dict["time_limit"]
    task_number = names_dict["task_number"]
    sample_number = names_dict["sample_number"]
    expression_type = names_dict["expression_type"]

    summary_to_compare_to, timelines_to_compare_to = glt_gurb.get_dataframe([log_file_to_compare_to_location], timelines=True)
    all_other_logs, all_other_timelines =  glt_gurb.get_dataframe([f"{'*'}.log"], timelines=True)

    df_to_create_columns = [
        "Best incumbent","Time best incumbent reached",	"Prevous best with normal settings at time", "Best time to reach in original"
    "is_better_than_previous_at_time",	"Optimal value",	"optimal Found at time",
    "Task",	"Sample"	"Expression_type",	"Gap",	"MIPFocus",	"Heuristics",	"Cuts",	"Presolve",	"PreSparsify",	"NodeMethod",	"Method" ,	"Aggregate",	"PreDual"
    ]
    df = pd.DataFrame(columns=df_to_create_columns)


def check_SV_constraints_and_model_feasilibity(
        cobra_model: Model,
):
    solution = cobra_model.optimize()
    if solution.status == "optimal":
        print("Model is feasible")
        return True
    else:
        print("Model is infeasible")
        return False



def check_all_json_models_in_dir_for_feasibility(
        dir_with_models: str,
) -> Dict[str, bool]:
    is_feasible = {}
    for file in listdir(dir_with_models):
        if file.endswith(".json") and file.startswith("model_json_"):
            try:
                ends = file.split(".json")[0]
                filename = filename.split("model_json_")[1]
                model = load_json_model(path.join(dir_with_models, file))
                is_feasible[filename] = check_SV_constraints_and_model_feasilibity(model)
                print(f"Model {file} checked")
            except:
                print(f"Model {file} could not be loaded")
                continue
    return is_feasible

