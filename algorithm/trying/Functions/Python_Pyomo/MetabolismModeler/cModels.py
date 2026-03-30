from typing import List, Tuple
# noinspection PyUnresolvedReferences
from pyomo.environ import (
    ConcreteModel,
    RangeSet,
    Var,
    Constraint,
    Objective,
    minimize,
    Binary,
    Reals,
    NonNegativeReals,
    NonNegativeIntegers,
    Integers,
    RealSet,
    BooleanSet,
)
import cobra
import numpy as np
from os.path import join
from os import getcwd, walk
from re import search
from glob import glob
from scipy.io import loadmat
from logging import getLogger, ERROR
from warnings import catch_warnings, filterwarnings
from json import load as json_load, dump as json_dump
from copy import deepcopy

__all__ = [ # create models
    "create_basic_pyomo_model_SV_constraints",
    "__LEGACY__load_irreversible_model_and_expression_from_mat_files",
    "load_reversible_model_and_expression_from_mat_files",
    "load_full_reversible_expression",
    "check_if_pass_task_and_get_amount_of_reactions_added",
]


def load_full_reversible_expression(main_folder_ori: str,
                                    expression_type: str = "local_threshold",
                                    ) -> List[List[float]]:
    main_folder = main_folder_ori
    if expression_type == "local_threshold":
        file_to_get = "full_expression_reversible.json"
    elif expression_type == "global_percentile_90":
        file_to_get = "full_expression_reversible_percentile_90.json"
    elif expression_type == "global_no_threshold":
        file_to_get = "full_expression_reversible_no_threshold.json"
    else:
        raise ValueError("Expression type not recognized")
    full_reversible_expression = None
    search_pattern = join(main_folder, file_to_get)
    files = glob(search_pattern, recursive=True)
    if files:
        for file in files:
            with open(file, "r") as f:
                full_reversible_expression = json_load(f)
    if full_reversible_expression is None:
        print("Full reversible expression was not found")
        raise ValueError("One of the necessary files was not found")
    else:
        print("Full reversible expression loaded successfully")
    return full_reversible_expression


def check_if_pass_task_and_get_amount_of_reactions_added(
    task_number: int, main_folder: str
) -> (bool, int):
    file_to_get = f"pass_task_and_reactions_added_task{task_number:03}.json"
    amount_of_reactions_added = None
    pass_task = None
    task_directory_name = f"Task_{task_number:03}"
    search_pattern = join(main_folder, task_directory_name, file_to_get)
    files = glob(search_pattern, recursive=True)
    if files:
        for file in files:
            with open(file, "r") as f:
                data = json_load(f)
                pass_task = data[0]
                amount_of_reactions_added = data[1]
    if amount_of_reactions_added is None or pass_task is None:
        if amount_of_reactions_added is None:
            print("Amount of reactions added was not found")
        if pass_task is None:
            print("Pass task was not found")
        raise ValueError("One of the necessary files was not found")
    else:
        print("Amount of reactions added and pass task loaded successfully")
    return pass_task, amount_of_reactions_added


def load_irrev_and_rev_model_from_task_specific_mat_files(
    task_number: int, main_folder: str
) -> (cobra.Model, cobra.Model, np.ndarray):
    filenames_to_get = [
        f"irrev_model_task{task_number:03}.mat",
        f"task_specific_model_task{task_number:03}.mat",
        f"rev2irrev_task{task_number:03}.json",
    ]
    irrev_model = None
    rev_model = None
    rev2irrev = None
    task_directory_name = f"Task_{task_number:03}"
    for file_to_get in filenames_to_get:
        search_pattern = join(main_folder, task_directory_name, file_to_get)
        files = glob(search_pattern, recursive=True)
        if files:
            for file in files:
                if file_to_get == f"irrev_model_task{task_number:03}.mat":
                    with catch_warnings():
                        original_logging_level = getLogger().getEffectiveLevel()
                        getLogger().setLevel(ERROR)
                        filterwarnings("ignore", "No defined compartments in model")
                        print(f"Loading {file}")
                        irrev_model = cobra.io.load_matlab_model(file)
                        getLogger().setLevel(original_logging_level)
                        print(f"{file} is loaded successfully")

                elif file_to_get == f"task_specific_model_task{task_number:03}.mat":
                    with catch_warnings():
                        original_logging_level = getLogger().getEffectiveLevel()
                        getLogger().setLevel(ERROR)
                        filterwarnings("ignore", "No defined compartments in model")
                        print(f"Loading {file}")
                        rev_model = cobra.io.load_matlab_model(file)
                        getLogger().setLevel(original_logging_level)
                        print(f"{file} is loaded successfully")
                elif file_to_get == f"rev2irrev_task{task_number:03}.json":
                    with open(file, "r") as f:
                        rev2irrev = json_load(f)

    if irrev_model is None or rev_model is None:
        if irrev_model is None:
            print("Irrev model was not found")
        if rev_model is None:
            print("Rev model was not found")
        raise ValueError("One of the necessary files was not found")
    else:
        print("Irrev and Rev model loaded successfully")
    return irrev_model, rev_model, rev2irrev


def create_irrev_expression_from_reversible(
    reversible_expression: List[List[float]],
    rev2irrev: List[List[int]],
    amount_of_reactions_added: int,
) -> List[List[float]]:
    irrev_expression = deepcopy(reversible_expression)
    length_irrev = len(rev2irrev)
    amount_of_erronious_added = 0
    part_to_check_for_erronious = rev2irrev[-amount_of_reactions_added:]
    for part in part_to_check_for_erronious:
        if len(part) > 1:
            amount_of_erronious_added += 1
    part_that_is_wrong = rev2irrev[:-amount_of_reactions_added]
    for idx in range(len(part_that_is_wrong)):
        if len(rev2irrev[idx]) > 1:
            rev2irrev[idx][1] += amount_of_erronious_added
    max_length = max(
        (rev2irrev[idx][1]) for idx in range(length_irrev) if len(rev2irrev[idx]) > 1
    )
    diff = max_length - length_irrev
    for idx in range(diff):
        irrev_expression.append([-1] * len(irrev_expression[0]))

    for idx in range(length_irrev):
        if len(rev2irrev[idx]) > 1:
            irrev_expression[rev2irrev[idx][1] - 1] = reversible_expression[idx]

    return irrev_expression


def create_task_specific_expression_from_full_reversible(
    full_reversible_expression: List[List[float]],
    rev2irrev: List[List[int]],
    amount_of_reactions_added: int,
) -> (np.ndarray, np.ndarray):

    reversible_expression = deepcopy(full_reversible_expression)
    for idx in range(amount_of_reactions_added):
        reversible_expression.append([-1] * len(reversible_expression[0]))

    irrev_expression = create_irrev_expression_from_reversible(
        reversible_expression, rev2irrev, amount_of_reactions_added
    )

    return irrev_expression, reversible_expression


def create_basic_pyomo_model_SV_constraints(
    stoichiometric_matrix: np.ndarray,
    cobra_model: cobra.Model,
    using_cobra: bool = True,
) -> (ConcreteModel, dict):
    model = ConcreteModel()
    num_metabolites = len(stoichiometric_matrix)  # Example number of metabolites
    num_reactions = len(stoichiometric_matrix[0])  # Example number of reactions

    model.Metabolites = RangeSet(1, num_metabolites)
    model.Rxns = RangeSet(1, num_reactions)
    stoich_dict = {
        met: {
            rxn: stoichiometric_matrix[met, rxn]
            for rxn in range(num_reactions)
            if stoichiometric_matrix[met, rxn] != 0
        }
        for met in range(num_metabolites)
    }
    model.V_flux = Var(model.Rxns, domain=Reals)
    # noinspection PyTypeChecker, PyUnresolvedReferences
    for rxn in model.Rxns:
        # noinspection PyTypeChecker, PyUnresolvedReferences
        model.V_flux[rxn].setlb(cobra_model.reactions[rxn - 1].lower_bound)
        # noinspection PyTypeChecker, PyUnresolvedReferences
        model.V_flux[rxn].setub(cobra_model.reactions[rxn - 1].upper_bound)

    def sv_zero_rule(model, met):
        if stoich_dict[met - 1] != {}:
            return (
                sum(
                    stoich_dict[met - 1][rxn] * model.V_flux[rxn + 1]
                    for rxn in stoich_dict[met - 1]
                )
                == 0
            )
        else:
            return Constraint.Skip

    model.SV_constraint = Constraint(model.Metabolites, rule=sv_zero_rule)

    return model, stoich_dict


def __LEGACY__load_irreversible_model_and_expression_from_mat_files(
    task_number: int, main_folder: str = "", directory_path: str = ""
) -> (cobra.Model, list, list):
    '''This code is legacy and should not be used. It is kept for reference purposes only.'''
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
            folder_to_search = getcwd()
        else:
            folder_to_search = main_folder
        all_dirs = [dirpath for dirpath, dirnames, files in walk(folder_to_search)]
        pattern = rf"Task0*{task_number}S"
        for dir_ in all_dirs:
            if search(pattern, dir_):
                directory_path = dir_
                break
    print("Loading files from: ", directory_path)
    for file_to_get in filenames_to_get:
        search_pattern = join(directory_path, "**", file_to_get)
        files = glob(search_pattern, recursive=True)
        if files:
            for file in files:
                if file_to_get == "newIrrevModelKapprox.mat":
                    with catch_warnings():
                        original_logging_level = getLogger().getEffectiveLevel()
                        getLogger().setLevel(ERROR)
                        filterwarnings("ignore", "No defined compartments in model")
                        model = cobra.io.load_matlab_model(file)
                        getLogger().setLevel(original_logging_level)

                elif file_to_get == "expValueMeanAdjustedTaskIrrevKapprox.mat":
                    expression_cells = loadmat(file)
                    expression = [
                        value[0]
                        for value in expression_cells[
                            "expValueMeanAdjustedTaskIrrevKapprox"
                        ]
                    ]

                elif file_to_get == "rev2irrevIrrevKapprox.mat":
                    rev2irrev_cells = loadmat(file)
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
    main_folder_1 = "E:\\Git\\Metabolic_Task_Score\\Data\\Base MAT Files\\Consensus Model 1.17\\No Thresh"
    # work_model, expression, rev2irrev =
    # # print(model)
    # # print(expression)
    # # print(rev2irrev)


def load_reversible_model_and_expression_from_mat_files(
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
            folder_to_search = getcwd()
        else:
            folder_to_search = main_folder
        all_dirs = [dirpath for dirpath, dirnames, files in walk(folder_to_search)]
        pattern = rf"Task0*{task_number}S"
        for dir_ in all_dirs:
            if search(pattern, dir_):
                directory_path = dir_
                break
    print("Loading files from: ", directory_path)
    for file_to_get in filenames_to_get:
        search_pattern = join(directory_path, "**", file_to_get)
        files = glob(search_pattern, recursive=True)
        if files:
            for file in files:
                if file_to_get == "newReversibleModelForKApproxmiation.mat":
                    with catch_warnings():
                        original_logging_level = getLogger().getEffectiveLevel()
                        getLogger().setLevel(ERROR)
                        filterwarnings("ignore", "No defined compartments in model")
                        model = cobra.io.load_matlab_model(file)
                        getLogger().setLevel(original_logging_level)

                elif (
                    file_to_get
                    == "expValueMeanAdjustedTaskReversibleKApproximation.mat"
                ):
                    expression_cells = loadmat(file)
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

if __name__ == "__main__":
    placeholder = 5
    # main_folder_1 = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing"
    # main_folder_1 = "E:\\Git\\Metabolic_Task_Score\\Data\\Base MAT Files\\Consensus Model 1.17\\No Thresh"
    # model, expression = load_reversible_model_and_expression_from_mat_files(
    #     placeholder, main_folder_1
    # )
    # print(expression)
    # json_file = "full_expression_reversible_no_threshold.json"
    # with open(json_file, "w") as f:
    #     json_dump(expression, f)
    # # print(model)
    # # print(expression)
    # # print(rev2irrev)