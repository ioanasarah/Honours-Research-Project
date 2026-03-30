from json import dump
import numpy as np
from pandas import DataFrame, concat
import pyomo.environ as pe
from pyomo.environ import Var
from cobra import Model
from cobra.io import load_matlab_model
from scipy.io import loadmat
from os import path, makedirs, getcwd
from .file_management import save_excel_or_CSV_if_file_already_exist
from re import findall

__all__ = [
    "save_files_as_json",
    "save_solver_results",
    "get_active_reactions_from_solved_model",
    "find_or_create_user_based_directory",
    "adjust_backward_reaction_names",
]


def adjust_backward_reaction_names(cobra_model: Model) -> Model:
    ### necessary since irreversbile model was created without reaction.names for the backward reactions, should not be necessary in the future
    for reaction in cobra_model.reactions:
        if reaction.id.endswith("_b") or reaction.id.endswith("_r"):
            forward_reaction_name_id = reaction.id[:-2]
            try:
                forward_reaction = cobra_model.reactions.get_by_id(
                    forward_reaction_name_id
                )
                reaction.name = forward_reaction.name
                forward_reaction.id = forward_reaction_name_id + "_f"
                if forward_reaction.lower_bound < 0:
                    reaction.bounds = (0, -forward_reaction.lower_bound)
                    forward_reaction.lower_bound = 0
            except:
                pass

    return cobra_model


def find_or_create_user_based_directory(
    output_folder: str, user_name: str, specific_folder=None, base_files: bool = False
) -> (str, str):
    # assume program is ran from a folder SRC inside a repository.
    # the SRC is inside Git\\Metabolic_Task_Score\\Data\\pyomo_python_files
    # the main_data_folder is then: Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing
    # user folder for output is then: Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\user_name\\output_folder
    temporary_path = getcwd()
    try:
        main_data_folder = path.abspath(
            path.join(getcwd(), "..", "..", "..", "Data", "pyomo_python_files")
        )
        user_folder = path.join(main_data_folder, user_name)
        makedirs(user_folder, exist_ok=True)
        output_folder = path.join(user_folder, output_folder)
        makedirs(output_folder, exist_ok=True)
        if base_files:
            main_data_folder = path.abspath(
                path.join(
                    getcwd(),
                    "..",
                    "..",
                    "..",
                    "Data",
                    "pyomo_python_files",
                    "Base_files",
                )
            )
        elif specific_folder is not None:
            main_data_folder = path.join(main_data_folder, specific_folder)
    except Exception as e:
        print(e)
        main_data_folder = temporary_path
        output_folder = temporary_path
    return main_data_folder, output_folder


def save_files_as_json(main_folder: str) -> None:
    workspace = loadmat(main_folder + "\\" + "matlab_workspace.mat")
    expression = workspace["expressionRxns"]
    significance = workspace["significance"]
    task_structure = workspace["taskStructure"]
    # save expression and significance as json
    # save task_structure as json
    list_of_expressions = []
    for expr in expression:
        list_of_expressions.append(expr.tolist())

    with open(main_folder + "\\" + "expressionRxns.json", "w") as f:
        dump(list_of_expressions, f)

    list_of_significance = []
    for sig in significance:
        list_of_significance.append(sig.tolist())
    with open(main_folder + "\\" + "significance.json", "w") as f:
        dump(list_of_significance, f)

    def deep_convert_np_to_list(obj):
        if isinstance(obj, np.ndarray):
            return deep_convert_np_to_list(obj.tolist())
        elif isinstance(obj, list):
            return [deep_convert_np_to_list(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(deep_convert_np_to_list(item) for item in obj)
        else:
            return obj

    # incase structure is weird, do the following after taking from mat:
    # task_struct = scipy.io.loadmat(main_folder + "\\" + "taskStructure.mat")
    # task_struct = task_struct["taskStructure"]

    # new_task_struct = deep_convert_np_to_list(task_struct)
    # new_task_struct = new_task_struct[0][0]
    # new_task_struct = new_task_struct[0]
    # Usage:
    list_of_task_structure = deep_convert_np_to_list(task_structure)
    with open(main_folder + "\\" + "taskStructure.json", "w") as f:
        dump(list_of_task_structure, f)

    dump(list_of_task_structure, open(main_folder + "\\" + "taskStructure.json", "w"))

if __name__ == "__main__":
    placeholder = 5
    # main_folder = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing"
    # save_files_as_json(main_folder)
    main_folder_1 = "E:\\Git\\Metabolic_Task_Score\\Data\\Base MAT Files\\Consensus Model 1.17\\No Thresh"
    # save_files_as_json(main_folder_1)


def save_solver_results(
    instance: pe.ConcreteModel,
    result,
    solver,
    cobra_model: Model,
    filename: str,
    main_folder: str,
    expression_array: np.ndarray = None,
    old_expression_array: np.ndarray = None,
    epsilon: float = 1e-3,
) -> None:
    var_values = [
        (v.name, v.value) for v in instance.component_data_objects(pe.Var, active=True)
    ]
    if len(var_values)/2 > len(cobra_model.reactions):
        var_values = var_values[:2*len(cobra_model.reactions)]
    var_values = [
        (name.split("[")[0], int(name.split("[")[1].split("]")[0]), value)
        for name, value in var_values
    ]
    df_vars = DataFrame(var_values, columns=["Variable", "Index", "Value"])

    min_indices = df_vars.groupby("Variable")["Index"].min().reset_index()
    min_indices.columns = ["Variable", "MinIndex"]
    df_vars = df_vars.merge(min_indices, on="Variable")
    df_vars["AdjustedIndex"] = df_vars.apply(
        lambda row: row["Index"] + 1 if row["MinIndex"] == 0 else row["Index"], axis=1
    )
    df_vars = df_vars.pivot(
        index="AdjustedIndex", columns="Variable", values="Value"
    ).reset_index()
    # df_vars = df_vars.pivot(index='Index', columns='Variable', values='Value').reset_index()

    if cobra_model is not None:
        df_vars["Reaction Name"] = [reaction.name for reaction in cobra_model.reactions]
        df_vars["Reaction ID"] = [reaction.id for reaction in cobra_model.reactions]
        df_vars["Reaction Formula"] = [
            str(reaction.build_reaction_string(use_metabolite_names=True))
            for reaction in cobra_model.reactions
        ]
        df_vars["Reaction Formula Using IDs"] = [
            str(reaction.build_reaction_string(use_metabolite_names=False))
            for reaction in cobra_model.reactions
        ]
    if expression_array is not None:
        df_vars["Expression"] = old_expression_array
        df_vars["Modified_Expression"] = expression_array
    df_vars.drop(columns=["AdjustedIndex"], inplace=True)
    df_vars.insert(0, "Zeros", 0)

    columns_order = [
        "Zeros",
        "V_flux",
        "Expression",
        "Reaction ID",
        "Reaction Name",
        "Reaction Formula",
        "Reaction Formula Using IDs",
    ]
    other_columns = [col for col in df_vars.columns if col not in columns_order]
    df_vars = df_vars[columns_order + other_columns]

    solver_output = solver._log
    pattern = r'\((\d+\.\d+) work units\)'
    matches = findall(pattern, solver_output)
    second_work_units = 999
    if matches is not None and len(matches) >= 2:
        second_work_units = matches[1]

    result_dict = {
        "Termination condition": result.Solver.Termination_condition,
        "Work units": float(second_work_units),
        "Time": result.Solver.Time,
        "Objective": instance.objective(),
    }
    # df_result = DataFrame.from_records([result_dict]).transpose()
    df_result = DataFrame(list(result_dict.items()), columns=["Settings", "Values"])
    # df_settings = DataFrame(df_result[0].values, columns=['Settings'])
    # df_settings = df_settings.reset_index(drop =True)
    # df_result.reset_index(drop = True)
    df_settings = df_result

    df_final = concat([df_vars, df_settings], axis=1)

    location_to_save = path.join(main_folder, filename)
    save_excel_or_CSV_if_file_already_exist(df_final, location_to_save)

    columns_to_check = ["V", "Y", "V_flux", "Y_rxn_active"]
    df_filtered = df_vars[abs(df_vars["Y_rxn_active"]) > 0]
    df_filtered = df_filtered.reset_index(drop=True)

    for column in columns_to_check:
        if column in df_vars.columns:
            # df_filtered = df_vars[abs(df_vars[column]) >= epsilon]
            # df_filtered = df_filtered.reset_index(drop=True)
            df_filtered = df_filtered[abs(df_filtered[column]) >= epsilon]
            df_filtered = df_filtered.reset_index(drop=True)
            df_final_filtered = concat([df_filtered, df_settings], axis=1)
            location_to_save_filtered = path.join(main_folder, f"filtered_{filename}")
            save_excel_or_CSV_if_file_already_exist(
                df_final_filtered, location_to_save_filtered
            )
            # create file with 2 columns and header for esher
            df_esher = df_filtered[["Reaction ID", "V_flux"]]
            location_to_save_esher = path.join(main_folder, f"esher_{filename}")
            with_header = True
            make_CSV = True
            save_excel_or_CSV_if_file_already_exist(
                df_esher, location_to_save_esher, with_header, make_CSV
            )
            break



if __name__ == "__main__":
    placeholder = 5
    # main_folder = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing"
    # save_solver_results(instance, result, 'solver_results.csv', main_folder)


# add function for base creation of SV=0 LP constraint from stoichiometry matrix


def return_new_model_with_only_included_reactions(
    reactions_for_filtering: list,
) -> Model:
    new_model = Model()
    new_model.add_reactions(reactions_for_filtering)
    return new_model


def get_active_reactions_from_solved_model(
    instance,
    expression_old: [],
    expression: [],
    cobra_model: Model,
    **settings_kwargs_dict: dict,
) -> (dict, tuple):
    reactions_included_this_task = (
        {}
    )  # task_number: {reaction_id: (reaction_name, flux, reaction_formula)}
    if "epsilon" in settings_kwargs_dict:
        epsilon = settings_kwargs_dict["epsilon"]
    else:
        epsilon = 1e-3
    var_values = [
        (v.name, v.value) for v in instance.component_data_objects(Var, active=True)
    ]
    if len(var_values)/2 > len(cobra_model.reactions):
        var_values = var_values[:2*len(cobra_model.reactions)]
    var_values = [
        (name.split("[")[0], int(name.split("[")[1].split("]")[0]), value)
        for name, value in var_values
    ]
    df_vars = DataFrame(var_values, columns=["Variable", "Index", "Value"])
    df_vars = df_vars.pivot(
        index="Index", columns="Variable", values="Value"
    ).reset_index()
    reactions_for_filtering_model = []
    for rxn in range(1, len(cobra_model.reactions) + 1):
        if abs(df_vars["V_flux"][rxn - 1]) >= epsilon:
            name = cobra_model.reactions[rxn - 1].name
            id = cobra_model.reactions[rxn - 1].id
            reactions_for_filtering_model.append(cobra_model.reactions[rxn - 1])
            formula = cobra_model.reactions[rxn - 1].build_reaction_string(
                use_metabolite_names=True
            )
            flux = df_vars["V_flux"][rxn - 1]
            reactions_included_this_task[id] = (
                0,
                flux,
                expression_old[rxn - 1],
                id,
                name,
                formula,
                expression[rxn - 1],
            )
    amount_of_reactions = len(reactions_included_this_task.keys())
    amount_of_expressionless_reactions = len(
        [value for key, value in reactions_included_this_task.items() if value[2] == -1]
    )
    amount_of_reactions_with_expression = (
        amount_of_reactions - amount_of_expressionless_reactions
    )
    print(
        f"\nAmount of reactions with expression: {amount_of_reactions_with_expression}"
    )
    print(
        f"Amount of reactions without expression: {amount_of_expressionless_reactions}"
    )
    print(f"Total amount of reactions: {amount_of_reactions}")

    new_model = return_new_model_with_only_included_reactions(
        reactions_for_filtering_model
    )

    return (
        reactions_included_this_task,
        (
            amount_of_reactions,
            amount_of_expressionless_reactions,
            amount_of_reactions_with_expression,
        ),
        new_model,
    )
