from cobra import Model, Reaction
from cobra.util import create_stoichiometric_matrix
import numpy as np
from copy import deepcopy
from cobra.io import (
    load_matlab_model,
    save_matlab_model,
)
from json import dump
from os import makedirs
from time import time
import file_management as fn
from SRC.Special_sparse_matrix_fastCC_code.special_fastcc_implementation import fastcc
from scipy.io import loadmat

__all__ = [
    "close_exchange_reactions",
    "remove_exchange_reactions",
    "prepare_data_for_analysis",
    "open_exchange_reactions_for_task",
    "create_reactions_to_add_and_to_create_rev2irrev",
    "add_reactions_to_model",
    "convert_to_irreversible",
    "create_task_specific_reversible_and_irreversible_models",
]


def close_exchange_reactions(model: Model) -> Model:
    print("Closing exchange reactions")
    # print(len(model.exchanges))
    # print(len(model.sinks))
    # print(len(model.demands))

    for reaction in model.exchanges:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    return model


def remove_exchange_reactions(model: Model) -> Model:
    print("Removing exchange reactions")
    rxns_to_remove = [rxn for rxn in model.exchanges]
    for rxn in rxns_to_remove:
        model.remove_reactions(rxn)
    return model


def prepare_data_for_analysis(
    cobra_model: Model, expression: list, significance: list
) -> (Model, list):
    # folder = r"E:\Git\Metabolic_Task_Score\Data\Main_files\For_running\models\Consensus Model 1.17\No Thresh\\"
    # cobra_model = loadmat(folder + "modelCellfie.mat")
    # expression = loadmat(folder +"expressionRxnsCellfie.mat")["expressionRxns"]
    # significance = loadmat(folder+ "significanceCellfie.mat")["significance"]
    for expr_samples in expression:
        for idx, sample_value in enumerate(expr_samples):
            if not sample_value > 0:
                expr_samples[idx] = -1

    expression_new = deepcopy(expression)
    for idx in range(len(expression)):
        for idx2 in range(len(expression[idx])):
            expression_new[idx][idx2] = expression[idx][idx2] * significance[idx][idx2]
            if not expression_new[idx][idx2] > 0:
                expression_new[idx][idx2] = -1

    return cobra_model, expression_new


def open_exchange_reactions_for_task(
    original_model: Model, task_number: int, task_structure: list
) -> (Model, bool):
    '''
    This function takes the task_structure and opens the exchange reactions for the task_number
    :param original_model:
    :param task_number:
    :param task_structure:
    :return:
    '''
    pass_task = True
    # model = original_model.copy()
    model = original_model
    model.S = create_stoichiometric_matrix(model)
    task = []
    for task in task_structure:
        if task[0][0] == [str(task_number)]:
            task = task[0]
            break
    task_name = task[1]
    # task_input_metabolites = [item for sublist in task[5] for subsublist in sublist for item in subsublist]
    # task_input_metabolites_lb = np.array([item for sublist in task[6] for item in sublist])
    # task_input_metabolites_ub = np.array([item for sublist in task[7] for item in sublist])
    #
    # task_output_metabolites = [item for sublist in task[8] for subsublist in sublist for item in subsublist]
    # task_output_metabolites_lb = np.array([item for sublist in task[9] for item in sublist])
    # task_output_metabolites_ub = np.array([item for sublist in task[10] for item in sublist])

    # a = [x2+100 for x in [2,3,4] for x2 in range(x) ]

    task_input_metabolites = task[5]
    task_input_metabolites = [
        item
        for sublist in task_input_metabolites
        for subsublist in sublist
        for item in subsublist
    ]
    task_input_metabolites_lb = task[6]
    task_input_metabolites_lb = [
        item for sublist in task_input_metabolites_lb for item in sublist
    ]
    task_input_metabolites_ub = task[7]
    task_input_metabolites_ub = [
        item for sublist in task_input_metabolites_ub for item in sublist
    ]

    task_output_metabolites = task[8]
    task_output_metabolites = [
        item
        for sublist in task_output_metabolites
        for subsublist in sublist
        for item in subsublist
    ]
    task_output_metabolites_lb = task[9]
    task_output_metabolites_lb = [
        item for sublist in task_output_metabolites_lb for item in sublist
    ]
    task_output_metabolites_ub = task[10]
    task_output_metabolites_ub = [
        item for sublist in task_output_metabolites_ub for item in sublist
    ]
    rxn_subs = []
    for input_metabolite, lb, ub in zip(
        task_input_metabolites, task_input_metabolites_lb, task_input_metabolites_ub
    ):
        added_exchange_counter = 1
        if added_exchange_counter == 1:
            exchange_rxn = Reaction(f"temporary_exchange_{input_metabolite}")
            try:
                exchange_rxn.add_metabolites(
                    {model.metabolites.get_by_id(input_metabolite): 1}
                )
            except:
                print(f"Could not find metabolite {input_metabolite}")
                pass_task = False
                continue
            if input_metabolite in task_output_metabolites:
                exchange_rxn.lower_bound = -1000
                exchange_rxn.upper_bound = 1000
            else:
                exchange_rxn.lower_bound = lb
                exchange_rxn.upper_bound = ub
            model.add_reactions([exchange_rxn])
            rxn_subs.append(f"temporary_exchange_{input_metabolite}")

    model_metabolites = [met.id.upper() for met in model.metabolites]
    task_inputs_upper = [input.upper() for input in task_input_metabolites]
    indices = [i for i, x in enumerate(model_metabolites) if x in task_inputs_upper]

    if len(indices) != len(task_input_metabolites):
        print(f"ERROR: Could not find all inputs in {task_name}")
        print(f"Could not find all inputs")
        pass_task = False
    if len(indices) != len(np.unique(indices)):
        print(
            f'The constraints on some input(s) in "[{task_name}]" are defined more than one time'
        )
        pass_task = False

    rxn_prods = []
    for output_metabolite, lb, ub in zip(
        task_output_metabolites, task_output_metabolites_lb, task_output_metabolites_ub
    ):
        if output_metabolite in task_input_metabolites:
            # find reaction created for input_metabolite
            for rxn_sub in rxn_subs:
                if output_metabolite in rxn_sub:
                    exchange_rxn = model.reactions.get_by_id(rxn_sub)
                    exchange_rxn.lower_bound = -ub
                    exchange_rxn.upper_bound = ub
                    break
            print(f"Output metabolite {output_metabolite} is also an input metabolite")
            continue

        added_exchange_counter = 1

        if added_exchange_counter == 1:
            exchange_rxn = Reaction(f"temporary_exchange_{output_metabolite}")
            try:
                exchange_rxn.add_metabolites(
                    {model.metabolites.get_by_id(output_metabolite): -1}
                )
            except:
                print(f"Could not find metabolite {output_metabolite}")
                pass_task = False
                continue
            exchange_rxn.lower_bound = lb
            exchange_rxn.upper_bound = ub
            model.add_reactions([exchange_rxn])
            rxn_prods.append(f"temporary_exchange_{output_metabolite}")

    model_metabolites = [met.id.upper() for met in model.metabolites]
    task_output_upper = [output.upper() for output in task_output_metabolites]
    indices = [i for i, x in enumerate(model_metabolites) if x in task_output_upper]

    if len(indices) != len(task_output_metabolites):
        print(f"ERROR: Could not find all inputs in {task_name}")
        print(f"Could not find all inputs")
        pass_task = False
    if len(indices) != len(np.unique(indices)):
        print(
            f'The constraints on some input(s) in "[{task_name}]" are defined more than one time'
        )
        pass_task = False

    if pass_task:
        solution = model.optimize()
        if solution.status != "optimal":
            print(f"Task {task_number} did not converge, namely {solution.status}")
            pass_task = False
        else:
            print(f"Task {task_number} passed")

    return model, pass_task


def create_reactions_to_add_and_to_create_rev2irrev(
    model_old: Model,
) -> (list, list):
    model = model_old
    reactions_to_modify = [
        rxn for rxn in model.reactions if rxn.lower_bound < 0 and rxn.upper_bound > 0
    ]
    to_create_rev2irrev = [[idx + 1] for idx in range(len(model.reactions))]
    rxn_subs = [0] * len(reactions_to_modify)
    idx_2_counter = 0
    length_of_model = len(model.reactions)
    for idx, rxn in enumerate(model.reactions):
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            new_reaction = Reaction(rxn.id + "_r")
            new_reaction.subsystem = rxn.subsystem
            new_reaction.gene_reaction_rule = rxn.gene_reaction_rule
            # new_reaction.gene_name_reaction_rule = rxn.gene_name_reaction_rule
            # new_reaction.genes = rxn.genes
            new_reaction.name = rxn.name
            new_reaction.bounds = (0, -rxn.lower_bound)
            new_reaction.add_metabolites(
                dict((met, -1 * coeff) for met, coeff in rxn.metabolites.items())
            )
            rxn_subs[idx_2_counter] = new_reaction
            to_create_rev2irrev[idx].append(idx_2_counter + 1 + length_of_model)
            idx_2_counter += 1

    return rxn_subs, to_create_rev2irrev


def add_reactions_to_model(
    model: Model,
    amount_of_temporary_reactions_added: int,
    rxn_subs_ori: list,
    to_create_rev2irrev_ori: list,
) -> (Model, list):
    rxn_subs = deepcopy(rxn_subs_ori)
    to_create_rev2irrev = deepcopy(to_create_rev2irrev_ori)
    length_of_model = len(model.reactions)
    for value in to_create_rev2irrev:
        if len(value) > 1:
            value[1] += amount_of_temporary_reactions_added

    for i in range(amount_of_temporary_reactions_added):
        to_create_rev2irrev.append(
            [1 + length_of_model - amount_of_temporary_reactions_added + i]
        )

    for rxn_new in rxn_subs:
        rxn = model.reactions.get_by_id(rxn_new.id[:-2])
        rxn.id = rxn.id + "_f"
        rxn.bounds = (0, rxn.upper_bound)
    model.add_reactions(rxn_subs)

    if amount_of_temporary_reactions_added > 0:
        for idx, rxn in enumerate(
            model.reactions[-amount_of_temporary_reactions_added:]
        ):
            if rxn.lower_bound < 0 and rxn.upper_bound > 0:
                new_reaction = Reaction(rxn.id + "_r")
                new_reaction.gene_reaction_rule = rxn.gene_reaction_rule
                # new_reaction.gene_name_reaction_rule = rxn.gene_name_reaction_rule
                # new_reaction.genes = rxn.genes
                new_reaction.name = rxn.name
                new_reaction.subsystem = rxn.subsystem
                new_reaction.bounds = (0, -rxn.lower_bound)
                new_reaction.add_metabolites(
                    dict((met, -1 * coeff) for met, coeff in rxn.metabolites.items())
                )
                rxn_subs.append(new_reaction)
                to_create_rev2irrev[
                    length_of_model - amount_of_temporary_reactions_added + idx
                ].append(length_of_model + idx + 1)

                rxn.id = rxn.id + "_f"
                rxn.bounds = (0, rxn.upper_bound)

    return model, to_create_rev2irrev


def convert_to_irreversible(model_old: Model, copy_bool: bool = False) -> (Model, list):
    print("Converting to irreversible")
    if copy_bool:
        model = model_old.copy()
    else:
        model = model_old
    # model.solver = "gurobi"
    reactions_to_modify = [
        rxn for rxn in model.reactions if rxn.lower_bound < 0 and rxn.upper_bound > 0
    ]
    rev2irrev = [[idx + 1] for idx in range(len(model.reactions))]

    rxn_subs = [0] * len(reactions_to_modify)
    length_of_model = len(model.reactions)
    for idx, rxn in enumerate(reactions_to_modify):
        # print(idx, rxn)
        new_reaction = Reaction(rxn.id + "_r")
        new_reaction.subsystem = rxn.subsystem
        new_reaction.bounds = (0, -rxn.lower_bound)
        new_reaction.add_metabolites(
            dict((met, -1 * coeff) for met, coeff in rxn.metabolites.items())
        )

        rxn.id = rxn.id + "_f"
        rxn.lower_bound = 0

        # model.add_reactions([new_reaction])
        old_reaction_index = model.reactions.index(model.reactions.get_by_id(rxn.id))
        # new_reaction_index = model.reactions.index(model.reactions.get_by_id(new_reaction.id))
        rev2irrev[old_reaction_index].append(idx + 1 + length_of_model)
        rxn_subs[idx] = new_reaction
    # print("Added new reactions", rev2irrev)
    # print("Added new reactions", rxn_subs)
    model.add_reactions(rxn_subs)

    return model, rev2irrev


def create_task_specific_reversible_and_irreversible_models(
    main_folder: str,
    do_fastCC_filtering: bool = False,
    create_new_model_files: bool = False,
    expression_type: str = "local_threshold",
):
    time_start = time()
    model = load_matlab_model(main_folder + "\\" + "modelCellfie.mat")
    print(f"Model loaded in {time() - time_start} seconds")
    expression, significance, task_structure = fn.load_json_data(main_folder, expression_type)
    model, expression = prepare_data_for_analysis(model, expression, significance)
    if expression_type == "local_threshold":
        dump(expression, open(main_folder + "\\" + "full_expression_reversible.json", "w"))
    elif expression_type == "global_percentile_90":
        dump(
            expression,
            open(main_folder + "\\" + "full_expression_reversible_percentile_90.json", "w"),
        )
    elif expression_type == "global_no_threshold":
        dump(
            expression,
            open(main_folder + "\\" + "full_expression_reversible_no_threshold.json", "w"),
        )

    model = close_exchange_reactions(model)
    # model = remove_exchange_reactions(model)
    rxn_subs, to_create_rev2irrev = create_reactions_to_add_and_to_create_rev2irrev(
        model
    )
    dump(
        to_create_rev2irrev,
        open(main_folder + "\\" + "full_to_create_rev2irrev.json", "w"),
    )
    irreversible_model = model.copy()
    irreversible_model, rev2irrev2 = add_reactions_to_model(
        irreversible_model, 0, rxn_subs, to_create_rev2irrev
    )
    amount_of_added_reactions_base = len(irreversible_model.reactions) - len(
        model.reactions
    )
    expression_irrev = deepcopy(expression)
    for idx in range(amount_of_added_reactions_base):
        expression_irrev.append(["Fill"] * len(expression_irrev[0]))
    dump(
        expression_irrev,
        open(main_folder + "\\" + "full_expression_irreversible.json", "w"),
    )

    for task_number in range(1010,1011):
        print("\n_________\n", f"Task {task_number}\n")
        time_start = time()
        task_specific_model = model.copy()
        print(f"Time to copy model and expression: {time() - time_start} seconds")

        task_specific_model, pass_task = open_exchange_reactions_for_task(
            task_specific_model, task_number, task_structure
        )
        amount_of_reactions_added = len(task_specific_model.reactions) - len(
            model.reactions
        )

        makedirs(main_folder + "\\" + f"Task_{task_number:03}", exist_ok=True)
        dump(
            [pass_task, amount_of_reactions_added],
            open(
                main_folder
                + "\\"
                + f"Task_{task_number:03}"
                + f"\\pass_task_and_reactions_added_task{task_number:03}.json",
                "w",
            ),
        )
        if pass_task:
            # if os.path.exists(main_folder + "\\" + f"Task_{task_number:03}" + f"\\rev2irrev_task{task_number:03}.json"):
            #     continue
            save_matlab_model(
                task_specific_model,
                main_folder
                + "\\"
                + f"Task_{task_number:03}"
                + f"\\task_specific_model_task{task_number:03}.mat",
            )
            if do_fastCC_filtering:
                time_start = time()
                fast_cc_reversible_model, indices_removed = fastcc(
                    task_specific_model, 1.0, 1e-4
                )
                print(f"Time for fastCC: {time() - time_start} seconds")
                dump(
                    indices_removed,
                    open(
                        main_folder
                        + "\\"
                        + f"Task_{task_number:03}"
                        + f"\\indices_removed_fastCC_task{task_number:03}.json",
                        "w",
                    ),
                )
                ### reactions to remove to adjust expression

                # save new model
                # save expression information as list of indices (or create new rev2irrev

                # create irreversible model from fast_cc_reversible_model

            task_specific_model, rev2irrev = add_reactions_to_model(
                task_specific_model,
                amount_of_reactions_added,
                rxn_subs,
                to_create_rev2irrev,
            )
            dump(
                rev2irrev,
                open(
                    main_folder
                    + "\\"
                    + f"Task_{task_number:03}"
                    + f"\\rev2irrev_task{task_number:03}.json",
                    "w",
                ),
            )
            save_matlab_model(
                task_specific_model,
                main_folder
                + "\\"
                + f"Task_{task_number:03}"
                + f"\\irrev_model_task{task_number:03}.mat",
            )

    # find model in basefiles, find taskstructure in basefiles, create a folder for each task in tasks


if __name__ == "__main__":
    placeholder = 5
    create_task_specific_reversible_and_irreversible_models(
        "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\Base_files",
        False,
        False,
        "global_no_threshold"
    )
    # prepare_data_for_analysis(None,None,None)

