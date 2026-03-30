import networkx as nx
from mip import Model, xsum, BINARY, INTEGER, minimize
from pyomo.opt import SolverFactory
from pyomo.environ import (
    ConcreteModel,
    Var,
    Binary,
    Constraint,
    Objective,
    minimize,
    maximize,
    Reals,
    ConstraintList,
    NonNegativeReals,
    NonNegativeIntegers,
    Binary,
    Set,
    Param,
    RangeSet,
    Integers
)
import matplotlib.pyplot as plt
from collections import deque
from netgraph import Graph
import threading
from cobra import Model as cobrapy_model
from cobra.util import create_stoichiometric_matrix
from cobra.io import load_json_model
import numpy as np
from typing import Dict, Tuple
import pandas as pd
import py4cytoscape as p4c
# # toy model
# input_metabolites = ["input_met_1", "input_met_2", "input_met_3"]
# output_metabolites = ["output_met_1", "output_met_2"]
# internal_metabolites = ["B1", "B2", "B3", "B4", "B5", "B6"]
# exchange_metabolites = input_metabolites + output_metabolites
# metabolites = internal_metabolites + exchange_metabolites
# to_remove_metabolites = ["B1", "B3", "B5"]
# not_removed_metabolites = [met for met in metabolites if met not in to_remove_metabolites]
#
# input_reactions = ["input_1", "input_2", "input_3"]
# output_reactions = ["output_1", "output_2"]
# internal_reactions = ["A1", "A2", "A3", "A4", "A5", "A6", "A7"]
# exchange_reactions = input_reactions + output_reactions
# reactions = internal_reactions + exchange_reactions
#
# input_edges = [("input_1", "input_met_1"), ("input_2", "input_met_2"), ("input_3", "input_met_3"),
#                ("input_met_1", "A1"), ("input_met_2", "A2"), ("input_met_3", "A3"), ("input_met_1", "A4")
#                ]
# output_edges = [("output_met_1", "output_1"), ("output_met_2", "output_2"),
#                 ("A6", "output_met_2"), ("A7", "output_met_1"), ("A7", "output_met_2")]
# # for internal edges:
# # nodes to be rmeoved (and their corresponding edges) B1, B3, B5
# # A1 b1, A1 b2, A2 B2, A2 B3, A3 B3, A3 B6, B1 A4, B2 A4, B2 A5, B3 A5, A4 B4, A4 B5, A5 B4, A5 B5, B4 A6, B5 A7, B6 A7
# internal_edges = [
#     ("A1", "B1"), ("A1", "B2"), ("A2", "B2"), ("A2", "B3"), ("A3", "B3"), ("A3", "B6"),
#     ("B1", "A4"), ("B2", "A4"), ("B2", "A5"), ("B3", "A5"), ("A4", "B4"), ("A4", "B5"),
#     ("A5", "B4"), ("A5", "B5"), ("B4", "A6"), ("B5", "A7"), ("B6", "A7")
# ]
# edges = internal_edges + input_edges + output_edges
# to_remove_edges = [edge for edge in edges if edge[0] in to_remove_metabolites or edge[1] in to_remove_metabolites]
# not_removed_edges = [edge for edge in edges if edge not in to_remove_edges]

# add code to decide levels based on edge connections starting from input nodes
levels_dict = {
    1: ["input_1", "input_2", "input_3"],
    2: ["input_met_1", "input_met_2", "input_met_3"],
    3: ["A1", "A2", "A3"],
    4: ["B1", "B2", "B3", "B6"],
    5: ["A4", "A5"],
    6: ["B4", "B5"],
    7: ["A6", "A7"],
    8: ["output_met_1", "output_met_2"],
    9: ["output_1", "output_2"]
}



def identify_stoichiometric_matrix_of_sample_task_filtered_model(
        model: Model,
        df_in: pd.DataFrame,
        task_number: int,
        sample_number: int,
) -> (pd.DataFrame, dict):

    df = select_sample_and_task_from_dataframe(df_in, task_number, sample_number)
    reactions = df["reaction_name"].unique()
    dtype = np.float64
    stoichiometric_matrix = np.zeros((len(model.metabolites), len(model.reactions)), dtype=dtype)
    m_ind = model.metabolites.index
    r_ind = model.reactions.index

    for reaction in model.reactions:
        if reaction.id in reactions:
            for metabolite, stoich in reaction.metabolites.items():
                if df[df["reaction_name"] == reaction.id]["flux"].iloc[0] > 0:
                    stoichiometric_matrix[m_ind(metabolite), r_ind(reaction)] = stoich
                elif df[df["reaction_name"] == reaction.id]["flux"].iloc[0] < 0:
                    stoichiometric_matrix[m_ind(metabolite), r_ind(reaction)] = -stoich

    exchange_reactions = [reaction for reaction in reactions if reaction.startswith("temporary_exchange_")]
    input_exchange_reactions = [reaction for reaction in exchange_reactions if (
        df[df["reaction_name"] == reaction]["reaction_formula"].str.startswith(" -->").iloc[0]
    )]
    ### TODO index of reactions doesnt follow the same order as shown, reaction at index (counting from 0 in both cases) 16 in excel and model is 27 in stoichiometric_matrix
    ### TODO additionally appears to be metabolites with no input or output values such as met at index 2, this is fixed by changing the stoichiometric matrix if flux is negative
    output_exchange_reactions = [reaction for reaction in exchange_reactions if reaction not in input_exchange_reactions]
    input_exchange_reaction_metabolites = [reaction_name.split("_")[2] for reaction_name in input_exchange_reactions]
    output_exchange_reaction_metabolites = [reaction_name.split("_")[2] for reaction_name in output_exchange_reactions]
    input_exchange_reaction_metabolite_indices = [i for i, m in enumerate(model.metabolites) if m.id in input_exchange_reaction_metabolites]
    output_exchange_reaction_metabolite_indices = [i for i, m in enumerate(model.metabolites) if m.id in output_exchange_reaction_metabolites]

    for idx, reaction in enumerate(input_exchange_reactions):
        stoichiometric_matrix = np.hstack([stoichiometric_matrix, np.zeros((stoichiometric_matrix.shape[0], 1))])
        stoichiometric_matrix[
            input_exchange_reaction_metabolite_indices[idx],
            np.shape(stoichiometric_matrix)[1]-1] = 1
    for idx, reaction in enumerate(output_exchange_reactions):
        stoichiometric_matrix = np.hstack([stoichiometric_matrix, np.zeros((stoichiometric_matrix.shape[0], 1))])
        stoichiometric_matrix[
            output_exchange_reaction_metabolite_indices[idx],
            np.shape(stoichiometric_matrix)[1]-1] = -1

    metabolite_indices = np.where(~np.all(stoichiometric_matrix == 0, axis=1))[0]
    stoichiometric_matrix = stoichiometric_matrix[~np.all(stoichiometric_matrix == 0, axis=1)]
    stoichiometric_matrix = stoichiometric_matrix[:, ~np.all(stoichiometric_matrix == 0, axis=0)]

    reaction_names = [0] * len(reactions)
    reaction_index_dict = {}
    reaction_name_index_dict = {}
    temp_reaction_counter = 0
    for idx3, reaction in enumerate(reactions):
        not_found = True
        idx3 = idx3 - temp_reaction_counter
        for reaction_model in model.reactions:
            if reaction_model.id == reaction:
                reaction_names[idx3] = reaction_model.name
                reaction_index_dict[idx3] = reaction_model.id
                reaction_name_index_dict[idx3] = reaction_model.name
                not_found = False
                break
        if not_found:
            idx_to_put_exchange_reaction = len(reactions) - len(exchange_reactions) + temp_reaction_counter
            print(f"Reaction {reaction} not found in model")
            reaction_names[idx_to_put_exchange_reaction] = reaction
            reaction_index_dict[idx_to_put_exchange_reaction] = reaction
            reaction_name_index_dict[idx_to_put_exchange_reaction] = reaction
            temp_reaction_counter += 1

    metabolite_ids = [metabolite.id for idx, metabolite in enumerate(model.metabolites) if idx in metabolite_indices]
    metabolite_names = [metabolite.name for metabolite in model.metabolites if metabolite.id in metabolite_ids]
    metabolite_index_dict = {idx: metabolite for idx, metabolite in enumerate(metabolite_ids)}
    metabolite_name_index_dict = {idx: metabolite for idx, metabolite in enumerate(metabolite_names)}


    return (stoichiometric_matrix, reaction_index_dict, reaction_name_index_dict, metabolite_index_dict, metabolite_name_index_dict)

def select_sample_and_task_from_dataframe(
        df_in: pd.DataFrame, task_number: int, sample_number: int
) -> pd.DataFrame:
    df = df_in.copy(deep=True)
    df = df[df["sample"] == sample_number]
    df = df[df["task"] == task_number]
    return df
def check_if_direct_path_possible(reaction_map: dict, stoichiometric_matrix: pd.DataFrame, reaction_index_dict: dict) -> bool:
    # each reaction should have an output metabolite that is the input metabolite in at least 1 reaction
    # each reaction should also have an input metabolite that is the output metabolite in at least 1 reaction
    all_output_metabolites = [output_metabolite for input_metabolite, output_metabolite in reaction_map.values()]
    all_input_metabolites = [input_metabolite for input_metabolite, output_metabolite in reaction_map.values()]

    for reaction, (input_metabolite, output_metabolite) in reaction_map.items():
        if input_metabolite not in all_output_metabolites:
            print(f"Reaction {reaction} has input metabolite {reaction_index_dict[input_metabolite]} that is not an output metabolite in any reaction")
            # return False
        if output_metabolite not in all_input_metabolites:
            print(f"Reaction {reaction} has output metabolite {reaction_index_dict[output_metabolite]} that is not an input metabolite in any reaction")
            # return False
    return True
def calculate_direct_path_using_least_connected_metabolites(
        model: Model, df_in: pd.DataFrame, task_number: int, sample_number: int
) -> pd.DataFrame:
    (filtered_stoichiometric_matrix, reaction_index_dict,
     reaction_name_index_dict, metabolite_index_dict,
     metabolite_name_index_dict)= identify_stoichiometric_matrix_of_sample_task_filtered_model(
        model, df_in, task_number, sample_number
    )
    # for each row in the stoichiometric matrix, caulcaute the amount of non-zero columns
    metabolite_number_of_reactions_dict = {
        idx: len(np.nonzero(row)[0]) for idx, row in enumerate(filtered_stoichiometric_matrix)
    }
    reaction_number_of_metabolites_dict = {
        idx: (np.sum(row > 0), np.sum(row < 0)) for idx, row in enumerate(filtered_stoichiometric_matrix.T)
    }
    df = select_sample_and_task_from_dataframe(df_in, task_number, sample_number)
    reactions = df["reaction_name"].unique()
    already_chosen_metabolites = []
    full_reaction_map = {}
    for idx, reaction in enumerate(reactions):
        metabolite_input_indices_in_reaction = np.where(filtered_stoichiometric_matrix[:, idx] < 0)
        metabolite_output_indices_in_reaction = np.where(filtered_stoichiometric_matrix[:, idx] > 0)

        chosen_metabolite = None
        # check if the metabolites are already chosen
        if len(metabolite_output_indices_in_reaction) == 0:
            pass
        elif len(set(metabolite_input_indices_in_reaction[0]).intersection(already_chosen_metabolites)) > 0:
            chosen_metabolite = list(set(metabolite_input_indices_in_reaction[0]).intersection(already_chosen_metabolites))[0]
        # otherwise choose the metabolite with the least amount of reactions
        else:
            minimum_number_of_connections = 10000000
            dict_subset = {key: value[0] for key, value in reaction_number_of_metabolites_dict.items() if key in metabolite_input_indices_in_reaction[0]}
            for key, value in dict_subset.items():
                if value < minimum_number_of_connections:
                    minimum_number_of_connections = value
                    chosen_metabolite = key
        input_metabolite = chosen_metabolite
        already_chosen_metabolites.append(input_metabolite)

        chosen_metabolite = None
        if len(metabolite_output_indices_in_reaction) == 0:
            pass
        elif len(set(metabolite_output_indices_in_reaction[0]).intersection(already_chosen_metabolites)) > 0:
            chosen_metabolite = list(set(metabolite_output_indices_in_reaction[0]).intersection(already_chosen_metabolites))[0]
        else:
            minimum_number_of_connections = 10000000
            dict_subset = {key: value[1] for key, value in reaction_number_of_metabolites_dict.items() if key in metabolite_output_indices_in_reaction[0]}
            for key, value in dict_subset.items():
                if value < minimum_number_of_connections:
                    minimum_number_of_connections = value
                    chosen_metabolite = key
        output_metabolite = chosen_metabolite
        already_chosen_metabolites.append(output_metabolite)
        full_reaction_map[reaction] = (input_metabolite, output_metabolite)

    # check if all reactions can be connected using the input and output metabolites
    if check_if_direct_path_possible(full_reaction_map, filtered_stoichiometric_matrix,reaction_index_dict):
        print("Direct path is possible")

    ### logic:
    # every node is at exactly one grid point
    # every grid point has exactly one node
    # every reaction is at an even grid point
    # every reaction must be placed
    # metabolites do not have to be placed
    # if placed metabolites are at odd grid points
    # input reactions are at the first column
    # output reactions are at the last column
    # every reaction has a flow through every node
    # distance between gridpoints if a node is placed at the grid point and an edge between the nodes exists
    # minimize the sum of the distances between the grid points

def formulate_shortest_distance_on_grid_MILP_2(
        edges_ : [(str, str)],
        input_reactions : [str],
        output_reactions : [str],
        metabolites : [str],
        reactions : [str],
        reachable_reactions : {str: [str]},
) -> Dict: #Dict[str: Tuple[int, int]]\


    print("Creating model")
    model = ConcreteModel()
    max_columns = len(reactions) + len(metabolites) - len(input_reactions) - len(output_reactions)
    max_rows = max((len(reactions), len(metabolites)))
    nodes = reactions + metabolites
    rows = range(max_rows)

    model.reactions = Set(initialize=reactions)
    model.input_reactions = Set(initialize=input_reactions)
    model.output_reactions = Set(initialize=output_reactions)
    model.metabolites = Set(initialize=metabolites)
    model.nodes = Set(initialize=nodes)
    model.edges = Set(initialize=edges_)
    model.even_columns = RangeSet(2, max_columns+2, 2)
    model.odd_columns = RangeSet(3, max_columns+3, 2)
    model.columns = RangeSet(1, max_columns+1)
    model.rows = Set(initialize=rows)

    model.reactions_at_grid_point = Var(model.reactions, model.even_columns, model.rows, domain=Binary)
    model.metabolites_at_grid_point = Var(model.metabolites, model.odd_columns, model.rows, domain=Binary)
    model.reaction_column_position = Var(model.reactions, domain=Integers)
    model.metabolite_column_position = Var(model.metabolites, domain=Integers)
    model.metabolite_row_position = Var(model.metabolites, domain=Integers)
    model.reaction_row_position = Var(model.reactions, domain=Integers)
    model.horizontal_distance = Var(model.edges, domain=Reals)
    model.true_horizontal_distance = Var(model.edges, domain=NonNegativeReals)
    model.horizontal_slack_variable = Var(model.edges, domain=Reals)
    model.vertical_distance = Var(model.edges, domain=Reals)
    model.true_vertical_distance = Var(model.edges, domain=NonNegativeReals)
    model.vertical_slack_variable = Var(model.edges, domain=Reals)
    model.metabolite_not_placed = Var(model.metabolites, domain=Binary)
    model.flow_from_reaction_through_edge = Var(model.reactions, model.edges, domain=Binary)

    big_M = max(max_columns+3, max_rows+3)

    def edge_flow_constraint_1(model, reaction, node1, node2):
        edge = (node1, node2)
        if len([model.flow_from_reaction_through_edge[reaction, edge_] for edge_ in model.edges if edge_[1] == node1]) > 0:
            return sum(model.flow_from_reaction_through_edge[reaction, edge_] for edge_ in model.edges if edge_[1] == node1) >=\
                model.flow_from_reaction_through_edge[reaction, edge]
        else:
            return Constraint.Skip

    model.edge_flow_constraint = Constraint(model.reactions, model.edges, rule=edge_flow_constraint_1)

    def edge_flow_inactive_metabolite_constraint(model, reaction, node1, node2):
        edge = (node1, node2)
        if node1 in model.metabolites:
            return model.flow_from_reaction_through_edge[reaction, edge] <= 1 - model.metabolite_not_placed[node1]
        elif node2 in model.metabolites:
            return model.flow_from_reaction_through_edge[reaction, edge] <= 1 - model.metabolite_not_placed[node2]
        return Constraint.Skip

    model.edge_flow_inactive_metabolite_constraint = Constraint(model.reactions, model.edges, rule=edge_flow_inactive_metabolite_constraint)

    def flow_receive_constraint(model, reaction, reaction2):
        if (reaction2 in reachable_reactions[reaction]):
            return sum(model.flow_from_reaction_through_edge[reaction, edge] for edge in model.edges if edge[1] == reaction2) >= 1
        else:
            return Constraint.Skip

    model.flow_receive_constraint = Constraint(model.reactions, model.reactions, rule=flow_receive_constraint)

    def calculate_horizontal_distance(model, node1, node2):
        edge = (node1, node2)
        if node1 in model.reactions and node2 in model.metabolites:
            return model.horizontal_distance[edge] >= (model.metabolite_column_position[node2] - model.reaction_column_position[node1]) - big_M * model.metabolite_not_placed[node2]
        elif node1 in model.metabolites and node2 in model.reactions:
            return model.horizontal_distance[edge] >= (model.reaction_column_position[node2] - model.metabolite_column_position[node1]) - big_M * model.metabolite_not_placed[node1]
        return model.horizontal_distance[edge] == 0

    def calculate_vertical_distance(model, node1, node2):
        edge = (node1, node2)
        if node1 in model.reactions and node2 in model.metabolites:
            return model.vertical_distance[edge] >= (model.metabolite_row_position[node2] - model.reaction_row_position[node1]) - big_M * model.metabolite_not_placed[node2]
        elif node1 in model.metabolites and node2 in model.reactions:
            return model.vertical_distance[edge] >= (model.reaction_row_position[node2] - model.metabolite_row_position[node1]) - big_M * model.metabolite_not_placed[node1]
        return model.vertical_distance[edge] == 0

    def calculate_horizontal_distance2(model, node1, node2):
        edge = (node1, node2)
        if node1 in model.reactions and node2 in model.metabolites:
            return model.horizontal_distance[edge] <= (model.metabolite_column_position[node2] - model.reaction_column_position[node1]) + big_M * model.metabolite_not_placed[node2]
        elif node1 in model.metabolites and node2 in model.reactions:
            return model.horizontal_distance[edge] <= (model.reaction_column_position[node2] - model.metabolite_column_position[node1])  + big_M * model.metabolite_not_placed[node1]
        return model.horizontal_distance[edge] == 0

    model.calculate_horizontal_distance2 = Constraint(model.edges, rule=calculate_horizontal_distance2)

    def calculate_vertical_distance2(model, node1, node2):
        edge = (node1, node2)
        if node1 in model.reactions and node2 in model.metabolites:
            return model.vertical_distance[edge] <= (model.metabolite_row_position[node2] - model.reaction_row_position[node1]) + big_M * model.metabolite_not_placed[node2]
        elif node1 in model.metabolites and node2 in model.reactions:
            return model.vertical_distance[edge] <= (model.reaction_row_position[node2] - model.metabolite_row_position[node1]) + big_M * model.metabolite_not_placed[node1]
        return model.vertical_distance[edge] == 0


    model.calculate_vertical_distance2 = Constraint(model.edges, rule=calculate_vertical_distance2)
    def metabolite_placement_rule(model, metabolite):
        return model.metabolite_not_placed[metabolite] == 1 - sum(model.metabolites_at_grid_point[metabolite, column, row]
                                                          for column in model.odd_columns for row in model.rows)

    model.metabolite_placement_constraint = Constraint(model.metabolites, rule=metabolite_placement_rule)
    def set_horizontal_slack_variable(model, node1, node2):
        edge = (node1, node2)
        return model.horizontal_slack_variable[edge] + model.horizontal_distance[edge] == 0 # set to 1 to incentivize positive distance

    def set_true_horizontal_distance1(model, node1, node2):
        edge = (node1, node2)
        return model.true_horizontal_distance[edge] >= model.horizontal_distance[edge]
    def set_true_horizontal_distance2(model, node1, node2):
        edge = (node1, node2)
        return model.true_horizontal_distance[edge] >= model.horizontal_slack_variable[edge]

    def set_vertical_slack_variable(model, node1, node2):
        edge = (node1, node2)
        return model.vertical_slack_variable[edge] + model.vertical_distance[edge] == 0

    def set_true_vertical_distance1(model, node1, node2):
        edge = (node1, node2)
        return model.true_vertical_distance[edge] >= model.vertical_distance[edge]

    def set_true_vertical_distance2(model, node1, node2):
        edge = (node1, node2)
        return model.true_vertical_distance[edge] >= model.vertical_slack_variable[edge]



    def set_reaction_column_position(model, reaction):
        return model.reaction_column_position[reaction] == sum(column * model.reactions_at_grid_point[reaction, column, row]
                                                               for column in model.even_columns
                                                               for row in model.rows)

    def set_metabolite_column_position(model, metabolite):
        return model.metabolite_column_position[metabolite] == sum(column * model.metabolites_at_grid_point[metabolite, column, row]
                                                                   for column in model.odd_columns
                                                                   for row in model.rows)

    def set_metabolite_row_position(model, metabolite):
        return model.metabolite_row_position[metabolite] == sum(row * model.metabolites_at_grid_point[metabolite, column, row]
                                                               for column in model.odd_columns
                                                               for row in model.rows)

    def set_reaction_row_position(model, reaction):
        return model.reaction_row_position[reaction] == sum(row * model.reactions_at_grid_point[reaction, column, row]
                                                               for column in model.even_columns
                                                               for row in model.rows)

    def reaction_at_exactly_one_grid_point(model, reaction):
        return sum(model.reactions_at_grid_point[reaction, column, row] for column in model.even_columns for row in model.rows) == 1

    def metabolite_at_maximally_one_grid_point(model, metabolite):
        return sum(model.metabolites_at_grid_point[metabolite, column, row] for column in model.odd_columns for row in model.rows) <= 1

    def at_most_one_reaction_per_even_column_gridpoint(model, column, row):
        return sum(model.reactions_at_grid_point[reaction, column, row] for reaction in model.reactions) <= 1

    def at_most_one_metabolite_per_odd_column_gridpoint(model, column, row):
        return sum(model.metabolites_at_grid_point[metabolite, column, row] for metabolite in model.metabolites) <= 1

    def input_reactions_before_all_other_reactions_rule(model, input_reaction, regular_reaction):
        if regular_reaction not in model.input_reactions:
            # Enforce the constraint only if the regular_reaction is not an input reaction
            return model.reaction_column_position[input_reaction] <= model.reaction_column_position[regular_reaction] - 1
        else:
            return Constraint.Skip
    model.input_reactions_before_all_other_reactions_constraint = Constraint(model.input_reactions, model.reactions, rule= input_reactions_before_all_other_reactions_rule)
    def output_reaction_link_rule(model, output_reaction, regular_reaction):
        if regular_reaction not in model.output_reactions:
            # Enforce the constraint only if the regular_reaction is not an output reaction
            return model.reaction_column_position[output_reaction] >= model.reaction_column_position[regular_reaction] + 1
        else:
            return Constraint.Skip  # Skip the constraint if the regular_reaction is in output_reactions

    model.output_reaction_link_constraint = Constraint(model.output_reactions, model.reactions, rule=output_reaction_link_rule)
    model.reaction_at_exactly_one_grid_point = Constraint(model.reactions, rule=reaction_at_exactly_one_grid_point)
    model.metabolite_at_maximally_one_grid_point = Constraint(model.metabolites, rule=metabolite_at_maximally_one_grid_point)
    model.at_most_one_reaction_per_even_column_gridpoint = Constraint(model.even_columns, model.rows, rule=at_most_one_reaction_per_even_column_gridpoint)
    model.at_most_one_metabolite_per_odd_column_gridpoint = Constraint(model.odd_columns, model.rows, rule=at_most_one_metabolite_per_odd_column_gridpoint)
    model.set_reaction_column_position = Constraint(model.reactions, rule=set_reaction_column_position)
    model.set_metabolite_column_position = Constraint(model.metabolites, rule=set_metabolite_column_position)
    model.set_horizontal_slack_variable = Constraint(model.edges, rule=set_horizontal_slack_variable)
    model.calculate_horizontal_distance = Constraint(model.edges, rule=calculate_horizontal_distance)
    model.set_true_horizontal_distance1 = Constraint(model.edges, rule=set_true_horizontal_distance1)
    model.set_true_horizontal_distance2 = Constraint(model.edges, rule=set_true_horizontal_distance2)

    # vertical distance
    model.set_vertical_slack_variable = Constraint(model.edges, rule=set_vertical_slack_variable)
    model.calculate_vertical_distance = Constraint(model.edges, rule=calculate_vertical_distance)
    model.set_true_vertical_distance1 = Constraint(model.edges, rule=set_true_vertical_distance1)
    model.set_true_vertical_distance2 = Constraint(model.edges, rule=set_true_vertical_distance2)
    model.set_metabolite_row_position = Constraint(model.metabolites, rule=set_metabolite_row_position)
    model.set_reaction_row_position = Constraint(model.reactions, rule=set_reaction_row_position)


    def objective_rule(model):
        return sum(model.true_horizontal_distance[edge]*2 for edge in model.edges) + sum(model.true_vertical_distance[edge] for edge in model.edges)
    model.objective = Objective(rule=objective_rule, sense=minimize)

    # def temp_objective_rule(model):
    #     return sum(model.metabolites_at_grid_point[metabolite, column, row] for metabolite in model.metabolites for column in model.odd_columns for row in model.rows)
    # model.objective = Objective(rule=temp_objective_rule, sense=minimize)

    opt = SolverFactory('gurobi', tee=True)
    opt.solve(model, tee = True)
    for metabolite in model.metabolites:
        print(f"Metabolite {metabolite} is NOT active: {model.metabolite_not_placed[metabolite].value}")

    for reaction in model.reactions:
        for column in model.even_columns:
            for row in model.rows:
                if model.reactions_at_grid_point[reaction, column, row].value > 0.5:
                    print(f"Reaction {reaction} is at grid point {column}, {row}")
                    # print(f"Reaction {model.reaction_column_position[reaction].value}")
    for metabolite in model.metabolites:
        for column in model.odd_columns:
            for row in model.rows:
                if model.metabolites_at_grid_point[metabolite, column, row].value != 0:
                    print(f"Metabolite {metabolite} is at grid point {column}, {row}")
                    # print(f"Metabolite {model.metabolite_column_position[metabolite].value}")

    for edge in model.edges:
        print(f"Edge {edge} has distance {model.true_horizontal_distance[edge].value}")
    print(edges_)
    for edge in edges_:
        for reaction in reactions:
            if model.flow_from_reaction_through_edge[reaction, edge[0], edge[1]].value > 0.5:
                print(f"Edge {edge} for reaction {reaction} has flow {model.flow_from_reaction_through_edge[reaction,edge[0], edge[1]].value}")
def visualize_network(
        stoichiometric_matrix: np.ndarray,
        reaction_index_dict: {int: str},
        reaction_name_index_dict: {str: int},
        metabolite_index_dict: {int: str},
        metabolite_name_index_dict: {str: int},
        cobra_model: Model,
        to_remove_metabolites: [str] = list,
) -> None:

    (input_metabolites,
     output_metabolites,
     metabolites,
     input_reactions,
     output_reactions,
     reactions,
     edges) = create_nodes_and_edges_from_stoichiometric_matrix(stoichiometric_matrix,
                reaction_index_dict,
                reaction_name_index_dict,
                metabolite_index_dict,
                metabolite_name_index_dict,
                cobra_model)
    print("created nodes and edges")

    node_positions = formulate_shortest_distance_on_grid_MILP_2(edges,
                                                                input_reactions,
                                                                output_reactions,
                                                                metabolites,
                                                                reactions)

    return

    levels_dict = {value: [] for key, value in levels_dict_reverse.items()}
    for key, value in levels_dict_reverse.items():
        levels_dict[value].append(key)

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes
    G.add_nodes_from(reactions)
    G.add_nodes_from(metabolites)
    for node in G.nodes:
        for level, nodes in levels_dict.items():
            if node in nodes:
                G.nodes[node]['subset'] = level
                break
    # check for any missing nodes
    for node in G.nodes:
        try:
            pass
        except KeyError:
            print(node)

    # Add edges
    G.add_edges_from(edges)

    # Define node colors based on type
    color_map_reactions = []
    color_map_metabolites = []
    edge_color_map = []

    for node in G:
        if node in to_remove_metabolites:
            color_map_metabolites.append('red')
        elif node in input_reactions or node in input_metabolites:
            if node in input_reactions:
                color_map_reactions.append('green')
            elif node in input_metabolites:
                color_map_metabolites.append('green')
            else:
                print(node)
        elif node in output_reactions or node in output_metabolites:
            if node in output_reactions:
                color_map_reactions.append('lightgreen')
            elif node in output_metabolites:
                color_map_metabolites.append('lightgreen')
            else:
                print(node)
        else:
            if node in reactions:
                color_map_reactions.append('lightblue')
            elif node in metabolites:
                color_map_metabolites.append('cyan')
            else:
                print(node)

    color_map_edges = []
    alpha_map_edges = []
    not_removed_metabolites = [met for met in metabolites if met not in to_remove_metabolites]
    to_remove_edges = [edge for edge in edges if edge[0] in to_remove_metabolites or edge[1] in to_remove_metabolites]
    not_removed_edges = [edge for edge in edges if edge not in to_remove_edges]

    for edge in edges:
        if edge in to_remove_edges:
            color_map_edges.append('red')
            alpha_map_edges.append(0.5)
        else:
            color_map_edges.append('blue')
            alpha_map_edges.append(1)

    # figure settings
    fig, ax = plt.subplots(2)
    ax1, ax2 = ax
    fig.set_size_inches(10, 10)
    fig.suptitle('Network Graph of A and B Nodes')
    ax1.set_title('Full Network')
    ax2.set_title('Network with Removed Nodes')

    pos = nx.multipartite_layout(G)  # positions for all nodes
    nx.draw_networkx_nodes(G, pos, nodelist=reactions, node_size=700, node_color=color_map_reactions, ax=ax1)
    nx.draw_networkx_nodes(G, pos, nodelist=metabolites, node_size=700, node_color=color_map_metabolites, edgecolors=edge_color_map, ax=ax1)
    nx.draw_networkx_edges(G, pos, edgelist=edges,
                           edge_color=color_map_edges,
                           alpha=alpha_map_edges,
                           arrowstyle='->', arrowsize=12,
                           ax=ax1)

    nx.draw_networkx_labels(G, pos,
                            font_size=10,
                            verticalalignment='bottom',
                            font_family='sans-serif', ax=ax1)

    # with removed nodes
    G2 = nx.DiGraph()
    G2.add_nodes_from(reactions)
    G2.add_nodes_from(not_removed_metabolites)
    for node in G2.nodes:
        for level, nodes in levels_dict.items():
            if node in nodes:
                G2.nodes[node]['subset'] = level
                break
    # Add edges
    G.add_edges_from(not_removed_edges)
    not_removed_edge_color = ['blue' for edge in not_removed_edges]
    not_removed_edge_alpha = [1 for edge in not_removed_edges]
    not_removed_metabolite_color = [color_map_metabolites[idx] for idx, met in enumerate(metabolites) if met in not_removed_metabolites]
    pos2 = {key: item for key, item in pos.items() if key in G2.nodes}
    nx.draw_networkx_nodes(G2, pos2, nodelist=reactions, node_size=700, node_color=color_map_reactions)
    nx.draw_networkx_nodes(G2, pos2, nodelist=not_removed_metabolites, node_size=700, node_color=not_removed_metabolite_color, edgecolors=edge_color_map)
    nx.draw_networkx_edges(G2, pos2, edgelist=not_removed_edges,
                           edge_color=not_removed_edge_color,
                           alpha=not_removed_edge_alpha,
                           arrowstyle='->', arrowsize=12)
    nx.draw_networkx_labels(G2, pos2,
                            font_size=10,
                            verticalalignment='bottom',
                            font_family='sans-serif')

    plt.show()




### calculate levels algorithmically
def calculate_node_level(metabolite_nodes: [str],
                         input_reactions: [str],
                         output_reactions: [str],
                         reaction_nodes: [str],
                         edges_: [(str, str)]
                         ) -> {int: [str]}:
    # metabolite_nodes = metabolites
    # reaction_nodes = reactions
    # input_reactions = input_reactions
    # output_reactions = output_reactions
    # edges_ = edges

    # approach using direct assignment
    # TODO the metabolites are not assigned ott he correct exchange reactions (and probably others too)
    node_levels = {}
    remaining_nodes = metabolite_nodes + reaction_nodes
    run = True
    if run:
        node_levels = {}
        nodes_with_level = []
        node_edges_link_dict = {node: [edge for edge in edges_ if edge[0] == node] for node in remaining_nodes}
        for node in input_reactions:
            node_levels[node] = 1
            remaining_nodes.remove(node)
            nodes_with_level.append(node)

        node = nodes_with_level[0]
        connected_node_edges = node_edges_link_dict[node]
        connected_edge = connected_node_edges[0]

        for node in nodes_with_level:
            # print(f"node: {node}")
            for connected_edge in node_edges_link_dict[node]:
                # TODO loops make for problems in not finding a first metabolite
                # print(f"connected edge: {connected_edge}")
                connected_node = connected_edge[1]
                # print(f"connected node: {connected_node}")
                # check if connected node has any inputs that are not assigned yet
                edges_with_connected_node = [edge for edge in edges_ if edge[1] == connected_node]
                can_assign = True
                level_to_assign = node_levels[node] + 1
                if connected_node not in remaining_nodes:
                    can_assign = False
                for edge in edges_with_connected_node:
                    if not can_assign:
                        break
                    elif edge is connected_edge:
                        pass
                    elif edge[0] not in nodes_with_level:
                        # print(f"node not in nodes with level: {edge}")
                        can_assign = False
                        break
                    elif edge[0] in nodes_with_level:
                        if node_levels[edge[0]] > level_to_assign:
                            level_to_assign = node_levels[edge[0]]
                if can_assign:
                    # print(f"assigning level {level_to_assign} to {connected_node}")
                    node_levels[connected_node] = level_to_assign
                    remaining_nodes.remove(connected_node)
                    nodes_with_level.append(connected_node)

    assert len(remaining_nodes) == 0, f"Remaining nodes: {remaining_nodes}"
    return node_levels

def calculate_least_amount_of_edges_using_gurobi(
        model: Model, df_in: pd.DataFrame, task_number: int, sample_number: int
) -> pd.DataFrame:

    (filtered_stoichiometric_matrix, reaction_index_dict,
     reaction_name_index_dict, metabolite_index_dict,
     metabolite_name_index_dict) = identify_stoichiometric_matrix_of_sample_task_filtered_model(
        model, df_in, task_number, sample_number
    )
    # # for each row in the stoichiometric matrix, caulcaute the amount of non-zero columns
    # metabolite_number_of_reactions_dict = {
    #     idx: len(np.nonzero(row)[0]) for idx, row in enumerate(filtered_stoichiometric_matrix)
    # }
    # reaction_number_of_metabolites_dict = {
    #     idx: (np.sum(row > 0), np.sum(row < 0)) for idx, row in enumerate(filtered_stoichiometric_matrix.T)
    # }

    df = select_sample_and_task_from_dataframe(df_in, task_number, sample_number)
    reactions = df["reaction_name"].unique()
    reaction_indices = [i for i, r in enumerate(model.reactions) if r.id in reactions]
    exchange_reactions = [reaction for reaction in reactions if reaction.startswith("temporary_exchange_")]
    input_exchange_reactions = [reaction for reaction in exchange_reactions if (
        df[df["reaction_name"] == reaction]["reaction_formula"].str.startswith(" -->").iloc[0]
    )]
    output_exchange_reactions = [reaction for reaction in exchange_reactions if reaction not in input_exchange_reactions]
    non_exchange_reactions = [reaction for reaction in reactions if reaction not in exchange_reactions]

    # create edges from the stoichiometric matrix
    opt_model = ConcreteModel()
    # opt_model.output_exchange_reactions = pe.Set(initialize=output_exchange_reactions)
    # opt_model.input_exchange_reactions = pe.Set(initialize=input_exchange_reactions)
    # opt_model.main_reactions = pe.Set(initialize=non_exchange_reactions)
    # opt_model.metabolites = pe.Set(initialize=range(np.shape(filtered_stoichiometric_matrix)[0]))
    # opt_model.edges = pe.Set(initialize = [(i, j) for i in range(np.shape(filtered_stoichiometric_matrix)[0])
    #                                        for j in range(np.shape(filtered_stoichiometric_matrix)[1])
    #                                        if filtered_stoichiometric_matrix[i, j] != 0])
    #
    # undirected_symmetrical_edges = list(set([(i, j) for i in range(np.shape(filtered_stoichiometric_matrix)[0])
    #                                        for j in range(np.shape(filtered_stoichiometric_matrix)[1])
    #                                        if filtered_stoichiometric_matrix[i, j] != 0] +
    #                                 [(j, i) for i in range(np.shape(filtered_stoichiometric_matrix)[0])
    #                                          for j in range(np.shape(filtered_stoichiometric_matrix)[1])
    #                                             if filtered_stoichiometric_matrix[i, j] != 0]))
    #
    # opt_model.undirected_edges = pe.Set(initialize=undirected_symmetrical_edges)
    #
    # # decision variables
    # opt_model.active_metabolite = pe.Var(opt_model.metabolites, within=pe.Binary)
    # opt_model.active_edges = pe.Var(opt_model.edges, within=pe.Binary)

    # # objective function
    # def objective_rule(model):
    #     return pe.summation(model.active_edges)
    #
    # opt_model.active_edges.pprint()
    # # opt_model.objective = pe.Objective(minimize=pe.summation(opt_model.active_edges), sense=pe.minimize)
    # opt_model.objective = pe.Objective(rule=objective_rule, sense=pe.minimize)

    # constraints




    # nodes = metabolites and reactions
    # reactions are always on
    # metabolites can be turned on and off
    # xij = 1 if edge is on, 0 if edge is off
    # xij <= node_i
    # xij <= node_j
    # xij >= node_i + node_j - 1
    # inputreactions always have x_input_metaboltie = 1
    # outputreactions always have x_metabolite_output = 1
    # every reaction ri has at least one input metabolite and one output metabolite
    # xij >= 1 for ri
    # xji >= 1 for ri
    # yij = 1 if edge is on, 0 if edge is off

    # Create a model
    opt_model = ConcreteModel()
    undirected_symmetrical_edges = list(set([(i, j+np.shape(filtered_stoichiometric_matrix)[0]) for i in range(np.shape(filtered_stoichiometric_matrix)[0])
                                             for j in range(np.shape(filtered_stoichiometric_matrix)[1])
                                             if filtered_stoichiometric_matrix[i, j] != 0] +
                                            [(j + np.shape(filtered_stoichiometric_matrix)[0], i) for i in range(np.shape(filtered_stoichiometric_matrix)[0])
                                             for j in range(np.shape(filtered_stoichiometric_matrix)[1])
                                             if filtered_stoichiometric_matrix[i, j] != 0]))
    # Define sets for nodes and edges
    opt_model.Reactions = Set(initialize= range(np.shape(filtered_stoichiometric_matrix)[0]))
    opt_model.Metabolites = Set(initialize= range(np.shape(filtered_stoichiometric_matrix)[1]))
    opt_model.Nodes = Set(initialize=range(np.shape(filtered_stoichiometric_matrix)[0] + np.shape(filtered_stoichiometric_matrix)[1]))
    opt_model.Edges = Set(initialize=undirected_symmetrical_edges)
    # Define the set of special nodes
    output_nodes_start = np.shape(filtered_stoichiometric_matrix)[0] - len(exchange_reactions) + len(input_exchange_reactions)
    output_nodes_end = np.shape(filtered_stoichiometric_matrix)[0]
    opt_model.output_nodes = Set(initialize=range(output_nodes_start, output_nodes_end))

    input_nodes_start = np.shape(filtered_stoichiometric_matrix)[0] - len(exchange_reactions)
    input_nodes_end = np.shape(filtered_stoichiometric_matrix)[0] - len(output_exchange_reactions)
    opt_model.input_nodes = Set(initialize=range(input_nodes_start, input_nodes_end))

    # Define variables
    opt_model.EdgeOn = Var(opt_model.Edges, domain=Binary)
    opt_model.Flow = Var( opt_model.Edges, domain=NonNegativeReals)
    opt_model.NodeOn = Var(opt_model.Nodes, domain=Binary)
    for i in opt_model.Nodes:
        if i < np.shape(filtered_stoichiometric_matrix)[0]:  # If the node is in the first part
            opt_model.NodeOn[i].setlb(1)  # Set lower bound to 1
            opt_model.NodeOn[i].setub(1)  # Set upper bound to 1




    # Define constraints
    def edge_activation_rule_1(model, i, j):
        return model.EdgeOn[i, j] <= model.NodeOn[i]
    opt_model.EdgeActivation = Constraint(opt_model.Edges, rule=edge_activation_rule_1)

    def edge_activation_rule_2(model, i, j):
        return model.EdgeOn[i, j] <= model.NodeOn[j]
    opt_model.EdgeActivation = Constraint(opt_model.Edges, rule=edge_activation_rule_2)

    def edge_activation_rule_3(model, i, j):
        return model.EdgeOn[i, j] >= model.NodeOn[i] + model.NodeOn[j] - 1
    opt_model.EdgeActivation = Constraint(opt_model.Edges, rule=edge_activation_rule_3)

    def link_edge_with_flow_rule(model, i, j):
        return model.Flow[i, j] <= 1000 * model.EdgeOn[i, j]
    opt_model.LinkEdgeWithFlow = Constraint(opt_model.Edges, rule=link_edge_with_flow_rule)

    def link_node_with_edge_rule(model, i):
        return sum(model.EdgeOn[i, j] for j in model.Nodes if (i, j) in model.Edges) == sum(model.EdgeOn[j, i] for j in model.Nodes if (j, i) in model.Edges)

    opt_model.LinkNodeWithEdge = Constraint(opt_model.Nodes, rule=link_node_with_edge_rule)


    def flow_conservation_rule(model, node):
        if node in model.input_nodes:
            return sum(model.Flow[node, j] for j in model.Nodes if (node, j) in model.Edges) >= 1
        elif node in model.output_nodes:
            return sum(model.Flow[i, node] for i in model.Nodes if (i, node) in model.Edges) >= 1
        else:
            inflow = sum(model.Flow[i, j] for i, j in model.Edges if j == node and (i, j) in model.Edges)
            outflow = sum(model.Flow[i, j] for i, j in model.Edges if i == node and (i, j) in model.Edges)
            return inflow == outflow

    opt_model.FlowConservation = Constraint(opt_model.Nodes, rule=flow_conservation_rule)
    def flow_activation_rule(model, i, j):
        return model.Flow[i, j] <= model.NodeOn[i] * model.NodeOn[j] * 1000  # assuming a large capacity of 1000

    opt_model.FlowActivation = Constraint(opt_model.Edges, rule=flow_activation_rule)

    def flow_capacity_rule(model, i, j):
        return model.Flow[i, j] <= model.NodeOn[i] * 1000  # assuming a large capacity of 1000

    opt_model.FlowCapacity = Constraint(opt_model.Edges, rule=flow_capacity_rule)

    # Define objective function
    opt_model.obj = Objective(expr=sum(opt_model.Flow[i, j] for i, j in opt_model.Edges), sense=minimize)

    # Solve the model
    solver = SolverFactory('gurobi')
    solver.solve(opt_model)

    # Print the solution
    for i, j in opt_model.Edges:
        if opt_model.Flow[i, j].value > 0:
            print(f"Flow from {i} to {j}: {opt_model.Flow[i, j].value}")
    for i in opt_model.Nodes:
        if opt_model.NodeOn[i].value > 0:
            print(f"Node {i} is on")
def create_nodes_and_edges_from_stoichiometric_matrix(
        stoichiometric_matrix: np.ndarray,
        reaction_index_dict: {int: str},
        reaction_name_index_dict: {str: int},
        metabolite_index_dict: {int: str},
        metabolite_name_index_dict: {str: int},
        cobra_model: cobrapy_model
) -> (
        [str], # input_metabolites
        [str], # output_metabolites
        [str], # metabolites
        [str], # input_reactions
        [str], # output_reactions
        [str], # reactions
        [(str, str)] # edges
):
    input_metabolites = []
    output_metabolites = []
    metabolites = []
    input_reactions = []
    output_reactions = []
    reactions = []
    edges = []
    transposed_stoichiometric_matrix = stoichiometric_matrix.transpose()
    for idx, row in enumerate(transposed_stoichiometric_matrix):
        idx += 1
        reactions.append(f"{reaction_index_dict[idx-1]}")
        if np.shape(np.where(row == 0))[1] == len(row)-1:
            if np.shape(np.where(row == -1))[1] == 1:
                output_reactions.append(f"{reaction_index_dict[idx-1]}")
                # output_metabolites.append(f"output_{metabolite_name_index_dict[idx-1]}_{np.where(row == -1)[0][0] + 1}")
                output_metabolite_compartment = metabolite_index_dict[np.where(row == -1)[0][0]].split("[")[1].split("]")[0]
                output_metabolites.append(f"{metabolite_name_index_dict[np.where(row == -1)[0][0]]}[{output_metabolite_compartment}]")
            elif np.shape(np.where(row == 1))[1] == 1:
                input_reactions.append(f"{reaction_index_dict[idx-1]}")
                # input_metabolites.append(f"input_{metabolite_name_index_dict[idx-1]}_{np.where(row == 1)[0][0] + 1}")
                input_metabolite_compartment = metabolite_index_dict[np.where(row == 1)[0][0]].split("[")[1].split("]")[0]
                input_metabolites.append(f"{metabolite_name_index_dict[np.where(row == 1)[0][0]]}[{input_metabolite_compartment}]")
        for idx2, value in enumerate(row):
            idx2 += 1
            metabolite_name = metabolite_name_index_dict[idx2-1]
            metabolite_compartments = metabolite_index_dict[idx2-1].split("[")[1].split("]")[0]
            metabolite_full_name = f"{metabolite_name}[{metabolite_compartments}]"
            if value > 0:
                print(f"Adding edge: {reaction_index_dict[idx-1]} -> {metabolite_full_name}")
                edges.append((f"{reaction_index_dict[idx-1]}", f"{metabolite_full_name}"))
            elif value < 0:
                print(f"Adding edge: {metabolite_full_name} -> {reaction_index_dict[idx-1]}")
                edges.append((f"{metabolite_full_name}", f"{reaction_index_dict[idx-1]}"))

    metabolites = [f"{metabolite_name_index_dict[idx-1]}[{metabolite_index_dict[idx-1].split('[')[1].split(']')[0]}]" for idx in range(1, len(transposed_stoichiometric_matrix[0]) + 1)]

    print(f"Input metabolites: {input_metabolites}")
    return input_metabolites, output_metabolites, metabolites, input_reactions, output_reactions, reactions, edges

# # Example data
# input_nodes = ['input_1', 'input_2']
# output_nodes = ['output_1', 'output_2']
# exchange_nodes = input_nodes + output_nodes
# regular_nodes = ['A1', 'A2', 'A3', 'A4', 'A5']
# A_nodes = regular_nodes + exchange_nodes
# B_nodes = ['B1', 'B2', 'B3', 'B4', 'B5', "input_met_1", "input_met_2", "output_met_1", "output_met_2"]
# edges = [
#     ('A1', 'B1'), ('B1', 'A2'), ('A1', 'B2'), ('B2', 'A2'),
#     ('B2', 'A3'), ('A3', 'B3'), ('B3', 'A4'), ('A3', 'B4'),
#     ('B4', 'A4'), ('B3', 'A5'), ('A1', 'B5'), ('B5', 'A4'),
#     ('B5', 'A3'),("input_met_2", "A2"),  ("input_met_1", "A1"),
#     ("A4", "output_met_1"),  ("A5", "output_met_2")
# ]
#
# exchange_edges = [
#     ("input_1", "input_met_1"), ("input_2", "input_met_2"),
#     ("output_met_1", "output_1"), ("output_met_2", "output_2")
# ]
# edges_reverse = [(tup[1], tup[0]) for tup in edges]
# edges_reverse = []
# edges = edges + edges_reverse + exchange_edges



if __name__ == "__main__":

    input_metabolites = ["input_met_1", "input_met_2", "input_met_3"]
    output_metabolites = ["output_met_1", "output_met_2"]
    internal_metabolites = ["B1", "B2", "B3", "B4", "B5", "B6", "B8"]
    exchange_metabolites = input_metabolites + output_metabolites
    metabolites = internal_metabolites + exchange_metabolites
    to_remove_metabolites = ["B1", "B3", "B5"]
    not_removed_metabolites = [met for met in metabolites if met not in to_remove_metabolites]

    input_reactions = ["input_1", "input_2", "input_3"]
    output_reactions = ["output_1", "output_2"]
    internal_reactions = ["A1", "A2", "A3", "A4", "A5", "A6", "A7"]
    exchange_reactions = input_reactions + output_reactions
    reactions = internal_reactions + exchange_reactions

    input_edges = [("input_1", "input_met_1"), ("input_2", "input_met_2"), ("input_3", "input_met_3"),
                   ("input_met_1", "A1"), ("input_met_2", "A2"), ("input_met_3", "A3"), ("input_met_1", "A4")
                   ]
    output_edges = [("output_met_1", "output_1"), ("output_met_2", "output_2"),
                    ("A6", "output_met_2"), ("A7", "output_met_1"), ("A7", "output_met_2")]
    # for internal edges:
    # nodes to be rmeoved (and their corresponding edges) B1, B3, B5
    # A1 b1, A1 b2, A2 B2, A2 B3, A3 B3, A3 B6, B1 A4, B2 A4, B2 A5, B3 A5, A4 B4, A4 B5, A5 B4, A5 B5, B4 A6, B5 A7, B6 A7
    internal_edges = [
        ("A1", "B1"), ("A1", "B2"), ("A2", "B2"), ("A2", "B3"), ("A3", "B3"), ("A3", "B6"),
        ("B1", "A4"), ("B2", "A4"), ("B2", "A5"), ("B3", "A5"), ("A4", "B4"), ("A4", "B5"),
        ("A5", "B4"), ("A5", "B5"), ("B4", "A6"), ("B5", "A7"), ("B6", "A7"),
        ("A1", "B8"), ("B8", "A2"),
    ]
    edges = internal_edges + input_edges + output_edges
    to_remove_edges = [edge for edge in edges if edge[0] in to_remove_metabolites or edge[1] in to_remove_metabolites]
    not_removed_edges = [edge for edge in edges if edge not in to_remove_edges]
    print(edges)
    ### calculate reachable reactions for all reactions recursively
    # reachable_reactions = calculate_reachable_reactions(edges, reactions)

    def calculate_reachable_reactions(edges: [(str, str)], reactions: [str]) -> {str: [str]}:
        # Create a set of all nodes
        nodes = list(set([edge[0] for edge in edges] + [edge[1] for edge in edges]))
        reachable_nodes = {node: set() for node in nodes}
        start_edges_dict = {node: [edge for edge in edges if edge[0] == node] for node in nodes}

        for reaction in reactions:
            queue = deque([reaction])
            visited = set([reaction])

            while queue:
                current_node = queue.popleft()
                for neighbor in start_edges_dict[current_node]:
                    if neighbor not in reachable_nodes[reaction]:
                        reachable_nodes[reaction].add(neighbor)

                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

                    for further_reach in reachable_nodes[neighbor]:
                        if further_reach not in reachable_nodes[reaction]:
                            reachable_nodes[reaction].add(further_reach)

        reachable_reactions = {
            reaction: [rxn for rxn in reachable_nodes[reaction] if rxn in reactions]
            for reaction in reactions
        }

        return reachable_reactions


    from collections import deque


    def calculate_reachable_reactions(edges: [(str, str)], reactions: [str]) -> {str: [str]}:
        # Create a set of all nodes
        nodes = list(set([edge[0] for edge in edges] + [edge[1] for edge in edges]))
        reachable_nodes = {node: set() for node in nodes}
        start_edges_dict = {node: [] for node in nodes}
        for edge in edges:
            start_edges_dict[edge[0]].append(edge[1])
        fully_explored_nodes = set()
        def recursive_search(node, visited):
            visited.add(node)
            for neighbor in start_edges_dict[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    recursive_search(neighbor, visited)
            return visited

        for reaction in reactions:
            if reaction in fully_explored_nodes:
                print(f"Reaction {reaction} was already fully explored", reachable_nodes[reaction])
                reachable_nodes[reaction].update(reachable_nodes[reaction])
                continue

            queue = deque([reaction])
            visited = set([reaction])
            while queue:
                current_node = queue.popleft()
                for neighbor in start_edges_dict[current_node]:
                    if neighbor in fully_explored_nodes:
                        reachable_nodes[reaction].update(reachable_nodes[neighbor])
                    else:
                        if neighbor not in reachable_nodes[reaction]:
                            print(f'Adding {neighbor} to reachable nodes of {reaction}')
                            reachable_nodes[reaction].add(neighbor)
                            queue.append(neighbor)
                            visited.add(neighbor)

                fully_explored_nodes.add(current_node)

        # Filter out non-reaction nodes from the final result
        reachable_reactions = {
            reaction: [rxn for rxn in reachable_nodes[reaction] if rxn in reactions]
            for reaction in reactions
        }

        return reachable_reactions


    def calculate_reachable_reactions(edges: [(str, str)], reactions: [str]) -> {str: [str]}:
        nodes = list(set([edge[0] for edge in edges] + [edge[1] for edge in edges]))
        reachable_nodes = {node: set() for node in nodes}
        start_edges_dict = {node: [] for node in nodes}
        for edge in edges:
            start_edges_dict[edge[0]].append(edge[1])
        fully_explored_nodes = {}
        def dfs(current_node):
            if current_node in fully_explored_nodes:
                return fully_explored_nodes[current_node]
            local_reachable = set()

            for neighbor in start_edges_dict.get(current_node, []):
                if neighbor not in local_reachable:
                    local_reachable.add(neighbor)
                    local_reachable.update(dfs(neighbor))

            fully_explored_nodes[current_node] = local_reachable
            return local_reachable

        for reaction in reactions:
            reachable_nodes[reaction].update(dfs(reaction))

        reachable_reactions = {
            reaction: [rxn for rxn in reachable_nodes[reaction] if rxn in reactions]
            for reaction in reactions
        }

        return reachable_reactions

    reachable_reactions = calculate_reachable_reactions(edges, reactions)
    print(reachable_reactions)
    # print(edges)
    # # simplified model for testing 1 input, 1 output, A1 > B1, A1 > B2, B2 > A2, A2 > B3, B1 > out, B3 > out
    # input_metabolites = ["input_met_1"]
    # output_metabolites = ["output_met_1"]
    # internal_metabolites = ["B1", "B2", "B3", "B4"]
    # internal_metabolites = ["B1", "B2"]
    # exchange_metabolites = input_metabolites + output_metabolites
    # metabolites = internal_metabolites + exchange_metabolites
    #
    # input_reactions = ["input_1"]
    # output_reactions = ["output_1"]
    # internal_reactions = ["A1", "A2", "A3", "A4"]
    # internal_reactions = ["A1", "A2"]
    # exchange_reactions = input_reactions + output_reactions
    # reactions = internal_reactions + exchange_reactions
    #
    # input_edges = []
    # output_edges = []
    # # simplified model for testing 1 input, 1 output, A1 > B1, A1 > B2, B2 > A2, A2 > B3, B1 > out, B3 > out
    # internal_edges = [
    #     ("A1", "B1"), ("B1","A2"),
    #     ("input_1", "input_met_1"), ("input_met_1", "A1"), ("input_met_1", "A3"),
    #     ("A1", "output_met_1"), ("A2", "output_met_1"), ("A4", "output_met_1"), ("output_met_1", "output_1"),
    #     ("A3", "B2"), ("A3", "B3"), ("B2", "A4"),
    #     ("B3", "A2"), ("input_1", "output_met_1"),
    #     ("A1", "B4"), ("B4", "A2")
    # ]
    # internal_edges = [
    #     ("input_1", "input_met_1"), ("input_met_1", "A1"),
    #     ("A1", "B1"), ("B1", "A2"),
    #     ("A1", "B2"), ("B2", "A2"),
    #     ("A2", "output_met_1"), ("output_met_1", "output_1"),
    # ]
    # edges = internal_edges + input_edges + output_edges
    # print(edges)
    #
    # node_positions = formulate_shortest_distance_on_grid_MILP_2(edges,
    #                                                             input_reactions,
    #                                                             output_reactions,
    #                                                             metabolites,
    #                                                             reactions,
    #                                                             reachable_reactions,
    #                                                             )

def create_edges_nodes_from_model(
        model: cobrapy_model,
) -> (
        [str], # edges
        [str], # nodes
):
    stoichiometric_matrix = create_stoichiometric_matrix(model)
    reaction_index_dict = {idx: reaction.id for idx, reaction in enumerate(model.reactions)}
    reaction_name_index_dict = {reaction.id: idx for idx, reaction in enumerate(model.reactions)}
    metabolite_index_dict = {idx: metabolite.id for idx, metabolite in enumerate(model.metabolites)}
    metabolite_name_index_dict = {metabolite.id: idx for idx, metabolite in enumerate(model.metabolites)}

    nodes = [metabolite.id for metabolite in model.metabolites] + [reaction.id for reaction in model.reactions]
    edges = []
    print(len(metabolite_index_dict), len(reaction_index_dict))
    transposed_stoichiometric_matrix = stoichiometric_matrix.transpose()
    for idx, row in enumerate(transposed_stoichiometric_matrix):
        reaction = reaction_index_dict[idx]
        for idx2, value in enumerate(row):
            metabolite = metabolite_index_dict[idx2]
            if value > 0:
                edges.append((reaction, metabolite))
            elif value < 0:
                edges.append((metabolite, reaction))

    return edges, metabolites

if __name__ == "__main__":

    # main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_no_threshold_fractional_options_test_enum"
    main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Base_files"
    cellfie_model_location = main_folder + "\\modelCellfie.json"
    cobra_model = load_json_model(cellfie_model_location)
    edges, nodes = create_edges_nodes_from_model(cobra_model)
    print('loaded edges and nodes')

    nodes_set = set([node for edge in edges for node in edge])
    nodes_df = pd.DataFrame({'id': list(nodes_set)})

    nodes_df['name'] = nodes_df['id']
    edges_df = pd.DataFrame(edges, columns=['source', 'target'])
    edges_df['interaction'] = 'pp'
    p4c.cytoscape_ping()
    p4c.create_network_from_data_frames(nodes_df, edges_df, title='My Network', collection='My Collection')
    # result = p4c.analyze_network()
    # p4c.save_session(main_folder+ '\\full_metabolic_network.cys')
    # create_nodes_and_edges_from_stoichiometric_matrix
