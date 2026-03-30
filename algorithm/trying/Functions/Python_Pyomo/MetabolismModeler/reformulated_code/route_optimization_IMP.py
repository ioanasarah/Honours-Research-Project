"""
This module contains the classes and functions used to create the route optimization models using Pyomo.
It takes the files in the combinations folder created in preprocessing or call the preprocessing functions to create
the necessary files. A dictionary with settings is used to create the specific files and models.

 - A "combinations" directory based on the expression data, sample information, and tasklist that contains:
        - A closed base model (with the same name as the model) in the "models" directory
        - A irreversible version of the base model in the "models" directory
        - Expression data in the form of a json file parsed to the base models:
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



                - GPR_method:       minSum (taking the sum of OR values),
                                    minMax (taking the maximum of OR values),
                                    minMean (taking the mean of OR values),
                                    minMedian (taking the median of OR values),
                                    trimmeanSum_V (taking the mean of AND values after trimming V%),
                                                and then taking the sum of OR values,
                                    trimmeanMax_V (taking the mean of AND values after trimming V%),
                                                and then taking the maximum of OR values
                                    trimminSum_V (taking the min of AND values after trimming V%),
                                                and then taking the sum of OR values
                                    trimminMax_V (taking the min of AND values after trimming V%),
                                                and then taking the maximum of OR values

                - promiscuity_method:   Y (yes, a gene being the determining factor for multiple reactions
                                                is reduced by 1/n)
                                        N (no)

                - model_type:       rev (reversible),
                                    irr (irreversible)
            - A list  linking sample IDs to corresponding numbers as .json
            - A file with information on the added reactions and bounds within the model as .json
            - A metadata/log file with information on the creation of the combination files and whenever anything has been run
            - (Optional only for reversible) A rev_to_irrev file with information on the added reactions and bounds
                    within the model as .json for the irreversible model
            - (Optional only for irreversible) A list linking each reaction's index with its
                    corresponding irreversible reaction index as .json

Initially created on 2024-10-16 by Jelle Bonthuis (MaCSBio)
Updated on {date} by {author}:
    - {list of changes}
"""
# noinspection PyUnresolvedReferences
from pyomo.environ import (ConcreteModel, Var, Objective,
                           Constraint, minimize, SolverFactory,
                           quicksum, Binary,
                           NonNegativeReals, Reals, RangeSet, value)
import os
import re
import time
import warnings
from json import load, dump
import pandas as pd
from cobra import Model, Reaction
from cobra.io import load_json_model, save_json_model, save_matlab_model
from cobra.util import create_stoichiometric_matrix
from pyomo.core.base.param import Param
from datetime import datetime

version = "1.0.0"


class ModelCreator:
    """
    This class contains the functions to create the model files and the model objects for the route optimization.
    It uses the settings dictionary to create the model files and objects.
    """

    def __init__(self,
                 settings: dict,
                 main_data_folder: str,
                 inputs_dict: (dict, None) = None
                 ):
        """
        The constructor for the ModelCreator class.
        :param settings: A dictionary with settings for the model creation.
        """
        # add settings
        self.settings = settings
        self.main_data_folder = main_data_folder
        self.load_settings(settings)

        if inputs_dict is None:
            self.load_data()

        self.cobra_stoichiometric_matrix = create_stoichiometric_matrix(self.cobra_model)
        self.initialize_model()

    # noinspection PyAttributeOutsideInit
    def load_settings(self, settings: dict):
        self.model_name = settings["__global_options"].get("model_version", None)
        self.model_type = settings["__global_options"].get("model_type", False)
        self.expression_dataset_name = settings["__global_options"].get("expression_dataset", None)
        self.tasklist_name = settings["__global_options"].get("tasklist", None)
        self.transformation = settings["__global_options"].get("transformation", None)
        self.thresholding = settings["__global_options"].get("thresholding", None)
        self.GPR_method = settings["__global_options"].get("GPR_method", None)
        self.promiscuity_method = settings["__global_options"].get("promiscuity_method", None)
        self.epsilon = settings["__global_options"].get("epsilon", None)
        self.expressionless_reaction_value = settings["__global_options"].get("expressionless_reaction_value", None)
        self.create_log_file = settings["__global_options"].get("create_log_file", None)
        self.solver_settings = settings["__global_options"].get("solver_settings", None)
        self.tasks_to_run = settings["__global_options"].get("tasks", None)
        self.samples_to_run = settings["__global_options"].get("samples", None)
        self.unbounded_lower_bound_value = settings["__global_options"].get("unbounded_lower_bound_value", None)
        self.unbounded_upper_bound_value = settings["__global_options"].get("unbounded_upper_bound_value", None)
        self.create_log_file = settings["__global_options"].get("create_log_file", False)
        self.name_of_run = settings["__global_options"].get("name_of_run", "no_run_name_set")
        self.create_log_file = settings["__global_options"].get("create_log_file", True)
        self.combinations_folder = f"{self.main_data_folder}/combinations/{self.model_name}_{self.expression_dataset_name}_{self.tasklist_name}"
        self.result_output_folder = settings["__global_options"].get("result_output_folder", f"{self.combinations_folder}\\Results")
        self.overwrite_run_files = settings["__global_options"].get("overwrite_run_files", False)
        self.skip_existing_files = settings["__global_options"].get("skip_existing_files", False)
        self.ignore_sample_warnings = settings["__global_options"].get("ignore_sample_warnings", False)
        self.minimize_filtered_model = settings["__global_options"].get("minimize_filtered_model", False)
        if self.overwrite_run_files and self.skip_existing_files:
            warnings.warn("Both overwrite_run_files and skip_existing_files are set to True, which will lead to only skipping existing files")

        for option in [self.model_name, self.model_type, self.expression_dataset_name, self.tasklist_name,
                       self.transformation, self.thresholding, self.GPR_method, self.promiscuity_method,
                       self.epsilon, self.expressionless_reaction_value, self.create_log_file, self.solver_settings,
                       self.tasks_to_run, self.samples_to_run, self.unbounded_lower_bound_value, self.unbounded_upper_bound_value]:
            if option is None or option is False:
                warnings.warn(f"Option {option} is not set in the settings dictionary")

        # create names for the files and locations
        self.expression_name = (f"{self.transformation}__{self.thresholding}__"
                                f"{self.GPR_method}__{self.promiscuity_method}__"
                                f"{self.model_type}__reaction_expression_dataset.feather")
        self.per_reaction_adjusted_gene_contribution_name = (
            f"{self.transformation}__{self.thresholding}__"
            f"{self.GPR_method}__{self.promiscuity_method}__"
            f"rev__per_reaction_adjusted_gene_contribution.parquet"
        )
        model_irreversibility = "irreversible" if self.model_type == "irr" else "reversible"
        self.cobra_model_name = f"model_{self.model_name}__{model_irreversibility}.json"
        self.task_structure_name = f"task_structure_{self.tasklist_name}.json"
        self.rev_to_irrev_name = f"{self.model_name}__rev_to_irrev_dict.json"
    # noinspection PyAttributeOutsideInit
    def load_data(self):
        """
        This function loads the data needed to create the model files and objects.
        """
        self.cobra_model = load_json_model(f"{self.combinations_folder}/{self.cobra_model_name}")
        self.expression_data = pd.read_feather(f"{self.combinations_folder}/{self.expression_name}")
        if self.expression_data.isnull().values.any():
            missing_data_value_string = (f"There are missing values in the expression data that "
                                         f"have been replaced with the expressionless_reaction_value: {self.expressionless_reaction_value}")
            warnings.warn(missing_data_value_string)
            self.expression_data = self.expression_data.replace(-1, self.expressionless_reaction_value)

        with open(f"{self.combinations_folder}/{self.task_structure_name}", "r") as file:
            self.task_structure = load(file)
        with open(f"{self.combinations_folder}/{self.rev_to_irrev_name}", "r") as file:
            self.rev_to_irrev = load(file)

        self.per_reaction_adjusted_gene_contribution = pd.read_parquet(
            f"{self.combinations_folder}/{self.per_reaction_adjusted_gene_contribution_name}"
        )
    # noinspection PyAttributeOutsideInit
    def initialize_model(self):
        print("Initializing model")
        self.current_task = None
        self.current_sample = None
        self.route_model = ConcreteModel()

        self.route_model.constant_reactions = RangeSet(0, len(self.cobra_model.reactions) - 1)
        self.route_model.metabolites = RangeSet(0, len(self.cobra_model.metabolites) - 1)
        rxn = self.cobra_model.metabolites[8308]
        index_met = self.cobra_model.metabolites.index(rxn)
        lbs = [reaction.lower_bound for reaction in self.cobra_model.reactions]
        ubs = [reaction.upper_bound for reaction in self.cobra_model.reactions]
        # noinspection PyTypeChecker
        reaction_lb_dict = {reaction: lb for reaction, lb in zip(self.route_model.constant_reactions, lbs)}
        # noinspection PyTypeChecker
        reaction_ub_dict = {reaction: ub for reaction, ub in zip(self.route_model.constant_reactions, ubs)}

        def bound_rule(model, rxn):
            return (reaction_lb_dict[rxn], reaction_ub_dict[rxn])

        if self.model_type == "irr":
            self.route_model.constant_reaction_fluxes = Var(self.route_model.constant_reactions,
                                                            bounds=bound_rule,
                                                            domain=NonNegativeReals)
        elif self.model_type == "rev":
            self.route_model.constant_reaction_fluxes = Var(self.route_model.constant_reactions,
                                                            bounds=bound_rule,
                                                            domain=Reals)
            raise NotImplementedError("Reversible models are not implemented yet, specifically requires additional constraints for minus flux values")
        # noinspection PyTypeChecker
        empty_expression_dict = {reaction: 1 for reaction in self.route_model.constant_reactions}
        self.route_model.reaction_expression = Param(self.route_model.constant_reactions, initialize=empty_expression_dict, mutable=True,
                                                     domain=NonNegativeReals)
        start_time = time.time()
        print("Adding constraints")
        stoichiometric_dict = {}
        for rxn_idx, reaction in enumerate(self.cobra_model.reactions):
            for met in reaction.metabolites:
                met_idx = self.cobra_model.metabolites.index(met)
                if met_idx not in stoichiometric_dict:
                    stoichiometric_dict[met_idx] = {}
                stoichiometric_dict[met_idx][rxn_idx] = self.cobra_stoichiometric_matrix[met_idx, rxn_idx]

        print("Time for creating stoichiometric_dict:", time.time() - start_time)
        start = time.time()

        for idx in self.route_model.metabolites:
            if idx not in stoichiometric_dict:
                stoichiometric_dict[idx] = {}
        self.stoichiometric_dict = stoichiometric_dict

        def sv_zero_rule(model, met):
            if stoichiometric_dict[met] != {}:
                return (
                        sum(stoichiometric_dict[met][rxn] * model.constant_reaction_fluxes[rxn]
                            for rxn in stoichiometric_dict[met]) == 0
                )
            else:
                return Constraint.Skip

        self.route_model.SV_constraint = Constraint(self.route_model.metabolites, rule=sv_zero_rule)
        self.route_model.constant_reaction_is_active = Var(self.route_model.constant_reactions, domain=Binary)

        if self.model_type == "irr":
            def reaction_active_rule(model, rxn):
                return model.constant_reaction_fluxes[rxn] - (model.constant_reaction_is_active[rxn] * self.epsilon) >= 0

            def reaction_inactive_rule(model, rxn):
                return model.constant_reaction_fluxes[rxn] - (model.constant_reaction_is_active[rxn] * reaction_ub_dict[rxn]) <= 0

            def reaction_irreversible_reaction_pair_rule(model, rxn):
                return (
                    model.constant_reaction_is_active[rxn] + model.constant_reaction_is_active[self.rev_to_irrev[rxn][1] - 1] <= 1
                    if rxn < len(self.rev_to_irrev) and len(self.rev_to_irrev[rxn]) > 1
                    else Constraint.Skip
                )

            self.route_model.reaction_active_constraint = Constraint(self.route_model.constant_reactions, rule=reaction_active_rule)
            self.route_model.reaction_inactive_constraint = Constraint(self.route_model.constant_reactions, rule=reaction_inactive_rule)
            self.route_model.reaction_irreversible_reaction_pair_constraint = Constraint(self.route_model.constant_reactions,
                                                                                         rule=reaction_irreversible_reaction_pair_rule)
        elif self.model_type == "rev":
            raise NotImplementedError("Reversible models are not implemented yet, specifically requires additional constraints for minus flux values")

        # noinspection PyTypeChecker,  PyTestUnpassedFixture, PyTestUnresolvedReferences, PyPropertyDefinition
        self.route_model.objective = Objective(expr=quicksum(
            self.route_model.constant_reaction_is_active[rxn] / self.route_model.reaction_expression[rxn].value for rxn in self.route_model.constant_reactions),
                                               sense=minimize)
        print("Time for adding constraints:", time.time() - start)

    def add_task_specific_reactions(
            self,
            task: int
    ) -> (list, list):
        task_id = str(task)
        reactions_to_add = []
        input_metabolites = self.task_structure[task_id]["input_metabolites"]
        output_metabolites = self.task_structure[task_id]["output_metabolites"]
        input_lower_bound = self.task_structure[task_id]["input_metabolite_lower_bounds"]
        input_upper_bound = self.task_structure[task_id]["input_metabolite_upper_bounds"]
        output_lower_bound = self.task_structure[task_id]["output_metabolite_lower_bounds"]
        output_upper_bound = self.task_structure[task_id]["output_metabolite_upper_bounds"]

        def replace_unbounded(bound_list, replacement_value):
            return [bound if bound != "unbounded" else replacement_value for bound in bound_list]

        input_lower_bound, output_lower_bound = [
            replace_unbounded(bound_list, self.unbounded_lower_bound_value)
            for bound_list in [input_lower_bound, output_lower_bound]
        ]

        input_upper_bound, output_upper_bound = [
            replace_unbounded(bound_list, self.unbounded_upper_bound_value)
            for bound_list in [input_upper_bound, output_upper_bound]
        ]

        output_metabolite_dict = {metabolite: {
            "lower_bound": bound[0],
            "upper_bound": bound[1]
        } for metabolite, bound in zip(output_metabolites, zip(output_lower_bound, output_upper_bound))}

        input_metabolite_dict = {metabolite: {
            "lower_bound": bound[0],
            "upper_bound": bound[1]
        } for metabolite, bound in zip(input_metabolites, zip(input_lower_bound, input_upper_bound))}

        metabolites_in_both = list(set(input_metabolites).intersection(output_metabolites))
        if len(metabolites_in_both) > 0:
            for metabolite in metabolites_in_both:
                input_metabolites.remove(metabolite)
                output_metabolites.remove(metabolite)
                if self.model_type == "irr":
                    if input_metabolite_dict[metabolite]["lower_bound"] == output_metabolite_dict[metabolite]["lower_bound"]:
                        lower_bound_input = 0
                        lower_bound_output = 0
                    elif input_metabolite_dict[metabolite]["lower_bound"] - output_metabolite_dict[metabolite]["lower_bound"] > 0:
                        lower_bound_input = (input_metabolite_dict[metabolite]["lower_bound"] -
                                             output_metabolite_dict[metabolite]["lower_bound"]
                                             )
                        lower_bound_output = 0
                    else:
                        lower_bound_input = 0
                        lower_bound_output = (output_metabolite_dict[metabolite]["lower_bound"] -
                                              input_metabolite_dict[metabolite]["lower_bound"]
                                              )

                    new_reaction = Reaction(f"temporary_exchange_{metabolite}_f")
                    new_reaction.name = f"temporary_exchange_{metabolite}_f"
                    new_reaction.subsystem = "temporary_exchange"
                    new_reaction.lower_bound = lower_bound_input
                    new_reaction.upper_bound = input_metabolite_dict[metabolite]["upper_bound"]
                    metabolite_object = self.cobra_model.metabolites.get_by_id(metabolite)
                    new_reaction.add_metabolites({metabolite_object: 1})
                    self.cobra_model.add_reactions([new_reaction])
                    reactions_to_add.append(new_reaction)

                    new_reaction = Reaction(f"temporary_exchange_{metabolite}_r")
                    new_reaction.name = f"temporary_exchange_{metabolite}_r"
                    new_reaction.subsystem = "temporary_exchange"
                    new_reaction.lower_bound = lower_bound_output
                    new_reaction.upper_bound = output_metabolite_dict[metabolite]["upper_bound"]
                    metabolite_object = self.cobra_model.metabolites.get_by_id(metabolite)
                    new_reaction.add_metabolites({metabolite_object: -1})
                    self.cobra_model.add_reactions([new_reaction])
                    reactions_to_add.append(new_reaction)
                elif self.model_type == "rev":
                    new_reaction = Reaction(f"temporary_exchange_{metabolite}")
                    new_reaction.name = f"temporary_exchange_{metabolite}"
                    new_reaction.subsystem = "temporary_exchange"
                    new_reaction.lower_bound = - output_metabolite_dict[metabolite]["upper_bound"]
                    new_reaction.upper_bound = input_metabolite_dict[metabolite]["upper_bound"]
                    metabolite_object = self.cobra_model.metabolites.get_by_id(metabolite)
                    new_reaction.add_metabolites({metabolite_object: 1})
                    self.cobra_model.add_reactions([new_reaction])
                    reactions_to_add.append(new_reaction)

        for metabolite in input_metabolites:
            new_reaction = Reaction(f"temporary_exchange_{metabolite}")
            new_reaction.name = f"temporary_exchange_{metabolite}"
            new_reaction.subsystem = "temporary_exchange"
            new_reaction.lower_bound = input_metabolite_dict[metabolite]["lower_bound"]
            new_reaction.upper_bound = input_metabolite_dict[metabolite]["upper_bound"]
            metabolite_object = self.cobra_model.metabolites.get_by_id(metabolite)
            new_reaction.add_metabolites({metabolite_object: 1})
            self.cobra_model.add_reactions([new_reaction])
            reactions_to_add.append(new_reaction)
        for metabolite in output_metabolites:
            new_reaction = Reaction(f"temporary_exchange_{metabolite}")
            new_reaction.name = f"temporary_exchange_{metabolite}"
            new_reaction.subsystem = "temporary_exchange"
            new_reaction.lower_bound = output_metabolite_dict[metabolite]["lower_bound"]
            new_reaction.upper_bound = output_metabolite_dict[metabolite]["upper_bound"]
            metabolite_object = self.cobra_model.metabolites.get_by_id(metabolite)
            new_reaction.add_metabolites({metabolite_object: -1})
            self.cobra_model.add_reactions([new_reaction])
            reactions_to_add.append(new_reaction)

        SV_indices = [self.cobra_model.metabolites.index(metabolite) for metabolite in input_metabolites + output_metabolites + metabolites_in_both]

        return reactions_to_add, SV_indices

    def create_task_specific_models(self, task: int):
        try:
            for met_index in self.route_model.metabolites:
                if met_index in self.route_model.SV_constraint: # todo added for problem with missing
                    self.route_model.SV_constraint[met_index].activate()
            self.route_model.del_component(self.route_model.temporary_reactions)
            self.route_model.del_component(self.route_model.temporary_exchange_SV_constraint)
            self.route_model.del_component(self.route_model.temporary_exchange_fluxes)
        except Exception as e:
            print(e, ", this is expected in the first run")

        cobra_reactions_to_remove = [reaction.id for reaction in self.cobra_model.reactions if reaction.id.startswith("temporary_exchange_")]
        if len(cobra_reactions_to_remove) > 0:
            self.cobra_model.remove_reactions(cobra_reactions_to_remove)

        (reactions_to_add, SV_indices) = self.add_task_specific_reactions(task)

        temp_stoichiometric_dict = {}
        for reaction in reactions_to_add:
            rxn_index = self.cobra_model.reactions.index(reaction)
            for metabolite in reaction.metabolites:
                met_index = self.cobra_model.metabolites.index(metabolite)
                if met_index not in temp_stoichiometric_dict:
                    temp_stoichiometric_dict[met_index] = {}
                temp_stoichiometric_dict[met_index][rxn_index] = reaction.metabolites[metabolite]

        self.temp_stoichiometric_dict = temp_stoichiometric_dict
        self.route_model.temporary_reactions = RangeSet(len(self.cobra_model.reactions) - len(reactions_to_add), len(self.cobra_model.reactions) - 1)
        lbs = [reaction.lower_bound for reaction in reactions_to_add]
        # noinspection PyTypeChecker
        reaction_lb_dict = {reaction: lb for reaction, lb in zip(self.route_model.temporary_reactions, lbs)}
        ubs = [reaction.upper_bound for reaction in reactions_to_add]
        # noinspection PyTypeChecker
        reaction_ub_dict = {reaction: ub for reaction, ub in zip(self.route_model.temporary_reactions, ubs)}

        def bound_rule(model, rxn):
            return reaction_lb_dict[rxn], reaction_ub_dict[rxn]

        if self.model_type == "irr":
            self.route_model.temporary_exchange_fluxes = Var(self.route_model.temporary_reactions,
                                                             bounds=bound_rule,
                                                             domain=NonNegativeReals)
        elif self.model_type == "rev":
            self.route_model.temporary_exchange_fluxes = Var(self.route_model.temporary_reactions,
                                                             bounds=bound_rule,
                                                             domain=Reals)

        non_zero_columns_per_index = {
            met_index: [rxn for rxn in self.stoichiometric_dict[met_index] if self.stoichiometric_dict[met_index][rxn] != 0]
            for met_index in SV_indices
        }
        for met_index in SV_indices:
            self.route_model.SV_constraint[met_index].deactivate()

        #     self.persistent_gurobi.remove_constraint(self.route_model.SV_constraint[met_index])

        def temporary_exchange_SV_rule(model, met):
            return (
                    sum(self.stoichiometric_dict[met][rxn] * model.constant_reaction_fluxes[rxn] for rxn in non_zero_columns_per_index[met])
                    + sum(temp_stoichiometric_dict[met][rxn] * model.temporary_exchange_fluxes[rxn] for rxn in temp_stoichiometric_dict[met])
                    == 0
            )

        self.route_model.temporary_exchange_SV_constraint = Constraint(SV_indices, rule=temporary_exchange_SV_rule)

    def set_sample_specific_expression(self, sample: int):
        self.expression_data.fillna(self.expressionless_reaction_value, inplace=True)
        self.expression_data = self.expression_data.replace(-1, self.expressionless_reaction_value)
        for rxn in self.route_model.constant_reactions:
            self.route_model.reaction_expression[rxn] = self.expression_data[f"sample_{sample}"].iloc[rxn]

    def check_if_task_exists(self, task: int) -> bool:
        return str(task) in self.task_structure

    def check_if_sample_exists(self, sample: int) -> bool:
        sample_id = f"sample_{sample}"
        return sample_id in self.expression_data.columns

    def check_if_task_is_viable(self) -> bool:
        solution_ = self.cobra_model.optimize()
        if solution_.status != "optimal":
            return False
        return True

    def set_log_file_name(self, task_number: int, sample_number: int):
        if self.create_log_file:
            log_file_iterator = 1
            log_file_name = f"{self.result_output_folder}\\{self.name_of_run}\\Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}_gurobi"
            log_file_name_to_use = f'{log_file_name}.log'
            if not self.overwrite_run_files:
                while True:
                    if os.path.exists(log_file_name):
                        log_file_name_to_use = f"{log_file_name}_Try{log_file_iterator}.log"
                        log_file_iterator += 1
                    else:
                        break
            self.persistent_gurobi.options["LogFile"] = log_file_name_to_use

    def check_model_viability(self) -> bool:
        if self.result.solver.termination_condition == "infeasible":
            warnings.warn(f"Model for task {self.current_task} and sample {self.current_sample} is not optimal")
            return False
        if value(self.route_model.objective) < 0:
            warnings.warn(f"Model for task {self.current_task} and sample {self.current_sample} has a negative objective value")
            return False
        return True

    def check_if_optimization_model_respects_SV_constraints(self):
        reactions_with_positive_flux = [
            idx for idx, rxn in enumerate(self.cobra_model.reactions) if
                not rxn.id.startswith("temporary_") and
                abs(self.route_model.constant_reaction_fluxes[idx].value) > self.epsilon
        ]
        reactions_with_positive_flux.extend(
            idx for idx, rxn in enumerate(self.cobra_model.reactions) if
                rxn.id.startswith("temporary_") and
                abs(self.route_model.temporary_exchange_fluxes[idx].value) > self.epsilon
        )
        metabolites_included_in_reactions = set()
        for reaction in reactions_with_positive_flux:
            for met in self.stoichiometric_dict:
                if reaction in self.stoichiometric_dict[met]:
                    metabolites_included_in_reactions.add(met)
        for met in metabolites_included_in_reactions:
            if abs(
                sum(self.stoichiometric_dict[met][rxn] * self.route_model.constant_reaction_fluxes[rxn].value
                    for rxn in reactions_with_positive_flux if rxn in self.stoichiometric_dict[met]) +
                sum(self.temp_stoichiometric_dict[met][rxn] * self.route_model.temporary_exchange_fluxes[rxn].value
                    for rxn in reactions_with_positive_flux if met in self.temp_stoichiometric_dict
                    and rxn in self.temp_stoichiometric_dict[met]
                    and rxn in self.route_model.temporary_exchange_fluxes)
            ) > self.epsilon:
                warning_text = f'''
                Metabolite {self.cobra_model.metabolites[met].name} has a non-zero sum of fluxes, 
                larger than epsilon: {self.epsilon}, sum of fluxes: {sum(self.stoichiometric_dict[met][rxn] * self.route_model.constant_reaction_fluxes[rxn].value
                                    for rxn in reactions_with_positive_flux if rxn in self.stoichiometric_dict[met])}'''
                warnings.warn(warning_text)
                return False
        return True

    # noinspection PyTestUnpassedFixture
    def create_new_filtered_cobra_model(self) -> (Model, dict):
        filtered_model = Model()
        reaction_flux_dict = {reaction: self.route_model.constant_reaction_fluxes[rxn].value
                              for rxn, reaction in enumerate(self.cobra_model.reactions)
                              if not reaction.id.startswith("temporary_exchange") and abs(self.route_model.constant_reaction_fluxes[rxn].value > self.epsilon)}
        reaction_flux_dict.update({reaction: self.route_model.temporary_exchange_fluxes[rxn].value
                                   for rxn, reaction in enumerate(self.cobra_model.reactions)
                                   if reaction.id.startswith("temporary_exchange")})
        for reaction in self.cobra_model.reactions:
            if reaction in reaction_flux_dict:
                new_reaction = Reaction(reaction.id)
                new_reaction.name = reaction.name
                new_reaction.subsystem = reaction.subsystem
                new_reaction.lower_bound = reaction.lower_bound
                new_reaction.upper_bound = reaction.upper_bound
                new_reaction.add_metabolites(reaction.metabolites)
                filtered_model.add_reactions([new_reaction])
            elif reaction.id.startswith("temporary_exchange_"):
                new_reaction = Reaction(reaction.id)
                new_reaction.name = reaction.name
                new_reaction.subsystem = reaction.subsystem
                new_reaction.lower_bound = reaction.lower_bound
                new_reaction.upper_bound = reaction.upper_bound
                new_reaction.add_metabolites(reaction.metabolites)
                filtered_model.add_reactions([new_reaction])
        return filtered_model, reaction_flux_dict

    @staticmethod
    def is_newly_created_model_feasible(filtered_model: Model) -> bool:
        solution_ = filtered_model.optimize()
        if solution_.status != "optimal":
            return False
        return True

    @staticmethod
    def create_minimized_filtered_model(
            filtered_model: Model,
            reaction_flux_dict: dict
    ) -> (Model, dict):
        solution = filtered_model.optimize('minimize')
        fluxes = solution.fluxes
        for idx, reaction in enumerate(reaction_flux_dict):
            if idx < len(fluxes):
                reaction_flux_dict[reaction] = fluxes.iloc[idx]
        reactions_to_remove = [reaction.id for reaction in reaction_flux_dict if abs(reaction_flux_dict[reaction]) <= 0]
        reactions_to_remove = [filtered_model.reactions.get_by_id(reaction) for reaction in reactions_to_remove]
        filtered_model.remove_reactions(reactions_to_remove)
        reaction_flux_dict = {reaction: reaction_flux_dict[reaction] for reaction in reaction_flux_dict if reaction not in reactions_to_remove}
        return filtered_model, reaction_flux_dict

    @staticmethod
    def save_with_tries(
            file_name: str,
            object_to_save,
            funct,
            overwrite_run_files = False,
            requires_file_writing=False
    ):
        file_name_to_use = f"{file_name}"
        iterator = 1
        ending = file_name.split(".")[-1]
        if not overwrite_run_files:
            while True:
                if os.path.exists(file_name_to_use):
                    file_name_to_use = f"{file_name.split('.')[0]}_{iterator}.{ending}"
                    iterator += 1
                elif iterator > 1000:
                    warnings.warn(f"Could not save file {file_name_to_use}")
                    break
                else:
                    break
        if requires_file_writing:
            with open(file_name_to_use, "w") as file:
                funct(object_to_save, file)
        else:
            funct(object_to_save, file_name_to_use)

    def save_filtered_model(self, task_number: int, sample_number: int, filtered_model: Model):
        filtered_model_name = (f"{self.result_output_folder}\\{self.name_of_run}\\"
                               f"model_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.json")
        self.save_with_tries(filtered_model_name, filtered_model, save_json_model, self.overwrite_run_files)
        filtered_model_matlab_name = (f"{self.result_output_folder}\\{self.name_of_run}\\"
                                      f"model_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.mat")
        self.save_with_tries(filtered_model_matlab_name, filtered_model, save_matlab_model,self.overwrite_run_files)

    def convert_reaction_formula(self, reaction_formulas: list) -> list:
        # add stoichiometry and metabolite names
        converted_formulas = []
        for formula in reaction_formulas:
            metabolites_in_formula = formula.split(" ")
            form = []
            for idx, met in enumerate(metabolites_in_formula):
                if met.startswith("MAM"):
                    met_in_model = self.cobra_model.metabolites.get_by_id(met)
                    new_met_name = met_in_model.name + f"[{met_in_model.compartment}]"
                    form.append(new_met_name)
                else:
                    form.append(met)
            form = " ".join(form)
            converted_formulas.append(form)
        return converted_formulas

    def get_objective_penalties(self, reaction_ids: list) -> list:
        objective_penalties = [self.route_model.objective.args[0].args[self.cobra_model.reactions.index(self.cobra_model.reactions.get_by_id(reaction))] for
                               reaction in reaction_ids if not reaction.startswith("temporary_exchange_")]
        for idx, penalty in enumerate(objective_penalties):
            objective_penalties[idx] = float(str(penalty).split("*")[0])
        return objective_penalties

    def get_work_units(self) -> (float, float):
        log_data = self.persistent_gurobi._log
        # find the last occurrence of (:f work units) where :f is a float
        # grab this float from the str log_data
        pattern = r'\(\d+\.\d+ work units\)'
        matches = re.findall(pattern, log_data)
        if len(matches) > 1:
            root_work_units = float(matches[-2].split(" ")[0][1:])
            optimization_work_units = float(matches[-1].split(" ")[0][1:])
            return optimization_work_units, root_work_units
        elif len(matches) == 1:
            return float(matches[0].split(" ")[0][1:]), 0
        else:
            return 0, 0

    # noinspection PyTestUnpassedFixture
    def create_filtered_excel_csv_json_files(self, task_number: int, sample_number: int, reaction_flux_dict: dict):
        reaction_expression_dict = {reaction: self.route_model.reaction_expression[self.cobra_model.reactions.index(reaction)].value for rxn, reaction in enumerate(reaction_flux_dict) if
                                    not reaction.id.startswith("temporary_exchange_")}
        reaction_expression_dict.update({reaction: -1 for rxn, reaction in enumerate(reaction_flux_dict) if reaction.id.startswith("temporary_exchange_")})
        reaction_ids = [reaction.id for reaction in reaction_flux_dict]
        reaction_formulas = [reaction.reaction for reaction in reaction_flux_dict]
        reaction_names = [reaction.name for reaction in reaction_flux_dict]
        converted_formulas = self.convert_reaction_formula(reaction_formulas)
        reaction_objective_penalties = self.get_objective_penalties(reaction_ids)
        reaction_objective_penalties.extend([0 for _ in range(len(reaction_ids) - len(reaction_objective_penalties))])
        reaction_activity = [
            self.route_model.constant_reaction_is_active[self.cobra_model.reactions.index(reaction)].value
            for rxn, reaction in enumerate(reaction_flux_dict) if
               not reaction.id.startswith("temporary_exchange_") and
               reaction in self.cobra_model.reactions
        ]
        reaction_activity.extend([1 for _ in range(len(reaction_ids) - len(reaction_activity))])

        reaction_ids_only_forward = [reaction.id.split("_r")[0]+"_f" if reaction.id.endswith("_r") else reaction.id for reaction in reaction_flux_dict]
        forward_reaction_indices = [self.cobra_model.reactions.index(
            self.cobra_model.reactions.get_by_id(reaction_id)) for reaction_id in reaction_ids_only_forward]
        per_reaction_gene_contribution = [
            self.per_reaction_adjusted_gene_contribution[f"sample_{sample_number}"].values[rxn] for rxn in forward_reaction_indices
            if rxn < len(self.per_reaction_adjusted_gene_contribution[f"sample_{sample_number}"].values)
        ]
        per_reaction_gene_contribution.extend([{} for _ in range(len(reaction_ids) - len(per_reaction_gene_contribution))])

        reaction_df = pd.DataFrame()
        reaction_df["Flux"] = reaction_flux_dict.values()
        reaction_df["Expression"] = reaction_expression_dict.values()
        reaction_df["ID"] = reaction_ids
        reaction_df["Name"] = reaction_names
        reaction_df["Formula"] = reaction_formulas
        reaction_df["Converted_formula"] = converted_formulas
        reaction_df["Objective_penalty"] = reaction_objective_penalties

        options = ["Termination condition", "Root work units", "Optimization", "Time", "Objective"]
        optimization_work_units, root_work_units = self.get_work_units()
        option_values = [
            self.result.Solver.termination_condition,
            root_work_units,
            optimization_work_units,
            self.time_for_current_task,
            value(self.route_model.objective),
        ]
        options.extend(["" for _ in range(len(converted_formulas) - len(option_values))])
        option_values.extend(["" for _ in range(len(converted_formulas) - len(option_values))])
        amount_of_empty_rows_to_add = len(option_values) - len(converted_formulas) #TODO when reactions are fewer
        if amount_of_empty_rows_to_add > 0:
            reaction_df = reaction_df.reindex(reaction_df.index.tolist() + list(
                range(len(reaction_df.index),
                      len(reaction_df.index)+amount_of_empty_rows_to_add))
                                )
        reaction_df["Options"] = options
        reaction_df["Option_values"] = option_values

        filename = f"{self.result_output_folder}\\{self.name_of_run}\\filtered_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.xlsx"
        self.save_with_tries(filename, reaction_df, lambda x, y: x.to_excel(y, index=False),self.overwrite_run_files)

        reaction_df_reaction_flux = reaction_df[["ID", "Flux"]]
        filename = f"{self.result_output_folder}\\{self.name_of_run}\\esher_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.csv"
        self.save_with_tries(filename, reaction_df_reaction_flux, lambda x, y: x.to_csv(y, index=False), self.overwrite_run_files)

        # TODO add contributing gene

        results_dict = {
            "Task": task_number,
            "Sample": sample_number,
            "Model_type": self.model_type,
            "Name_of_run": self.name_of_run,
            "Combination_name": self.combinations_folder.split("\\")[-1].split("/")[-1],
            "Objective": value(self.route_model.objective),
            "Termination_condition": self.result.Solver.termination_condition,
            "Solution_is_feasible":  self.solution_is_feasible,
            "Solution_respects_SV_constraints": self.solution_respects_SV_constraints,
            "Newly_created_model_is_feasible": self.newly_created_model_is_feasible,
            "Root_work_units": root_work_units,
            "Optimization_work_units": optimization_work_units,
            "Time": self.time_for_current_task,
            "Expressionless_reaction_value": self.expressionless_reaction_value,
            "Epsilon": self.epsilon,
            "Unbounded_lower_bound_value": self.unbounded_lower_bound_value,
            "Unbounded_upper_bound_value": self.unbounded_upper_bound_value,
            "Solver_settings": self.solver_settings,
            "Transformation": self.transformation,
            "Thresholding": self.thresholding,
            "GPR_method": self.GPR_method,
            "Promiscuity_method": self.promiscuity_method,
            "Tasklist": self.tasklist_name,
            "Model_version": self.model_name,
            "Expression_dataset": self.expression_dataset_name,
            "Reactions": reaction_ids,
            "Reaction_fluxes": list(reaction_flux_dict.values()),
            "Reaction_expression": list(reaction_expression_dict.values()),
            "Reaction_objective_penalties": reaction_objective_penalties,
            "Reaction_names": reaction_names,
            "Reaction_formulas": reaction_formulas,
            "Converted_formulas": converted_formulas,
            "Version_of_code": version,
            "Optimization values": reaction_activity,
            "Contributing_gene_per_reaction":per_reaction_gene_contribution,
            "Date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        filename = f"{self.result_output_folder}\\{self.name_of_run}\\results_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.json"
        self.save_with_tries(filename, results_dict, dump, self.overwrite_run_files, requires_file_writing=True)

    def save_results(self, task_number: int, sample_number: int):
        (filtered_model, reaction_flux_dict) = self.create_new_filtered_cobra_model()
        self.newly_created_model_is_feasible = self.is_newly_created_model_feasible(filtered_model)
        if not self.newly_created_model_is_feasible:
            warnings.warn(f"Model for task {task_number} and sample {sample_number} is not feasible")
        if self.minimize_filtered_model:
            (filtered_model, reaction_flux_dict) = self.create_minimized_filtered_model(filtered_model, reaction_flux_dict)
        self.save_filtered_model(task_number, sample_number, filtered_model)
        self.create_filtered_excel_csv_json_files(task_number, sample_number, reaction_flux_dict)

    def create_unfeasible_results(self, task_number: int, sample_number: int):
        results_dict = {
            "Task": task_number,
            "Sample": sample_number,
            "Model_type": self.model_type,
            "Name_of_run": self.name_of_run,
            "Combination_name": self.combinations_folder.split("\\")[-1].split("/")[-1],
            "Objective": 0,
            "Termination_condition": self.result.Solver.termination_condition,
            "Solution_is_feasible": self.solution_is_feasible,
            "Solution_respects_SV_constraints": self.solution_respects_SV_constraints,
            "Newly_created_model_is_feasible": self.newly_created_model_is_feasible,
            "Root_work_units": 0,
            "Optimization_work_units": 0,
            "Time": self.time_for_current_task,
            "Expressionless_reaction_value": self.expressionless_reaction_value,
            "Epsilon": self.epsilon,
            "Unbounded_lower_bound_value": self.unbounded_lower_bound_value,
            "Unbounded_upper_bound_value": self.unbounded_upper_bound_value,
            "Solver_settings": self.solver_settings,
            "Transformation": self.transformation,
            "Thresholding": self.thresholding,
            "GPR_method": self.GPR_method,
            "Promiscuity_method": self.promiscuity_method,
            "Tasklist": self.tasklist_name,
            "Model_version": self.model_name,
            "Expression_dataset": self.expression_dataset_name,
            "Reactions": [],
            "Reaction_fluxes": [],
            "Reaction_expressions": [],
            "Reaction_objective_penalties": [],
            "Reaction_names": [],
            "Reaction_formulas": [],
            "Converted_formulas": [],
            "Version_of_code": version,
            "Optimization values": [],
            "Contributing_gene_per_reaction":[],
            "Date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        return results_dict

    def main_optimization_loop(self):
        print("Starting optimization loop")
        start_time = time.time()
        self.persistent_gurobi = SolverFactory("gurobi")
        for key, value_ in self.solver_settings.items():
            self.persistent_gurobi.options[key] = value_

        next_task = True
        sample_counter = 0
        task_counter = -1
        if not os.path.exists(self.result_output_folder):
            os.makedirs(self.result_output_folder)
        if not os.path.exists(f"{self.result_output_folder}\\{self.name_of_run}"):
            os.makedirs(f"{self.result_output_folder}\\{self.name_of_run}")

        while True:
            if sample_counter == len(self.samples_to_run) or next_task:
                task_counter += 1
                if task_counter == len(self.tasks_to_run):
                    break
                next_task = False
                sample_counter = 0
                task_number = self.tasks_to_run[task_counter]
                if not self.check_if_task_exists(task_number):
                    next_task = True
                    continue
                self.create_task_specific_models(task_number)
                if not self.check_if_task_is_viable():
                    warning_message = f"Task {task_number} is not viable (cannot perform cobra.optimize)."
                    warnings.warn(warning_message)
                    next_task = True
                    continue

            sample_number = self.samples_to_run[sample_counter]
            sample_counter += 1
            if not self.check_if_sample_exists(sample_number) and not self.ignore_sample_warnings:
                warning_message = f"Sample {sample_number} does not exist in the expression data."
                warnings.warn(warning_message)
                continue
            if self.skip_existing_files:
                if os.path.exists(f"{self.result_output_folder}\\{self.name_of_run}\\results_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.json"):
                    warning_message = (f"Sample {sample_number} and task {task_number} already exists in the result folder."
                                       f"Skip_existing_files is set to True, so this sample will be skipped.")
                    warnings.warn(warning_message)
                    continue
            print(f"Starting task {task_number} and sample {sample_number}")
            # noinspection PyTypeChecker, PyUnboundLocalVariable,PyTestUnpassedFixture
            self.set_log_file_name(task_number, sample_number)
            self.set_sample_specific_expression(sample_number)
            self.route_model.del_component(self.route_model.objective)
            # noinspection PyTypeChecker,  PyTestUnpassedFixture
            self.route_model.objective = Objective(expr=quicksum(
                self.route_model.constant_reaction_is_active[rxn] / self.route_model.reaction_expression[rxn].value for rxn in
                self.route_model.constant_reactions), sense=minimize)
            next_task = False
            self.result = self.persistent_gurobi.solve(self.route_model, tee=True)
            self.time_for_current_task = time.time() - start_time
            print(f"Task {task_number} and sample {sample_number} took {self.time_for_current_task} seconds")
            start_time = time.time()
            self.solution_is_feasible = self.check_model_viability()
            self.solution_respects_SV_constraints = self.check_if_optimization_model_respects_SV_constraints()
            if self.solution_is_feasible and self.solution_respects_SV_constraints:
                time_to_save = time.time() - start_time
                self.save_results(task_number, sample_number)
                print(f"Task {task_number} and sample {sample_number} has been saved in {time_to_save} seconds")
            else:
                results_dict = self.create_unfeasible_results(task_number, sample_number)
                filename = f"{self.result_output_folder}\\{self.name_of_run}\\results_Task_{task_number:03}_Sample{sample_number:03}_{self.model_type}.json"
                self.save_with_tries(filename, results_dict, dump, self.overwrite_run_files, requires_file_writing=True)

if __name__ == "__main__":
    start_task =2
    end_task = 2

    methods_dict_for_improved = {
        "__global_options": {
            "epsilon": 1e-4,
            "expressionless_reaction_value": 1.14,
            "create_log_file": True,
            "tasks": list(range(int(start_task), int(end_task + 1))),
            "samples": [1],  # list(range(1, 1010)),
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
            "tasklist": "IMP_tasklist",
            "model_version": "in_house",
            "expression_dataset": "DCM_test",
            "alternative_reaction_notation": False,
            "result_output_folder": r"C:\Users\inapa\PycharmProjects\Project_test_Ioana\Data\Results_test\\",
            "name_of_run": "test",
            "create_csv_files": True,
            "overwrite_run_files": True,
            # if True, will not try to create new file names but instead overwrite existing files
            "skip_existing_files": False,  # if True, will skip the task and sample if the files already exist
            "version_of_code": version,
            "ignore_sample_warnings": True,
            "minimize_filtered_model": True,
        },
    }
    main_data_folder = r"C:\Users\inapa\PycharmProjects\Project_test_Ioana\Data\Main_files"
    models_object = ModelCreator(methods_dict_for_improved, main_data_folder)
    models_object.main_optimization_loop()
