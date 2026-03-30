import numpy as np
from Functions.Python_Pyomo.MetabolismModeler.cModels import (
    create_basic_pyomo_model_SV_constraints,
)
from cobra.util import create_stoichiometric_matrix
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
    value,
)
from pyomo.opt import SolverFactory
from cobra.core import Model
from collections.abc import Iterable
from math import ceil
import logging
__all__ = ["create_init_like_model_irrev", "create_init_like_model_rev"]


def create_init_like_model_irrev(
    cobra_model: Model,
    expression_original: list,
    rev2irrev: list,
    reversible_expression: list,
    **kwargs,
) -> (ConcreteModel, list):
    """
    Create a Pyomo model from a COBRA model with irreversible reaction constraints.

    Args:
        cobra_model (Model): The COBRA model.
        expression_original (list): List of expression values.
        rev2irrev (list): Mapping of reversible to irreversible reactions.
        reversible_expression (list): List of reversible expression values.
        **kwargs: Additional keyword arguments for customization.

    Returns:
        tuple: Pyomo model and processed expression list.
    """
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    expression = expression_original.copy()
    epsilon = kwargs.get("epsilon", 1e-3)
    expression = expression_original.copy()
    median_expression = np.median([value for value in expression if value != -1])
    median_expression = kwargs.get("expressionless_reaction_value", median_expression)

    # Replace -1 values with median expression for valid optimization
    expression = [value if value != -1 else median_expression for value in expression]

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    if kwargs.get("Fix temporary exchange reactions to active", False):
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].lower_bound > 0:
                model.Y_rxn_active[rxn].fix(1)

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0
            # -1

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
            model.V_flux[rxn] - model.Y_rxn_active[rxn] * 1000 <= 0
        )  # should be UB of V

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1
            if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else Constraint.Skip
        )

    def objective_rule(model):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )
    model.objective = Objective(rule=objective_rule, sense=minimize)

    return model, expression



def bottleneck_maximization(
        cobra_model: Model,
        expression_original: list,
        rev2irrev: list,
        reversible_expression: list,
        **kwargs,
) -> (ConcreteModel, list):
    """
    """
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    epsilon = kwargs.get("epsilon", 1e-3)
    expression = expression_original.copy()

    # if expression is higher than 200, it would imply that the expression is not log transformed
    # thus the value chosen should be higher
    if np.max(expression) >= 200:
        median_expression = kwargs.get("adjusted_expressionless_for_bottleneck", 1000)
    else:
        median_expression = kwargs.get("adjusted_expressionless_for_bottleneck", 10)

    # Replace -1 values with median expression for valid optimization
    expression = [value if value != -1 else median_expression for value in expression]

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    model.minimum_expression = Var(domain=NonNegativeReals)
    reaction_upper_bounds = {}
    for rxn in model.Rxns:
        reaction_upper_bounds[rxn] = cobra_model.reactions[rxn - 1].upper_bound

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0
        # -1

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
                model.V_flux[rxn] - model.Y_rxn_active[rxn] * reaction_upper_bounds[rxn] <= 0
        )  # should be UB of V

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1
            if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else Constraint.Skip
        )

    def minimum_expression_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.minimum_expression <= expression[rxn - 1] * model.Y_rxn_active[rxn] + (1 - model.Y_rxn_active[rxn]) * 1000
    def maximize_minimum_expression(model):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return model.minimum_expression

    print("Adding constraints")
    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )
    model.minimum_expression_constraint = Constraint(model.Rxns, rule=minimum_expression_constraint_rule)
    model.objective = Objective(rule=maximize_minimum_expression, sense=maximize)

    if kwargs.get("solver", False):
        solver = kwargs.get("solver", SolverFactory('gurobi'))
        if "solver_settings" in kwargs:
            for key, value in kwargs["solver_settings"].items():
                solver.options[key] = value
    else:
        solver = kwargs.get("solver", SolverFactory('gurobi'))

    results = solver.solve(model, tee=True)

    # if not optimimal, return the model and the expression
    if results.solver.termination_condition != "optimal":
        return model, expression

    ######## step 2

    print("Minimum expression: ", model.minimum_expression.value)
    minimized_expression = model.minimum_expression.value * kwargs.get("adjust_min_expression_to", 1)

    #### new model
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    epsilon = kwargs.get("epsilon", 1e-3)

    expression = expression_original.copy()
    expressionless_reaction_value = minimized_expression #TODO have to fix this? doesnt seem to work
    expression = [value if value != -1 else expressionless_reaction_value for value in expression]

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    model.minimum_expression = Var(domain=NonNegativeReals)
    reaction_upper_bounds = {}
    for rxn in model.Rxns:
        reaction_upper_bounds[rxn] = cobra_model.reactions[rxn - 1].upper_bound

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0
        # -1

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
                model.V_flux[rxn] - model.Y_rxn_active[rxn] * reaction_upper_bounds[rxn] <= 0
        )  # should be UB of V

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1
            if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else Constraint.Skip
        )

    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )

    def expression_must_be_higher_than_minimum_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.Y_rxn_active[rxn] * expression[rxn - 1] >= minimized_expression - (1 - model.Y_rxn_active[rxn]) * 1000

    model.expression_higher_than_minimum_constraint = Constraint(model.Rxns, rule=expression_must_be_higher_than_minimum_rule)

    def new_objective(model):
        return sum(model.Y_rxn_active[rxn]/ expression[rxn - 1] for rxn in model.Rxns)

    model.objective = Objective(rule=new_objective, sense=minimize)
    print("Added constraints to second step of bottleneck maximization")
    return model, expression

def bottleneck_maximization_old(cobra_model: Model, expression_original: list, rev2irrev: list, reversible_expression: list, **kwargs) -> (ConcreteModel, list):
    """
    Perform bottleneck maximization using a two-step approach:
    1. Maximize the minimum expression level across all active reactions.
    2. Minimize the total number of active reactions while ensuring that the minimum expression level is met.
    """

    # logging.basicConfig(level=logging.INFO)
    # logger = logging.getLogger(__name__)
    # Step 1: Creating the stoichiometric matrix and initializing the Pyomo model.
    # logger.info("Creating stoichiometric matrix and Pyomo model.")
    print("Creating stoichiometric matrix and Pyomo model.")
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(stoichiometric_matrix, cobra_model)

    # Step 2: Process expression levels
    expression = expression_original.copy()
    epsilon = kwargs.get("epsilon", 1e-3)
    min_expr = kwargs.get("min_expr", 10)  # Target minimum expression level
    expressionless_value = kwargs.get("expressionless_reaction_value", min_expr)
    expression = [value if value != -1 else expressionless_value for value in expression]  # Replace missing values
    max_expression = max(expression)
    # logger.info("Processed expression list: %s", expression)
    # print("Processed expression list: ", expression)

    # Ensure binary variable for reaction activity exists
    if not hasattr(model, 'Y_rxn_active'):
        model.Y_rxn_active = Var(model.Rxns, domain=Binary)

    # Step 3: Maximize minimum expression (first optimization step)
    # logger.info("Starting Step 1: Maximize min_expr.")
    print("Starting Step 1: Maximize min_expr.")
    model.min_expr = Var(domain=Reals, bounds=(min_expr, None))  # Variable for minimum expression

    def min_expression_constraint_rule(model, rxn):
        big_m = max_expression + 1  # Large value for handling inactive reactions
        return model.min_expr <= expression[rxn - 1] * model.Y_rxn_active[rxn] + (1 - model.Y_rxn_active[rxn]) * big_m
        # Ensures that the minimum expression meets the threshold for active reactions.

    model.min_expression_constraint = Constraint(model.Rxns, rule=min_expression_constraint_rule)

    def objective_step1(model):
        return model.min_expr  # Objective to maximize the minimum expression level across all active reactions

    model.objective_step1 = Objective(rule=objective_step1, sense=maximize)

    # Solve the model for the first objective
    solver = SolverFactory('gurobi')
    results = solver.solve(model, tee=True)
    # logger.info("Step 1 solver status: %s", results.solver.status)
    # logger.info("Step 1 solver termination condition: %s", results.solver.termination_condition)
    print("Step 1 solver status: ", results.solver.status)
    print("Step 1 solver termination condition: ", results.solver.termination_condition)

    # Check if the problem was infeasible or unbounded
    if results.solver.termination_condition == 'infeasibleOrUnbounded':
        print("Step 1: The model is infeasible or unbounded.")
        # logger.error("Step 1: The model is infeasible or unbounded.")
        raise ValueError("Step 1 failed: Model infeasible or unbounded.")

    optimized_min_expr = value(model.min_expr)  # Get the optimized minimum expression value
    # logger.info("Optimized min_expr: %f", optimized_min_expr)
    print("Optimized min_expr: ", optimized_min_expr)

    # Deactivate the first objective before moving to Step 2
    model.objective_step1.deactivate()

    # Step 4: Apply fractional INIT approach (second optimization step)
    # logger.info("Starting Step 2: Apply fractional INIT approach.")
    print("Starting Step 2: Apply fractional INIT approach.")
    flex_min_expr = kwargs.get("flex_min_expr", optimized_min_expr)

    def expression_must_be_higher_than_minimum_rule(model, rxn):
        big_m = max_expression + 1  # Large value to handle inactive reactions
        return model.Y_rxn_active[rxn] * expression[rxn - 1] >= flex_min_expr - (1 - model.Y_rxn_active[rxn]) * big_m
        # Ensures that active reactions meet the flex_min_expr threshold, and inactive reactions bypass this constraint.

    # Add constraints for the second step
    model.expression_update_constraint = Constraint(model.Rxns, rule=expression_must_be_higher_than_minimum_rule)

    def objective_step2(model):
        return sum(model.Y_rxn_active[rxn] for rxn in model.Rxns)  # Objective to minimize the total number of active reactions

    model.objective_step2 = Objective(rule=objective_step2, sense=minimize)

    # Solve the model for the second objective
    results = solver.solve(model, tee=True)
    # logger.info("Step 2 solver status: %s", results.solver.status)
    # logger.info("Step 2 solver termination condition: %s", results.solver.termination_condition)
    print("Step 2 solver status: ", results.solver.status)
    print("Step 2 solver termination condition: ", results.solver.termination_condition)

    # Check if the problem was infeasible or unbounded
    if results.solver.termination_condition == 'infeasibleOrUnbounded':
        # logger.error("Step 2: The model is infeasible or unbounded.")
        print("Step 2: The model is infeasible or unbounded.")
        raise ValueError("Step 2 failed: Model infeasible or unbounded.")

    return model, expression  # Return the optimized model and processed expression list
def create_init_like_model_rev(
    cobra_model: Model, expression_original: list, **kwargs
) -> (ConcreteModel, list):
    # Constant part
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )

    epsilon = 1e-3
    expression = expression_original.copy()
    median_expression = np.median([value for value in expression if value != -1])
    if "epsilon" in kwargs:
        epsilon = kwargs["epsilon"]
    if "expressionless_reaction_value" in kwargs:
        median_expression = kwargs["expressionless_reaction_value"]
    # Variable definition
    print("Median expression: ", median_expression)
    print(f"Amount of -1 values in expression: {expression.count(-1)}")
    expression = [
        value if value != -1 else median_expression for value in expression
    ]  # Expressions include -1 values which are not valid in the optimization

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    model.S_counter_flux = Var(model.Rxns, domain=Reals)

    def reaction_constraint_S_V_linkage(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] + model.S_counter_flux[rxn] == 0

    def reaction_constraint_V_Y(model, rxn):
        # noinspection PyTypeChecker
        return (model.Y_rxn_active[rxn] * 1000) - model.V_flux[
            rxn
        ] >= 0  # should be UB of V

    def reaction_constraint_S_Y(model, rxn):
        # noinspection PyTypeChecker
        return (model.Y_rxn_active[rxn] * 1000) - model.S_counter_flux[
            rxn
        ] >= 0  # should be UB of V

    def objective_rule(model):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint_S_V_linkage = Constraint(
        model.Rxns, rule=reaction_constraint_S_V_linkage
    )
    model.reaction_constraint_V_Y = Constraint(model.Rxns, rule=reaction_constraint_V_Y)
    model.reaction_constraint_S_Y = Constraint(model.Rxns, rule=reaction_constraint_S_Y)
    model.objective = Objective(rule=objective_rule, sense=minimize)

    return model, expression

def create_irrev_init_like_linear_penalty_model(
        cobra_model: Model,
        expression_original: list,
        rev2irrev: list,
        reversible_expression: list,
        **kwargs,
) -> (ConcreteModel, list):
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )

    epsilon = kwargs.get("epsilon", 1e-3)
    expression = expression_original.copy()
    median_expression = np.median([value for value in expression if value != -1])
    median_expression = kwargs.get("expressionless_reaction_value", median_expression)

    # Replace -1 values with median expression for valid optimization
    expression = [value if value != -1 else median_expression for value in expression]

    # Additional parameters
    weight_factor = kwargs.get("weight_factor", 1.0)
    reaction_activation_penalty = kwargs.get("reaction_activation_penalty", 0.0)
    custom_flux_bounds = kwargs.get("custom_flux_bounds", {})
    objective_scaling = kwargs.get("objective_scaling", 1.0)
    constraint_relaxation = kwargs.get("constraint_relaxation", 0.0)

    # Model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    if kwargs.get("Fix temporary exchange reactions to active", False):
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].lower_bound > 0:
                model.Y_rxn_active[rxn].fix(1)

    def reaction_constraint_rule(model, rxn):
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= -constraint_relaxation

    def reaction_constraint_rule_lower(model, rxn):
        ub = custom_flux_bounds.get(rxn, 1000)
        return model.V_flux[rxn] - model.Y_rxn_active[rxn] * ub <= constraint_relaxation

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1:
            return model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1
        return Constraint.Skip

    def objective_rule(model):
        expression_max = max(expression) + 1
        return objective_scaling * sum(
            model.Y_rxn_active[rxn] * (expression_max - expression[rxn - 1]) * weight_factor
            - reaction_activation_penalty * model.Y_rxn_active[rxn]
            for rxn in model.Rxns
        )

    # Add constraints and objective
    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(model.Rxns, rule=reaction_constraint_rule_lower)
    model.reaction_irreversible_reaction_pair = Constraint(model.Rxns, rule=reaction_irreversible_reaction_pair_rule)
    model.objective = Objective(rule=objective_rule, sense=minimize)

    return model, expression

def create_not_fancy_linear_init_penalty_no_fancy_irrev(
        cobra_model: Model,
        expression_original: list,
        rev2irrev: list,
        reversible_expression: list,
        **kwargs,
) -> (ConcreteModel, list):
    """
    Create a Pyomo model from a COBRA model with irreversible reaction constraints.

    Args:
        cobra_model (Model): The COBRA model.
        expression_original (list): List of expression values.
        rev2irrev (list): Mapping of reversible to irreversible reactions.
        reversible_expression (list): List of reversible expression values.
        **kwargs: Additional keyword arguments for customization.

    Returns:
        tuple: Pyomo model and processed expression list.
    """
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    epsilon = kwargs.get("epsilon", 1e-3)
    expression = expression_original.copy()
    median_expression = np.median([value for value in expression if value != -1])
    median_expression = kwargs.get("expressionless_reaction_value", median_expression)

    # Replace -1 values with median expression for valid optimization
    expression = [value if value != -1 else median_expression for value in expression]

    # Calculate linear expression penalty
    max_expression = np.max(expression) + 1
    scaling_factor = kwargs.get("scaling_factor", 1.0)
    expression_penalty = [((max_expression - value) * scaling_factor) for value in expression]
    if kwargs.get("set_max_expression_to_integer", False):
        rounded_up = ceil(max_expression)
        minimum_multiplier = rounded_up / (max_expression)
        expression_penalty = [value * minimum_multiplier for value in expression_penalty]
    print(f"epsilon: {epsilon}")
    if kwargs.get("round_expression_factor", False):
        expression_penalty = [round(value, kwargs["round_expression_factor"]) for value in expression_penalty]

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    if kwargs.get("Fix temporary exchange reactions to active", False):
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].lower_bound > 0:
                model.Y_rxn_active[rxn].fix(1)

    adjust_max_bounds = kwargs.get("adjust_max_bounds_factor", False)
    max_bound_big_M = 1000
    if adjust_max_bounds:
        print(f"Adjusting max bounds by {adjust_max_bounds}")
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].upper_bound >= 1000:
                model.V_flux[rxn].setub = 1000 * adjust_max_bounds
        max_bound_big_M = 1000 * adjust_max_bounds

    for keys, values in kwargs.items():
        if not isinstance(values, Iterable):
            print(keys, values)

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
                model.V_flux[rxn] - model.Y_rxn_active[rxn] * max_bound_big_M <= 0
        )

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1 if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else Constraint.Skip
        )

    def objective_rule(model):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return sum(model.Y_rxn_active[rxn] * expression_penalty [rxn - 1] for rxn in model.Rxns )
        # return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )
    model.objective = Objective(rule=objective_rule, sense=minimize)

    return model, expression

def create_irrev_minimze_reactions_without_expression(
        cobra_model: Model,
        expression_original: list,
        rev2irrev: list,
        reversible_expression: list,
        **kwargs,
) -> (ConcreteModel, list):
    """
    Create a Pyomo model from a COBRA model with irreversible reaction constraints.

    Args:
        cobra_model (Model): The COBRA model.
        expression_original (list): List of expression values.
        rev2irrev (list): Mapping of reversible to irreversible reactions.
        reversible_expression (list): List of reversible expression values.
        **kwargs: Additional keyword arguments for customization.

    Returns:
        tuple: Pyomo model and processed expression list.
    """
    stoichiometric_matrix = create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    epsilon = kwargs.get("epsilon", 1e-3)
    expression = expression_original.copy()
    median_expression = np.median([value for value in expression if value != -1])
    median_expression = kwargs.get("expressionless_reaction_value", median_expression)

    # Replace -1 values with median expression for valid optimization
    expression = [value if value != -1 else median_expression for value in expression]

    # Calculate linear expression penalty
    max_expression = np.max(expression) + 1
    scaling_factor = kwargs.get("scaling_factor", 1.0)
    expression_penalty = [((max_expression - value) * scaling_factor) for value in expression]
    if kwargs.get("set_max_expression_to_integer", False):
        rounded_up = ceil(max_expression)
        minimum_multiplier = rounded_up / (max_expression)
        expression_penalty = [value * minimum_multiplier for value in expression_penalty]
    print(f"epsilon: {epsilon}")
    if kwargs.get("round_expression_factor", False):
        expression_penalty = [round(value, kwargs["round_expression_factor"]) for value in expression_penalty]

    # model formulation
    model.Y_rxn_active = Var(model.Rxns, domain=Binary)
    if kwargs.get("Fix temporary exchange reactions to active", False):
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].lower_bound > 0:
                model.Y_rxn_active[rxn].fix(1)

    adjust_max_bounds = kwargs.get("adjust_max_bounds_factor", False)
    max_bound_big_M = 1000
    if adjust_max_bounds:
        print(f"Adjusting max bounds by {adjust_max_bounds}")
        for rxn in model.Rxns:
            if cobra_model.reactions[rxn - 1].upper_bound >= 1000:
                model.V_flux[rxn].setub = 1000 * adjust_max_bounds
        max_bound_big_M = 1000 * adjust_max_bounds

    for keys, values in kwargs.items():
        if not isinstance(values, Iterable):
            print(keys, values)

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0
        # flux  0.002 - 1 * 0.001 >= 0

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
                model.V_flux[rxn] - model.Y_rxn_active[rxn] * max_bound_big_M <= 0
        )

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1 if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else Constraint.Skip
        )

    def objective_rule(model):
        # noinspection PyTypeChecker, PyUnresolvedReferences
        return sum(model.Y_rxn_active[rxn] for rxn in model.Rxns )
        # return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint = Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )
    model.objective = Objective(rule=objective_rule, sense=minimize)

    return model, expression
