import pyomo.environ as pe
import numpy as np
from .cModels import *
import cobra

__all__ = ["create_init_like_model_irrev_legacy", "create_init_like_model_rev_legacy"]


def create_init_like_model_irrev_legacy(
    cobra_model: cobra.Model,
    expression_original: list,
    rev2irrev: list,
    reversible_expression: list,
    **kwargs,
) -> (pe.ConcreteModel, list):
    # Constant part
    stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(cobra_model)
    model, stoich_dict = create_basic_pyomo_model_SV_constraints(
        stoichiometric_matrix, cobra_model
    )
    # model = sap.add_loopless(cobra_model, model)
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
    model.Y_rxn_active = pe.Var(model.Rxns, domain=pe.Binary)
    if "Fix temporary exchange reactions to active" in kwargs:
        if kwargs["Fix temporary exchange reactions to active"] == True:
            for rxn in cobra_model.reactions:
                if rxn.lb > 0:
                    model.Y_rxn_active[rxn.id] = 1

    def reaction_constraint_rule(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0

    def reaction_constraint_rule_lower(model, rxn):
        # noinspection PyTypeChecker
        return (
            model.V_flux[rxn] - model.Y_rxn_active[rxn] * 1000 <= 0
        )  # should be UB of V

    def reaction_irreversible_reaction_pair_rule(model, rxn):
        # noinspection PyTypeChecker
        return (
            model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn - 1][1]] <= 1
            if rxn <= len(rev2irrev) and len(rev2irrev[rxn - 1]) > 1
            else pe.Constraint.Skip
        )

    def objective_rule(model):
        # noinspection PyTypeChecker
        return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint = pe.Constraint(model.Rxns, rule=reaction_constraint_rule)
    model.reaction_constraint_lower = pe.Constraint(
        model.Rxns, rule=reaction_constraint_rule_lower
    )
    model.reaction_irreversible_reaction_pair = pe.Constraint(
        model.Rxns, rule=reaction_irreversible_reaction_pair_rule
    )
    model.objective = pe.Objective(rule=objective_rule, sense=pe.minimize)

    return model, expression


def create_init_like_model_rev_legacy(
    cobra_model: cobra.Model, expression_original: list, **kwargs
) -> (pe.ConcreteModel, list):
    # Constant part
    stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(cobra_model)
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
    model.Y_rxn_active = pe.Var(model.Rxns, domain=pe.Binary)
    model.S_counter_flux = pe.Var(model.Rxns, domain=pe.Reals)

    def reaction_constraint_S_V_linkage(model, rxn):
        # noinspection PyTypeChecker
        return model.V_flux[rxn] + model.s_counter_flux[rxn] == 0

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
        # noinspection PyTypeChecker
        return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)

    print("Adding constraints")
    model.reaction_constraint_S_V_linkage = pe.Constraint(
        model.Rxns, rule=reaction_constraint_S_V_linkage
    )
    model.reaction_constraint_V_Y = pe.Constraint(
        model.Rxns, rule=reaction_constraint_V_Y
    )
    model.reaction_constraint_S_Y = pe.Constraint(
        model.Rxns, rule=reaction_constraint_S_Y
    )
    model.objective = pe.Objective(rule=objective_rule, sense=pe.minimize)

    return model, expression
