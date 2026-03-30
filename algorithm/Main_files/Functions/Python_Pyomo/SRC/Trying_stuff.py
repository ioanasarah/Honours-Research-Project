import pyomo.environ as pe
import numpy as np
import pyomo.opt as po
from SRC import Helper_functions as hf
import cobra
import json

user = "Jelle"
main_data_folder = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\INIT_runs_for_testing"
main_data_folder, output_folder  = hf.find_or_create_user_based_directory("test_runs", user)
using_toy_model = False
epsilon = 1e-3

# [11, 6, 15, 24, 25, 88, 93, 85, 1, 96, 86, 2, 3, 4, 5, 84, 10] # test set in ascending order of solving time
# loop that puts in this as tasknumber
# create formulation for reversible with V S and Y variable
# compare solving time between reversible and irreversible formulations

version_name = "INIT_like"
task_number = 77
cobra_model = None
if not using_toy_model:
    # cobra_model, expression, rev2irrev = hf.load_irreversible_model_and_expression_from_mat_files(task_number, main_data_folder)
    main_folder = "E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\Base_files"
    expression, significance, task_structure = hf.load_json_data(main_folder)

    main_folder = f"E:\\Git\\Metabolic_Task_Score\\Data\\pyomo_python_files\\Base_files\\Task_{task_number:03}"
    expression = json.load(open(f"{main_folder}\\expression_irreversible_task{task_number:03}.json", "r"))
    # expression = json.load(open(f"{main_folder}\\task_specific_expression_task001.json", "r"))
    sample_number = 1
    expression =[value[sample_number-1] for value in expression]
    expression = np.array(expression).astype(float)


    # cobra_model = cobra.io.load_matlab_model(f"{main_folder}\\task_specific_model_task001.mat")
    # print(cobra_model.optimize())
    cobra_model = cobra.io.load_matlab_model(f"{main_folder}\\irrev_model_task{task_number:03}.mat")
    # print(cobra_model.optimize())

    rev2irrev = json.load(open(f"{main_folder}\\rev2irrev_task{task_number:03}.json", "r"))

    median_expression = np.median(expression)
    median_expression = 0.3

    expression_old = expression.copy()
    print("Median expression: ", median_expression)
    expression = [value if value != -1 else median_expression for value in expression] # Expressions include -1 values which are not valid in the optimization
                                                                                    # they indicate missing knowledge about reaction/gene linkage
    stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(cobra_model)
    print(len(stoichiometric_matrix))
    print(len(stoichiometric_matrix[0]))
else:
    # Define the stoichiometric matrix (example)
    stoichiometric_matrix = np.array([
        #22  23  23  4   1   1   2   2

        #10  10  10  10  0   0   0   0
        [1, -1,  0,  0, -1,  0,  0,  0], #= 0 # input
        [0,  1, -1,  0,  0,  0,  0,  0],
        [0,  0,  1, -1,  0,  0,  0,  1],  # output
        [0,  0,  0,  0,  1, -1,  0,  0],
        [0,  0,  0,  0,  0,  1, -1,  0],
        [0,  0,  0,  0,  0,  0,  1, -1]
    #min 1/22 1/23 1/23 1/4 1/1 1/1 1/2 1/2
        # max v1 = 10
               # minimize sum of Yi
        # Yi = binary
        # Vi - Yi * 0.001   >= 0
        # Vi - Yi * 1000    <= 0 if vi >0 then Yi = 1
        # M = 10000
        # epsilon = 0.0001

        # if v can be -1000, 1000

        # Si + Vi = 0 -10  10
        # Vi - Yi * 1000 <= 0
        # Si - Yi * 1000 <= 0

        # Zi >=  Vi
        # Zi >= -Vi
        # Zi -   Yi * 0.001   >= 0 if vi = 0 then zi = 0 then yi = 0
        # Zi -   Yi * 1000    <= 0 if vi > 0 then Yi


    ])
    # A goes in, C goes out
    # flux rxn1 = 5
    # flux rxn2 = 3
    # flux rxn3 = 3
    # flux rxn4 = 5
    # flux rxn5 = 2
    expression = [1] * len(stoichiometric_matrix[0])

# Define the model
model = pe.ConcreteModel()

num_metabolites = len(stoichiometric_matrix)  # Example number of metabolites
num_reactions = len(stoichiometric_matrix[0])  # Example number of reactions
print(num_metabolites, num_reactions, len(expression))
metabolites = list(range(1, num_metabolites + 1))
reactions = list(range(1, num_reactions + 1))

print("Adding decision variables and S matrix")
model.Metabolites = pe.RangeSet(1, num_metabolites)
model.Rxns = pe.RangeSet(1, num_reactions)
stoich_dict = {met: {rxn: stoichiometric_matrix[met, rxn] for rxn in range(num_reactions) if stoichiometric_matrix[met, rxn] != 0} for met in range(num_metabolites)}
#
# # Define the decision variabes
model.V_flux = pe.Var(model.Rxns, domain=pe.Reals)
model.Y_rxn_active = pe.Var(model.Rxns, domain=pe.Binary )
#
if using_toy_model:
    # set model.V_flux bounds to (0, 1000)
    # noinspection PyTypeChecker,PyUnresolvedReferences
    def function_for_bounds():
        for rxn in model.Rxns:
            model.V_flux[rxn].setlb(0)
            model.V_flux[rxn].setub(1000)
        model.V_flux[1].setlb(1)
        model.V_flux[1].setub(1)

    function_for_bounds()

else:
    # noinspection PyTypeChecker,PyUnresolvedReferences
    def function_for_bounds():
        for rxn in model.Rxns:
            model.V_flux[rxn].setlb(cobra_model.reactions[rxn - 1].lower_bound)
            model.V_flux[rxn].setub(cobra_model.reactions[rxn - 1].upper_bound)

    function_for_bounds()


def sv_zero_rule(model, met):
    # noinspection PyTypeChecker
    return sum(stoich_dict[met - 1][rxn] * model.V_flux[rxn + 1] for rxn in stoich_dict[met - 1]) == 0


def reaction_constraint_rule(model, rxn):
    # noinspection PyTypeChecker
    return model.V_flux[rxn] - (model.Y_rxn_active[rxn] * epsilon) >= 0

def reaction_constraint_rule_lower(model, rxn):
    # noinspection PyTypeChecker
    return model.V_flux[rxn] - model.Y_rxn_active[rxn] * 1000 <= 0 # should be UB of V

def reaction_irreversible_reaction_pair_rule(model, rxn):
    # noinspection PyTypeChecker
    return model.Y_rxn_active[rxn] + model.Y_rxn_active[rev2irrev[rxn-1][1]] <= 1 if rxn <= len(rev2irrev) and len(rev2irrev[rxn-1]) > 1 else pe.Constraint.Skip

def objective_rule(model):
    # noinspection PyTypeChecker
    return sum(model.Y_rxn_active[rxn]/expression[rxn - 1] for rxn in model.Rxns)


print("Adding constraints")
model.sv_zero_constraint = pe.Constraint(model.Metabolites, rule=sv_zero_rule)
print("added SV = 0 constraint")
model.reaction_constraint = pe.Constraint(model.Rxns, rule=reaction_constraint_rule)
model.reaction_constraint_lower = pe.Constraint(model.Rxns, rule=reaction_constraint_rule_lower)
model.reaction_irreversible_reaction_pair = pe.Constraint(model.Rxns, rule=reaction_irreversible_reaction_pair_rule)
model.objective = pe.Objective(rule=objective_rule, sense=pe.minimize)

# model formulation
# model.Y_rxn_active = pe.Var(model.Rxns, domain=pe.Binary)
# model.S_counter_flux = pe.Var(model.Rxns, domain=pe.Reals)

#
# def reaction_constraint_S_V_linkage(model, rxn):
#     # noinspection PyTypeChecker
#     return model.V_flux[rxn] + model.S_counter_flux[rxn] == 0
#
#
# def reaction_constraint_V_Y(model, rxn):
#     # noinspection PyTypeChecker
#     return (model.Y_rxn_active[rxn] * 1000) - model.V_flux[rxn] >= 0  # should be UB of V
#
#
# def reaction_constraint_S_Y(model, rxn):
#     # noinspection PyTypeChecker
#     return (model.Y_rxn_active[rxn] * 1000) - model.S_counter_flux[rxn] >= 0  # should be UB of V
#
#
# def objective_rule(model):
#     # noinspection PyTypeChecker
#     return sum(model.Y_rxn_active[rxn] / expression[rxn - 1] for rxn in model.Rxns)
#

# print("Adding constraints")
# model.reaction_constraint_sv_zero = pe.Constraint(model.Metabolites, rule=sv_zero_rule)
# model.reaction_constraint_S_V_linkage = pe.Constraint(model.Rxns, rule=reaction_constraint_S_V_linkage)
# model.reaction_constraint_V_Y = pe.Constraint(model.Rxns, rule=reaction_constraint_V_Y)
# model.reaction_constraint_S_Y = pe.Constraint(model.Rxns, rule=reaction_constraint_S_Y)
# model.objective = pe.Objective(rule=objective_rule, sense=pe.minimize)
# model.pprint()

solver = po.SolverFactory('gurobi')
instance = model.create_instance()
result = solver.solve(instance, tee=True)
# hf.save_solver_results(instance, result, cobra_model, f"{version_name}_solution_Task_{task_number}", output_folder, expression)
hf.save_solver_results(instance, result, cobra_model, f"{version_name}_solution_Task_{task_number:03}",
                    output_folder, expression, expression_old)
