import Legacy
from Legacy import Model_formulations as modform_legacy, cModels as cM
from Legacy.sResults import find_or_create_user_based_directory


user = "Jelle"
main_data_folder, output_folder = find_or_create_user_based_directory(
    "test_runs",
    user,
    "INIT_runs_for_testing",
    # base_files= True,
)

Task = 1
name_of_run = f"INIT_compare_rev_vs_irrev{Task}"
# test_array = [10]
test_array = [Task]
sample_array = [1]
methods_dict = {
    "INIT_irrev": (
        "irreversible",
        modform_legacy.create_init_like_model_irrev_legacy,
        {
            "epsilon": 1e-3,
            "expressionless_reaction_value": 0.3,
            "Fix temporary exchange reactions to active": False,
        },
    ),
    # "INIT_rev_eps_5": ("reversible",
    #             modform.create_init_like_model_rev_legacy,
    #             {"epsilon": 1e-5,
    #              "expressionless_reaction_value": 0.3},
    #             ),
    #
    # "INIT_rev_eps_3": ("reversible",
    #              modform.create_init_like_model_rev,
    #              {"epsilon": 1e-5,
    #               "expressionless_reaction_value": 0.3},
    #              ),
}

cM.run_test_suite_different_model_formulations_legacy(
    name_of_run, test_array, sample_array, methods_dict, main_data_folder, output_folder
)

# 250 tasks and 332 samples
###
# how much does fastCC filtering affect performance (fastcc per task takes about 300-600s)
# how much does reversible vs irreversible affects (creating irreversible takes 40s)
# how much does epsilon e-4 e-5 e-6
# how much does model formulation matter:
#           for irreversible V, Y; V, Z, Y; V, S, Y
#           for reversible V, Z, Y; V, S, Z, Y
#  temporary exchange reactions (task specific input/output) if reaction.V.lb > 0, rxn_Y_active lb, ub = 1
#  how much difference in time for same task with different samples
#
#


# median value should be dependent on sample

# add in value for sample number and run into funciton and saving as well

# add passTask and check if folder exists, otherwise skip
