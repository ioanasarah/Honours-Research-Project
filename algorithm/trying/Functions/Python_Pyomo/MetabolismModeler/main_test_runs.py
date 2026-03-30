from Functions.Python_Pyomo.MetabolismModeler.sResults import (
    find_or_create_user_based_directory,
)
from Functions.Python_Pyomo.MetabolismModeler import Analysis as an
from Functions.Python_Pyomo.MetabolismModeler import Model_formulations as modform
import cProfile

user = "Jelle"
main_data_folder, output_folder = find_or_create_user_based_directory(
    "test_runs",
    user,
    "INIT_runs_for_testing",
    base_files=True,
)


Task = ""

name_of_run = f"testing_performance_code"
# test_array = [10]
# test_array = list(range(2,10)) + list(range(50,100)) + list(range(240, 260)) + list(range(320, 330))
# sample_array = list(range(1, 50))
task_array = [1010]
sample_array = [1, 99]# list(range(1,332))

methods_dict_for_improved = {
    "__global_options": {
        "transformation": "NT",
        "thresholding": "NT",
        "GPR_method": "minSum",
        "promiscuity_method": "Y",
        "model_type": "irr",
        "check_if_input_files_correct": True,
        "tasklist": "metabolic_tasks_2024_v1.01",
        "model_version": "consensus_2024",
        "expression_dataset": "DCM_magnet",
        "alternative_reaction_notation": False,
    },
}
'''
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
                - transformation:   NT (no transformation),
                                    L2T (log2 transformation),
                                    L10T (log10 transformation),
                                    LNT (ln transformation),
                                    ZTL (z-score transformation local) only using data in sample,
                                    ZTG (z-score transformation global) using all data in expression
'''

methods_dict = {
    "__global_options": {
        "transformation": "NT",
        "thresholding": "NT",
        "GPR_method": "minSum",
        "promiscuity_method": "Y",
        "check_if_input_files_correct": True,
        "tasklist": "metabolic_tasks_2024_v1.01",
        "model_version": "consensus_2024",
        "expression_dataset": "DCM_magnet",
        "expression_type":
            # "local_threshold",  # "global_percentile_90", "global_no_threshold", "local_thereshold"
            "global_no_threshold",  # "global_percentile_90", "global_no_threshold", "local_thereshold"
    },
    "INIT_irrev": (
        "irreversible",
        modform.create_init_like_model_irrev,
        {
            "epsilon": 1e-4,
            "expressionless_reaction_value": 0.3,
            "create_log_file": True,
            "Fix temporary exchange reactions to active": False,
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
        },
    ),
}

methods_dict_bottleneck = {
    "__global_options": {
        "expression_type":
        # "local_threshold",  # "global_percentile_90", "global_no_threshold", "local_thereshold"
            "global_no_threshold",  # "global_percentile_90", "global_no_threshold", "local_thereshold"
    },
    "bottleneck_1_INIT_irrev": (
        "irreversible",
        modform.bottleneck_maximization,
        {
            "epsilon": 1e-4,
            "expressionless_reaction_value": 0.3,
            "create_log_file": True,
            "Fix temporary exchange reactions to active": False,
            "adjust_min_expression_to": 1,

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
        },
    ),
}
# cProfile.run("an.run_test_suite_different_model_formulations(name_of_run, task_array, sample_array, methods_dict, main_data_folder, output_folder)")

an.run_test_suite_different_model_formulations(
    name_of_run, task_array, sample_array, methods_dict, main_data_folder, output_folder
)
# an.run_test_suite_different_model_formulations(
#     name_of_run, task_array, sample_array, methods_dict_bottleneck, main_data_folder, output_folder
# )

# test_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
# sample_array = [2, 10, 62, 104] #45 did not converge !!! TODO why
# sample_array = [1]
# # test_array = list(range(27,351))
# # test_array = [1,2,23,24,29,30,31,32,41, 104, 105, 106, 120, 251, 261, 269 , 299, 300, 301, 328]
# # test_array = [108, 160, 219, 232, 281, 282, 333, 334, 335]
# # test_array = [106]
# # test_array = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110]
# test_array = [ 106]
# epsilon_options = [1e-3, 1e-4, 1e-5, 1e-6]
# feasibility_tol_options = [1e-5, 1e-6,1e-7, 1e-8, 1e-9]
# intfeas_tol_options = [1e-7, 1e-8, 1e-9]
# seeds = [1, 2, 3, 4, 5]
# already_done = [
#                 [1e-6, 1e-5, 1e-5, 1], [1e-6, 1e-5, 1e-6, 1], [1e-6, 1e-5, 1e-7, 1],
#                 [1e-6, 1e-5, 1e-8, 1], [1e-6, 1e-5, 1e-9, 1],
#
# ]
# for seed in seeds:
#     for epsilon in epsilon_options:
#         for feasibility_tol in feasibility_tol_options:
#             for intfeas_tol in intfeas_tol_options:
#                 if [epsilon, intfeas_tol, feasibility_tol, seed] in already_done:
#                     continue
#                 methods_dict = {
#                     "__global_options": {
#                         "expression_type":
#                             "global_no_threshold",
#                             # "local_threshold",  # "global_percentile_90", "global_no_threshold", "local_threshold"
#                     },
#                     f"INIT_irrev{epsilon}{intfeas_tol}{feasibility_tol}": (
#                         "irreversible",
#                         modform.create_init_like_model_irrev,
#                         {
#                             "epsilon": epsilon,
#                             "expressionless_reaction_value": 0.3,
#                             "create_log_file": True,
#                             "Fix temporary exchange reactions to active": False,
#                             "solver_settings": {
#                                 "TimeLimit": 3600,
#                                 "IntFeasTol": intfeas_tol,
#                                 "FeasibilityTol": feasibility_tol,
#                                 "Seed": seed,
#                                 # "MIPFocus": 2,
#                                 # "Heuristics": 0.01,
#                                 # "Cuts": 3,
#                                 # "Presolve": 2,
#                                 # "ScaleFlag": 3,
#                             },
#                         },
#                     ),
#                 }
#                 an.run_test_suite_different_model_formulations(
#                     name_of_run, test_array, sample_array, methods_dict, main_data_folder, output_folder
#                 )
#
#
