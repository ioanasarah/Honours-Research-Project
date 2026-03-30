#Convert to MAT

import cobra 

model = cobra.io.load_json_model(r"E:\Downloads\model_json_INIT_irrev_Task_010_Sample001_global_no_threshold.json")
cobra.io.save_matlab_model(model, r"E:\Downloads\model_json_INIT_irrev_Task_010_Sample001_global_no_threshold.mat")