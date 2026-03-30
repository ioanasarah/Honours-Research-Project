from cobra.io import (save_json_model,
                      save_matlab_model,
                      write_sbml_model,
                      load_json_model,
                      load_matlab_model,
                      read_sbml_model,
                      )



newModel_P3_3_3 = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3_3.3.mat")

newModel_P3 = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3.mat")

newModel_P3_3_2= load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3_3.2.mat") # error buffer

model_P3_IMP = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/model_P3_IMP.mat") # error buffer

save_json_model(newModel_P3_3_3,"C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3_3_3.json" )

newModel_P3_3_IMP = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3_3_IMP.mat")

save_json_model(newModel_P3_3_IMP,"C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/newModel_P3_3_IMP.json" )

model_w_alpha = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/model_with_alpha.mat")

model_w_alpha_ok = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/MSP project period jan 2025/prostaglandins/model_with_alpha_ok.mat")

concensus = load_matlab_model("C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/SysBio_COBRA_v1.17_consensus.mat")

save_json_model(concensus,"C:/Users/inapa/Documents/Repertoire etudiant_IMP/Projets/Projet 16 Metabolic task list/SysBio_COBRA_v1.17_consensus.json")

model_dummy=load_matlab_model("C:/Users/inapa/PycharmProjects/Metabolic_Task_Score/Data/Main_files/For_running/models/model_with_FA_dummy_reactions/model_1_17_with_dummy_FA_reactions.mat")
