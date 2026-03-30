Metabolic tasks scripts:
The scripts up to the 2022_2 versions function as their equivalents from the CobraToolbox. 
The only added functionnalities are:
- fixing an important issue with getting essential tasks.
- the 'saveUsed' parameter, which allows to return all the reactions carying 
flux to perform a task during the essential reaction check. 
- the ability to use the EQU field from the RAVEN toolbox when designing tasks.
- the ability to allow any input/output by using "ALLMETS"

Starting from the 2023_3 version, the checkMetabolicTasksHumanGEM_2023_2 function now works by specifying every argument (aside from the model) as fields of a param structure. This should avoid accidents caused by the function's hight number of arguments. It also allows for the use of a relaxation parameter to define the strictness of the essentiality check.

Cellfie_2022:
Computes the Cellfie scores for metabolic tasks from gene expression data, based on the essential reactions needed for the task in a reference model.

cobra_to_metNames_tasks.R: R script which can be used to convert our COBRA tasklists to a tasklist using metNames instead of IDs, useful for visual inspection. Could be modified to allow conversion other COBRA tasklists