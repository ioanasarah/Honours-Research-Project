# Honours Research Project: Modelling paracetamol toxicity
The basis of my Honour’s Research Project lies in constraint-based metabolic modelling, in which an algorithm maps gene expression data onto the Human1 Genome Scale Model to identify the most “optimal” pathway via a penalty minimisation calculation. This in silico modelling approach allows for the effects of medications to be modelled, providing a deeper understanding of their metabolic  and toxicity pathways without causing harm to humans or animals. This project has focused on mapping the pathways of paracetamol degradation and metabolic reactions affected by paracetamol toxicity.

## Explanation of folder contents

- [`escher_maps/`](./escher_maps)  
  Relevant Escher maps for visualizations of paracetamol-specific tasks.

- [`paracetamol_specific_expression_dataset/`](./paracetamol_specific_expression_dataset)  
  Gene expression datasets containing 7 samples, one healthy and 6 others at different levels of overdose. T0 is the healthy sample, T30min is 30 minutes of continuous paracetamol exposure, T3h is 3 hours of continuous paracetamol exposure, T6h is 6 hours of continuous paracetamol exposure,	T12h is 12 hours of continuous paracetamol exposure, T24h is 24 hours of continuous paracetamol exposure. T4d is 4 days of continuous paracetamol exposure, only the survivor cells from this cohort were analysed, post 21 days of outgrowth. (data taken from https://doi.org/10.1038/s41598-018-37940-6)
  

- [`paracetamol_specific_model/`](./paracetamol_specific_model)  
  Human1 GEM modified with more metabolites, reactions, and new gene rules to reflect a more accurate representation of the metabolism of paracetamol.

- [`tasklist_paracetamol_specific/`](./tasklist_paracetamol_specific)  
  Task definitions for paracetamol degradation and toxicity pathways, but also metabolic reactions which are known to be affected at toxic doses of the drug.

## Explanation of code
- [`adding_cobra_reactions.py`](./code/adding_cobra_reactions.py)  
Extends the Human1 GEM by modifying existing reactions and adding new ones to ensure the biologically accurate mapping of acetaminophen toxicity and glutathione metabolism. The code loads an existing .mat model using COBRApy, modifies it, and saves it.   

- [`comparing_expression_values.py`](./code/comparing_expression_values.py)  
Combines expression values for a specific reaction across multiple samples and tasks for comparative analysis. Allows for the observation of changes in gene expression of one reaction in different conditions/samples.  

- [`creating_csv_for_escher.py`](./code/creating_csv_for_escher.py)  
Converts excel based expression data into a simplified .CSV format compatible with Escher for mapping of pathways. 

## Main files for running algorithm 
- [`preprocessing.py`](./algorithm/Main_files/preprocessing.py)  
Carries out preprocessing of the data in preparation for [route_optimization.py](./algorithm/Main_files/route_optimization.py) which maps the gene expression data onto the GEM. 
This file sets the model, transcriptomic dataset, method for processing, and location of the storage of results for each run. 

- [`route_optimization.py`](./algorithm/Main_files/route_optimization.py)  
This file takes the preprocessed files from the combination folder created by [preprocessing.py](./algorithm/Main_files/preprocessing.py), and maps the gene expression dataset onto the reactions in the model and predicts the optimal pathway which can be taken for the tasks in the tasklist file. 
