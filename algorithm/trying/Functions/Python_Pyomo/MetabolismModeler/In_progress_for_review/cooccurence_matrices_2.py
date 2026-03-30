import json
import os
import re
import cobra
from cobra import Model, Reaction, Metabolite
import numpy as np
from cobra.io import load_json_model, load_matlab_model

def create_solutions_of_other_samples_in_solution_matrix(
		solution_sets: dict[int, list[str]],
		expression_matrix: [list[float]],  # size 2d reaction * sample
		cobra_model: cobra.Model,  # contains reactions in same order as expression_matrix
		objective_function_per_sample: str,
		expressionless_reaction_value_per_sample: dict[int, float],
):
	#objectivce_values_dict = {}
	#print("went in def")
	union_of_reactions_dict = {}
	reaction_indices = {}
	sorted_reaction_indices = {}
	for idx, sample_number in enumerate(solution_sets.keys()):
		list_of_reactions = solution_sets[sample_number]
		#print(list_of_reactions)
		for reaction_idx, reaction in enumerate(list_of_reactions):
			list_of_reactions = [
				re.sub(r"_f", "", reaction) if re.search(r"_f", reaction) else 
				re.sub(r"_b", "", reaction) if re.search(r"_b", reaction) else 
				re.sub(r"_r", "", reaction) if re.search(r"_r", reaction) else 
				reaction
				for reaction in list_of_reactions if not re.search(r"temporary_exchange", reaction)
			]
			print(list_of_reactions)


#            if re.search(r"_f", list_of_reactions[reaction_idx]):
#                list_of_reactions[reaction_idx] = re.sub(r"_f", "", list_of_reactions[reaction_idx])
#            elif re.search(r"_b", list_of_reactions[reaction_idx]):
#                list_of_reactions[reaction_idx] = re.sub(r"_b", "", list_of_reactions[reaction_idx])
#            elif re.search(r"_r", list_of_reactions[reaction_idx]):
#                list_of_reactions[reaction_idx] = re.sub(r"_r", "", list_of_reactions[reaction_idx])
#             elif re.search(r"temporary_exchange", list_of_reactions[reaction_idx]):
#            	list_of_reactions.remove(list_of_reactions[reaction_idx])

		list_of_reactions = list(set(list_of_reactions)) # removing duplicates
		ordered_list_of_reactions = [(idx, item) for idx, item in enumerate(list_of_reactions, start=1)]  #assigning indexes for each reaction for each sample
		#print(ordered_list_of_reactions)
		union_of_reactions_dict[sample_number] = ordered_list_of_reactions # dict where keys are sample numbers and items are (idx,reaction) -- for all samples
		#print(union_of_reactions_dict)
		#reaction_indices[sample_number] = [idx for idx, reaction in enumerate(cobra_model."metabolites") if reaction.id in union_of_reactions or "metabolites"."id" in union_of_reactions.values()]  #dict where sample number: [idx of reactions in cobra model]

		list_of_reaction = [val[1] for val in ordered_list_of_reactions]
		for idx, reaction in enumerate(cobra_model.reactions):
			# rxn1
			# rxn2
			if reaction.id in list_of_reaction: # true for rxn2
				for idx, rxn in union_of_reactions_dict[sample_number]: # {rxn2, rxn2090, rxn99239}
					if rxn == reaction.id:
						print(idx)
						#reaction_indices[sample_number].append(idx)		
						#print(reaction_indices)
				#print(cobra_model[reaction][id])
			#	reaction_indices[sample_number] = idx for idx, reaction
		# cobra error dict has no attribute reactions
		#not sure if this part of the code works 
		sorted_reaction_indices[sample_number] = [sorted(reaction_indices) for sample_number in reaction_indices.keys()]
		#print(sorted_reaction_indices)
	#print(union_of_reactions_dict)
 # makign a dict where sample_number:list of indexes of reactions in cobra model(ordered)


 #old code:
	#reaction_indices_tuple = [(reaction, idx) for reaction, idx in zip(union_of_reactions, reaction_indices)]
	#print(len(reaction_indices_tuple))
	# ([rxn1, 5000), (rxn2, 6000), (rxn3, 2304)]
   # reaction_indices_tuple = sorted(reaction_indices_tuple, key=lambda x: x[1])
	# sorting tuple based on index -- we would have indexes in chronological order + reactions in a tuple
	#sorted_reactions_dictionary = {tup[0]: idx for idx, tup in enumerate(reaction_indices_tuple)}
	#print(len(sorted_reactions_dictionary))
	#print(len(union_of_reactions))
	#print(sorted_reactions_dictionary)
	#print(union_of_reactions)



	# print(sorted_reactions_dictionary)
	# {"rxn3": [2304], "rnx1": [5000], ...}
	# only search for relevant reactions?
	# [
	#  [4,5,1],
	#  [5,10,10],
	#  [1,10,10]
	#  ]  # returnlist = sum([1/x for x in [2,3,4]]) list comprehension
	#             # returnlist = []
	#             # for x in [2,3,4]:
	#             #     returnlist.append(1/x)
	#             # returnlist = sum(returnlist)
	
	#             mydict = {idx:f"{x} is the number" for idx, x in enumerate([2,3,4])}
	#             # {0: "2 is the number", 1: "3 is the number", 2: "4 is the number"}

	 #    expression_matrix =[
	 #    [1, 2, 3],
	 #    [4, 5, 6],
	 #    [7, 8, 9],
	 #    [10, 11, 12],
	 # ]

	matrix_with_values = np.zeros((len(solution_sets.keys()), len(solution_sets.keys())))
	for row_sample_number in solution_sets.keys():
		# expression_of_sample = expression_matrix[:, row_sample_number]
		expression_of_sample = [value for sublist_ in expression_matrix for idx,
									 value in enumerate(sublist_) if idx == row_sample_number-1]
		expressionless_reaction_value = expressionless_reaction_value_per_sample[row_sample_number]
		expression_of_sample = np.where(expression_of_sample == -1,
										expressionless_reaction_value,
										expression_of_sample)
		for col_sample_number in solution_sets.keys():
			reactions_in_solution = solution_sets[col_sample_number]
			expression_of_reactions = [expression_of_sample[sorted_reactions_dictionary[reaction]] for reaction in reactions_in_solution if reaction in sorted_reactions_dictionary.keys()]
			temporary_exchange_reactions = [reaction for reaction in reactions_in_solution if reaction.startswith("temporary_exchange")]
			amount_to_add_to_sum = len(temporary_exchange_reactions) * (1/expressionless_reaction_value)
			objective_value = sum([1/expression for expression in expression_of_reactions]) + amount_to_add_to_sum
			matrix_with_values[col_sample_number, row_sample_number] = objective_value

	print(matrix_with_values)



if __name__ == "__main__":

	path = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Base_files"
	expression_file_name = "full_expression_reversible_no_threshold.json"
	expression = json.load(open(os.path.join(path, expression_file_name)))
	cobra_model = load_json_model(os.path.join(path, "modelCellfie.json"))
	objective_string = "INIT_irrev"
	solution_sets = {1: [
		"MAR04388_f",
		"MAR04379",
		"MAR04358",
		"MAR04363_f",
		"MAR04368_f",
		"MAR04391_f",
		"MAR06412",
		"MAR07747",
		"MAR03957",
		"MAR04145",
		"MAR04152_f",
		"MAR04209",
		"MAR04410_f",
		"MAR04589_f",
		"MAR06413",
		"MAR06414_f",
		"MAR06914",
		"MAR06916",
		"MAR06918",
		"MAR06921",
		"MAR08746",
		"MAR04855_f",
		"MAR04898_f",
		"MAR05043",
		"MAR06328_f",
		"MAR07638",
		"MAR01485",
		"MAR20069",
		"temporary_exchange_MAM01965[c]",
		"temporary_exchange_MAM02630[c]",
		"temporary_exchange_MAM02751[c]",
		"temporary_exchange_MAM01285[c]",
		"temporary_exchange_MAM02039[c]",
		"temporary_exchange_MAM01371[c]",
		"temporary_exchange_MAM02040[c]",
		"temporary_exchange_MAM01596[c]",
		"MAR04365_r",
		"MAR04373_r",
		"MAR04375_r",
		"MAR04381_r",
		"MAR04280_r",
		"MAR04002_r",
		"MAR06409_r",
		"MAR04141_r",
		"MAR04408_r",
		"MAR04458_r",
		"MAR04652_r",
		"MAR04922_r"
		], 2: ["MAR04922_r"]}
	expressionless_reaction_dict = {i+1: 1.14 for i in range(len(solution_sets.keys()))}
	
	create_solutions_of_other_samples_in_solution_matrix(
		solution_sets,
		expression,
		cobra_model,
		objective_string,
		expressionless_reaction_dict
	)