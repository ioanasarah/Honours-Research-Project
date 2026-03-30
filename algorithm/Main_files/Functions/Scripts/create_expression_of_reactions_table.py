# create reaction in model csv file

import cobra
import json
import csv


main_location = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Base_files\\"
model_location = main_location + "modelCellfie.json"
expression_location = main_location + "full_expression_reversible_no_threshold.json"


model = cobra.io.load_json_model(model_location)
with open(expression_location, "r") as f:
    expression_data = json.load(f)

reactions = [reaction.id for reaction in model.reactions]
samples = list(range(1, 323))

csv_file = "reaction_data_no_thresh.csv"
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the header (sample numbers)
    writer.writerow(["Reaction/Sample"] + samples)
    
    # Write the data row by row (reaction name + values)
    for i, reaction in enumerate(reactions):
        writer.writerow([reaction] + expression_data[i])


json_file = "reaction_data_no_thresh.json"
reaction_data = {}

# Create a JSON structure where each reaction is a dictionary of sample-to-value mappings
for i, reaction in enumerate(reactions):
    reaction_data[reaction] = {f"Sample_{sample}": expression_data[i][j] for j, sample in enumerate(samples)}

# Write the data to a JSON file
with open(json_file, 'w') as file:
    json.dump(reaction_data, file, indent=4)
