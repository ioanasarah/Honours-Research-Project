import cobra
from cobra.io import load_matlab_model, save_json_model
import os


# Directory containing the .mat files
mat_files_directory = 'E:/Git/MaCSBio-GEM/General/Functions/Metabolic_tasks/Task1Sample1FastCCIrrevLoopLaw/'

# Create the "newModels" directory
new_models_directory = os.path.join(mat_files_directory, 'newModels')
os.makedirs(new_models_directory, exist_ok=True)

fluxTables_directory = os.path.join(mat_files_directory, 'fluxTables')
os.makedirs(fluxTables_directory, exist_ok=True)

LastValidSOlutions_directory = os.path.join(mat_files_directory, 'LastValidSolutions')
os.makedirs(LastValidSOlutions_directory, exist_ok=True)

ActiveReactionsExcelFiles_directory = os.path.join(mat_files_directory, 'ActiveReactionsExcelFiles')
os.makedirs(ActiveReactionsExcelFiles_directory, exist_ok=True)

# Find all newModel .mat files
file_list = [file for file in os.listdir(mat_files_directory) if file.startswith('newModel') and file.endswith('.mat')]
# Move the .mat files to the "newModels" directory
for file in file_list:
    old_path = os.path.join(mat_files_directory, file)
    new_path = os.path.join(new_models_directory, file)
    os.rename(old_path, new_path)

# Find all newModel .mat files
file_list = [file for file in os.listdir(mat_files_directory) if file.startswith('fluxTable') and file.endswith('.mat')]
for file in file_list:
    old_path = os.path.join(mat_files_directory, file)
    new_path = os.path.join(fluxTables_directory, file)
    os.rename(old_path, new_path)

# Find all newModel .mat files
file_list = [file for file in os.listdir(mat_files_directory) if file.startswith('lastValidSolution') and file.endswith('.mat')]
for file in file_list:
    old_path = os.path.join(mat_files_directory, file)
    new_path = os.path.join(LastValidSOlutions_directory, file)
    os.rename(old_path, new_path)


file_list = [file for file in os.listdir(mat_files_directory) if file.startswith('ActiveReactionsTask') and file.endswith('.xlsx')]
for file in file_list:
    old_path = os.path.join(mat_files_directory, file)
    new_path = os.path.join(ActiveReactionsExcelFiles_directory, file)
    os.rename(old_path, new_path)



# Create the "ModelsJSON" directory
json_files_directory = os.path.join(mat_files_directory, 'ModelsJSON')
os.makedirs(json_files_directory, exist_ok=True)

# Loop over each file in the "newModels" directory
#print((os.listdir(new_models_directory))
for file in os.listdir(new_models_directory):
    # Extract the name without the extension
    name = os.path.splitext(file)[0]
    print("\n")
    print(file)
    
    # Load the MATLAB model
    model = load_matlab_model(os.path.join(new_models_directory, file))
    
    # Create the JSON file in the "ModelsJSON" directory
    json_file = name + '.json'
    json_path = os.path.join(json_files_directory, json_file)
    save_json_model(model, json_path)



