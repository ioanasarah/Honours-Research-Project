# Which functions are used in this script? and in the functions it uses?

# Look through the names of the functions of the folder
# Iterately check if used in function

# List of all functions name
# Main function
#   Open the function
#   Add used function into list (if find, add)
# Go through the used functions
#   Open the function
#   Add used function into list (if find, add)
import os

currentFolder = os.getcwd()

main_function = "CellFie_AverageGeneExpression_2023_clean"
mat_files_directory = "/home/jcornelis/Documents/MaCSBio-GEM-GIT/General/Functions/Metabolic_tasks/Jelle/Main Functions/"

lst_all_functions = [fcts.split('.')[0] for fcts in os.listdir(mat_files_directory) if fcts.endswith(".m")] 

def look_through(fct, lst_mem=[]):
    with open(fct +".m", "r") as f:
        lst_used_functions = set()
        code = f.read()
        code = code[code.index("\n"):]
        for function in lst_all_functions:
            if code.find(function) != -1:
                lst_used_functions.add(function)
        lst_used_functions2 = lst_used_functions.copy()

        if lst_mem == lst_used_functions:
            return lst_used_functions
        for el in lst_used_functions:
            a = look_through(el, lst_used_functions2)
            for i in a:
                lst_used_functions2.add(i)
        return sorted(lst_used_functions2)
  
def get_code(lst_fct):
    # Saves code of all - quite storage intensive?
    # other way? 
    #   -specify some functions which should not :(
    #   -shrink it afterwards? zip, etc
    os.mkdir(currentFolder + "/LoggedFunctions")
    for fct in lst_fct:
        os.chdir(mat_files_directory)
        with open(fct +".m", "r") as f:
            temp_code = f.read()
        os.chdir(currentFolder+ "/LoggedFunctions")
        with open(fct + ".txt", "w") as wr:
            wr.write(temp_code)  
      
test_function = "createCoreMILPProblem"
test_functio2 = "createKMILPProblemAverageExpressionCoreMILP"
test_function3 = "cellfie"

os.chdir(mat_files_directory)
get_code(look_through(main_function))

# os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/Comparisons")
# with open("theFunction.txt", "w") as wr:
#     for el in a:
#         wr.write(el + "\n")   
    
# Comparison should be
# Are there the same functions?
# Are these functions the same

# Removes all spaces, enter not to mess up if space difference
# Compare line by line 
# https://www.geeksforgeeks.org/compare-two-files-line-