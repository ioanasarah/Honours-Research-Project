# Used to compare the code of the functions used for the run

import os 
import difflib

def compare_run(folder1, folder2):
    # Comparison should be
    # Are there the same functions?
    # Are these functions the same
    lst_all_functions1 = [fcts.split('.')[0] for fcts in os.listdir(folder1) if fcts.endswith(".txt")] 
    lst_all_functions2 = [fcts.split('.')[0] for fcts in os.listdir(folder2) if fcts.endswith(".txt")] 
    
    lst_similar = []
    lst_different = []
    lst_notPresentIn2 = []
    lst_notPresentIn1 = []
    for el in lst_all_functions1:
        if el in lst_all_functions2:    
            #check if text = text
            # if not, comparison
            with open(folder1 + '/' + el + '.txt', 'r') as file_1, open(folder2 + '/' + el + '.txt', 'r') as file_2:
                file_1_text = file_1.readlines()
                file_2_text = file_2.readlines()
                if file_1_text == file_2_text:
                    lst_similar.append(el)
                else:
                    lst_different.append(el)
                    diff = difflib.unified_diff(
                        file_1_text, file_2_text, fromfile="file1.txt",
                        tofile="file2.txt", lineterm='')
                    with open('difference '+el +'.txt', "w") as f:
                        for line in diff:
                            f.writelines(line)
        else:
            lst_notPresentIn2.append(el)
            
    for el in lst_all_functions2:
        if el not in lst_all_functions1:
            lst_notPresentIn1.append(el)
            
    
    with open('differenceAll.txt', "w") as f:
        f.write("Functions that are changed :" + str(lst_different) + "\n")
        f.write("Functions that were added :" + str(lst_notPresentIn1) + "\n")
        f.write("Functions that were removed :" + str(lst_notPresentIn2) + "\n")
    return lst_similar, lst_different, lst_notPresentIn1, lst_notPresentIn2


folder1 = ""
folder2 = ""

compare_run(folder1, folder2)