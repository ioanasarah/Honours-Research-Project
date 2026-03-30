# Used to compare the presence of reactions between several different solution for a certain task

 
# Adapt the previously built comparison of the solution to include more than 2
# How to compare?
# Presence and absence of the reaction
import pandas as pd
import numpy as np
import os
import sys
import re
import matplotlib.pyplot as plt
import matplotlib as mpl

# Find the list of all the reactions present for the particular task (in all the solutions)
# With index of how many times
# If index = total number of solutions, do not include it in
# How to select all the run for a task??
print("Python file running")

tasks = ["Task20Sample1"] #"TestingT1S1a=0.5b=0", "Task1Sample1", "Task11Sample1","Task1Sample1TestingFastCC",

inter_tasks = ["Task108Sample1TestingFastCC", "Task1Sample1TestingFastCC"]
os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/MILPs")
for task in tasks:
    print(task)
    allReactions = False
    print(os.getcwd())
    mat_files_directory = os.getcwd() + "\\" + task

    list_tests = [dirs for dirs in os.listdir(mat_files_directory) if dirs.startswith('RunNewK')]
    AllIDS = {}
    for test in list_tests:
        print(test)
        myPath = mat_files_directory + "\\" + test + "\\"
        # Look for ActiveReactionsTask_IMAT_Min_ConvergedSol.xlsx
        try:
            print(myPath)
            try:
                a = pd.read_excel(myPath + "ActiveReactionsTask_IMAT_Min_ConvergedSol.xlsx",na_values=np.nan,header=None, engine="openpyxl")
            except:
                a = pd.read_excel(myPath + "ActiveReactionsExcelFiles\\ActiveReactionsTask_IMAT_Min_ConvergedSol.xlsx",na_values=np.nan,header=None, engine="openpyxl")
        except:
            try:
                a = pd.read_excel(myPath + "ActiveReactionsTask_ConvergedSolution.xlsx",na_values=np.nan,header=None, engine="openpyxl")
            except:
                a = pd.read_excel(myPath + "ActiveReactionsExcelFiles\\ActiveReactionsTask_ConvergedSolution.xlsx",na_values=np.nan,header=None,engine="openpyxl")
            print("Warning: IMAT not used for", test)
        for el in a.loc[:,3]:
            if el not in AllIDS:
                AllIDS[el] = [test]            
            else:
                i = AllIDS.get(el)
                i.append(test)
                AllIDS[el] = i

    list_IDS = AllIDS.keys()    

    if not allReactions:
        temp_lst = []
        for el in list_IDS:
            if len(AllIDS[el]) == len(list_tests):
                # Contained in all, can be removed
                pass    
            else:
                temp_lst.append(el)
        list_IDS = temp_lst        
        print(list_IDS)    
    matrix = []
    for id in list_IDS:
        lst = []
        for el in list_tests:
            if el in AllIDS[id]:
                lst.append(1) 
            else:
                lst.append(0)
        matrix.append(lst)

    mymatrix = np.array(matrix)
    mymatrix = mymatrix.transpose()

    fig, ax = plt.subplots()
    im = ax.imshow(mymatrix)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(list_IDS)), labels=list_IDS, fontsize=5)
    ax.set_yticks(np.arange(len(list_tests)), labels=list_tests, fontsize=5)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(list_tests)):
        for j in range(len(list_IDS)):
            text = ax.text(j, i, mymatrix[i, j],
                        ha="center", va="center", color="w")

    name = task[:15] + ".jpg"
    #os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/Comparisons")
    #plt.savefig(name, bbox_inches="tight")
    plt.show()

# Per solution, create a list with presence and absence of reactions

# Build a graph out of it 
# Plot it, plt

# Notes: Gathering of similar reactions??
