# used to compare the reactions of two solutions (integrating the flux in the comparison) 

import pandas as pd
import numpy as np
import os
import math
# Add an extra column to sort it out - 1, Intact, 2, changed, 3, deleted, 4 added


# From 2 excels to a common excel with the similarities and differences
os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/Comparisons/pathwayComparison")

t = pd.read_excel("Task1.0_LP_ORmedian.xlsx",na_values=np.nan,header=None) #Task1.5_MILP_ORlow_TRORoneOverMedianExp.xlsx
c = pd.read_excel("Task1.5_MILP_ORlow_TRORoneOverMedianExp.xlsx",na_values=np.nan,header=None)

flux = t.loc[:,1]
IDs = t.loc[:,3]
dic_t ={}

for i in range(len(t.loc[:,1])):
    dic_t[t.loc[i,3].strip("''")] = [None, t.loc[i,1] , t.loc[i,2], t.loc[i,4], t.loc[i,5]]
    
dic_c = {}
for i in range(len(c.loc[:,1])):
    dic_c[c.loc[i,3].strip("''")] = [None, c.loc[i,1], c.loc[i,2], c.loc[i,4], c.loc[i,5]]

# Find the differences between the excels
lst_changed = []
lst_added = []
lst_removed = []
for el in dic_t:
    try:
        #mult = 10**20
        #if math.ceil(dic_t[el][1] * mult)//mult != math.ceil(dic_c[el][1] * mult)//mult:
        if dic_t[el][1] != dic_c[el][1]:
            dic_c[el][0] = "Changed"
            dic_t[el][0] = "Changed"
            lst_changed.append(el)
        else:
            dic_c[el][0] = "Intact"
            dic_t[el][0] = "Intact"            
    except:
        dic_t[el][0] = "Removed"
        lst_removed.append(el)

for el in dic_c:
    if dic_c[el][0] == None:
        dic_c[el][0] = "Added"
        lst_added.append(el)
        

# Convert it back to excel, with modification pointed out
# Create a dataframe with all the reactions in order
# 1. Intact
# 2. Changed
# 3. Added
# 4. Removed

# Apply the corresponding color scheme

df = pd.DataFrame(columns=["State", "Flux", "changedFlux", "Expression", "ID", "Name", "Equation"])
for el in dic_t:
    if dic_t[el][0] == "Intact":
        state = 0
        changedFlux = dic_t[el][1]
    elif dic_t[el][0] == "Changed":
        state = 1
        changedFlux = dic_c[el][1]
    elif dic_t[el][0] == "Removed":
        state = 2
        changedFlux = 0
    
    entry = pd.DataFrame.from_dict({
        "State": [state],
        "Flux": [dic_t[el][1]],
        "changedFlux": [changedFlux],
        "Expression":  [dic_t[el][2]],
        "ID": [el],
        "Name": [dic_t[el][3]],
        "Equation": [dic_t[el][4]]
    })    
    df = pd.concat([df, entry], ignore_index=True)

for el in dic_c:
    if dic_c[el][0] == "Added":
        entry = pd.DataFrame.from_dict({
        "Flux": [0],
        "changedFlux": [dic_c[el][1]],
        "Expression":  [dic_c[el][2]],
        "ID": [el],
        "Name": [dic_c[el][3]],
        "Equation": [dic_c[el][4]]
        })    
        df = pd.concat([df, entry], ignore_index=True)
        
def coloredFlux(val):
    lst_color = []
    for el in range(len(df)):
        lst_color.append('color: black')
    for c in lst_changed:
        lst_color[list(df.loc[:, "ID"]).index(c)] = 'color: orange' 
    for el in lst_added:
        lst_color[list(df.loc[:, "ID"]).index(el)] = 'color: green' 
        
    for el in lst_removed:
        lst_color[list(df.loc[:, "ID"]).index(el)] = 'color: red' 
        
    #print(lst_color)
    return lst_color

df.style.apply(coloredFlux).to_excel('myComparisonFlux.xlsx')
