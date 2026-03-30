import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def plotTimeVsConvergence(lst_times, lst_gap):
    # Get gap and total time
    plt.scatter(lst_times, lst_gap)
    plt.yscale("log")
    plt.xlabel('Time (in seconds)')
    plt.ylabel('Convergence gap')
    i = 0
    for x,y in zip(lst_times,lst_gap):
        plt.annotate(lst_tasks[i],
                    (x,y),
                    xytext=(0,0),
                    size = 5,  
                    ha='center')
        i += 1
    plt.savefig("TimeVsConvergence.jpg") 
    plt.show()

def considerNaNs(lst_tasks, lst_reactions, allReactions=False, length=-0.1):
    taskReactionDic = {}
    for i in range(len(lst_tasks)):
        taskReactionDic[lst_tasks[i]] = lst_reactions[i]

    allTasks = []
    if allReactions:
        for el in list(range(349)):
           allTasks.append(el+1)
    else:
        # find highest number
        # construct list of that size
        highest = max(lst_tasks)
        for el in list(range(highest)):
            allTasks.append(el+1)
    # check if in the initial list
    
    lst_ReactionsNaN = []
    for i in allTasks:   
        if i in lst_tasks:
            lst_ReactionsNaN.append(taskReactionDic[i])
        else:
            lst_ReactionsNaN.append(length)
    return [allTasks, lst_ReactionsNaN]

def plotTaskVsConvergence(lst_tasks, lst_gap):
    result_list = [i for _, i in sorted(zip(lst_tasks, lst_gap))]
    lst_tasks = sorted(lst_tasks)
    allTasks, lst_gapNaN = considerNaNs(lst_tasks,result_list,False)
    
    plt.bar(allTasks, lst_gapNaN)
    plt.yscale("log")
    plt.xlabel('Task Number')
    plt.ylabel('Convergence gap')
    plt.show()
    
def plotTaskVsInitialReactions(lst_tasks, lst_reactions):
    result_list = [i for _, i in sorted(zip(lst_tasks, lst_reactions))]
    lst_tasks = sorted(lst_tasks)
    plt.bar(lst_tasks, result_list)
    plt.xlabel('Task Number')
    plt.ylabel('Number of reactions in initial model')
    plt.savefig("TaskVsInitialReactions.jpg") 
    plt.show()    

def plotTaskVsInitialndiMATReactions(lst_tasks, lst_initialReactions, lst_finalReactions):
    # Adapt so that there are 2 bars per task to compare initial and final per task
    # Or not use bars but scatter with subplot and 2 different colors
    # deal with the NANs
    result_list = [i for _, i in sorted(zip(lst_tasks, lst_initialReactions))]
    result_list2 = [i for _, i in sorted(zip(lst_tasks, lst_finalReactions))]
    lst_tasks = sorted(lst_tasks)

    ##
    allTasks, lst_initialReactionsNaN = considerNaNs(lst_tasks, result_list,False, -5)
    _, lst_iMATReactionsNaN = considerNaNs(lst_tasks, result_list2,False, -5)
    
    plt.bar(allTasks, lst_initialReactionsNaN, label="Reactions in final model")
    plt.bar(allTasks, lst_iMATReactionsNaN, label="Reactions after iMAT")
        
    plt.xlabel('Task Number')
    plt.ylabel('Number of reactions in final model')
    plt.legend(loc="upper left")
    #plt.savefig("TaskVsExpNdAllReactions.jpg")
    # plt.xlabel('Task Number')
    # plt.ylabel('Number of reactions in initial model')
    #plt.savefig("TaskVsInitialReactions.jpg") 
    plt.show()

def plotIterationVsConvergence(lst_it, lst_gap):
    plt.scatter(lst_it, lst_gap)
    plt.xlabel('Number of search iterations')
    plt.ylabel('Convergence gap')
    plt.savefig("IterationVsConvergence.jpg") 
    plt.show()
    
def plotNumberReactionsVsConvergence(lst_gap,  lst_reactions):
    plt.scatter(lst_gap, lst_reactions)
    plt.ylabel('Convergence gap')
    plt.xlabel('Number of reactions in initial model')
    i = 0
    for x,y in zip(lst_gap,lst_reactions):
        plt.annotate(lst_tasks[i],
                    (x,y),
                    xytext=(0,0), 
                    ha='center')
        i += 1
    plt.savefig("InitialReactionVsConvergence.jpg") 
    plt.show()
    
def plotNumberReactionsVsTime(lst_reactions, lst_times):
    plt.scatter(lst_reactions, lst_times)
    plt.ylabel('Time (in seconds)')
    plt.xlabel('Number of reactions in initial model')
    plt.savefig("InitialReactionVsTime.jpg") 
    plt.show()
    
def plotNumberInitialReactionsVsReactionsiMAt(lst_InitialReactions, lst_ReactionsiMAT):
    plt.scatter(lst_InitialReactions, lst_ReactionsiMAT)
    plt.ylabel('Number of reactions in final model (after iMAt)')
    plt.xlabel('Number of reactions in initial model')
    plt.savefig("InitialReactionVsiMATReactions.jpg") 
    plt.show()

##
os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/MILPs")

mat_files_directory = os.getcwd() 
print("Python file running")
print(mat_files_directory)

lst_times = []
lst_gap = []
lstTasks = []
lstNumberIteration =[]
lstNumberReactionsInitial = []
lstNumberReactionsFinal = []
lstNumberReactionsiMAT = []
lstNumberExpReactionsFinal = []
lstNumberExpReactionsiMAT = []

lst_kAll = []
lst_colorsAll = []
lst_lastIterationTime = []


task_list = [tasks for tasks in os.listdir(mat_files_directory) if tasks.startswith('Task') ] #or tasks.startswith('Testing')
for el in task_list:
    task_directory = mat_files_directory + "/" + el
    os.chdir(task_directory)
    print(el)
    a = re.findall(r'Task\d+', el)
    lstTasks.append(int(a[0].strip("Task"))) # !!
    # Get which tasks and Sample??
    model = el.split("_")[0]
    mydic = {}
    mydic2 = {}
    # Get the tested parameter + value
    # Move in the file
    dir_list = [dirs for dirs in os.listdir(task_directory) if dirs.startswith('RunNewKApproxHTv2_home') ] #
    for tests in dir_list:
        myTestedParameter = tests.split("_")[1]
        os.chdir(task_directory + "/" + tests)
        # Open file myLog
        # Retrieve important info
        # K infeasible, k feasible, gap
        # K infeasible= 4.468333, K feasible = 4.439662 and last Ktest= 4.441096 
        # Gap = 0.028671 
        mylist = {}
        mylist2 = {}
        i = 0
        j = 0
        notTesting = True
        if notTesting:
            try:
                with open("myLog.txt", "r") as f:
                    content = f.readlines()
                    for line in content:
                        if line.startswith("K infeasible"):
                            a = re.split(", | and ", line)
                            k_infeas = a[0].strip("K infeasible= ")
                            k_feas = a[1].strip("K feasible = ")
                            mylist["k bracket [feasible -- infeasible]"] = [float(k_feas), float(k_infeas)]
                            # 
                        elif line.startswith("Gap"):
                            gap = line.strip("Gap = ").strip()
                            mylist["Gap"] = float(gap)
                            lst_gap.append(float(gap))
                            # 
                        elif line.startswith('Task Score'):
                            score = line.strip("Task Score is ").strip()
                            score = float(score)
                            mylist["Task Score"] = score
                            #"Task Score"
                        elif line.startswith("Amount of reactions"):
                            r = line.strip("Amount of reactions in original and pruned model: ").strip()
                            allreactionsBoth = r
                            allReactions = r.split(" -- ")[0].strip()
                            allReactionsIMAT = r.split(" -- ")[1].strip()
                            #
                            lstNumberReactionsFinal.append(int(allReactions))
                            lstNumberReactionsiMAT.append(int(allReactionsIMAT))
                            mylist["Amount of reactions [Initial -- IMAT]"] = [int(allReactions), int(allReactionsIMAT)]
                        elif line.startswith("Amount of Expression"):
                            r = line.strip("Amount of Expression carrying reactions in original and pruned model: ").strip()
                            expReactionsBoth = r 
                            expReactions = r.split(" -- ")[0].strip()
                            expReactionsIMAT = r.split(" -- ")[1].strip()
                            # 
                            lstNumberExpReactionsFinal.append(int(expReactions))
                            lstNumberExpReactionsiMAT.append(int(expReactionsIMAT))

                            mylist["Amount of reactions w/ expression [Initial -- IMAT]"] = [int(expReactions), int(expReactionsIMAT)]
                        elif line.startswith("Total Time "):
                            time = line.strip("Total Time required is ")
                            time = time.strip()
                            mylist["Total time"] = float(time)
                            lst_times.append(float(time))

                        elif line.startswith("Feasible Space: "):
                            i += 1
                            k = line.split(" = ")[1].strip()
                            m = "k = " + k
                            name = "Feasible space {}".format(i)
                            mylist2[name] = m
                            
                        elif line.startswith("Infeasible Space: "):
                            j += 1
                            k = line.split(" = ")[1].strip()
                            m = "k = " + k
                            name = "Infeasible space {}".format(j)
                            mylist2[name] = m
                        elif line.startswith("Number of reactions in FastCC Irreversible Model: "):
                            numberReactions = int(line.strip("Number of reactions in FastCC Irreversible Model: "))
                            lstNumberReactionsInitial.append(numberReactions)
                    lstNumberIteration.append(i+j)

            except FileExistsError:
                print("Error: No myLog.txt")
            except FileNotFoundError:
                print("Error: No myLog.txt")
            
            mydic[myTestedParameter] = mylist
            mydic2[myTestedParameter] = mylist2  
         
        ## GEt data from LogFile.txt
        # Move to logFiles
        os.chdir(task_directory + "/" + tests + "/logFiles")
        # Open it
        try: 
            with open("logFile.txt", "r") as f:
                lines = f.read()
                eventPlot = True
                if eventPlot:
                    a = re.findall(r'{\[ ?([0-9]+.[0-9]+)\]}', lines)
                    lst_a = a
                    # for feas space, if tested value
                    a = sorted(a)
                    #a = list(set(a))
                    lst_k = []
                    lst_colors = []
                    for el in a:
                        # initial always infeasible??
                        nbr = float(el)
                        lst_k.append(nbr)
                        if nbr > score:
                            lst_colors.append("red")
                        else:
                            lst_colors.append("green")
                    lst_colors = [tuple(lst_colors)]
                    
                    lst_kAll.append(lst_k)
                    lst_colorsAll.append(lst_colors)
                    
                    #plt.scatter(lst_k, list(range(len(lst_k))))
                    #plt.show()
                    
                    # Possible recenter the x axis around the k score to set log
                    
                iterationPlot = False
                if iterationPlot:
                    a = re.findall(r'in ([0-9]+.[0-9]+) seconds', lines)
                    lst_iterationTime = []
                    lst_iterationValue = []
                    lst_col = []
                    i = 0
                    for el in a:
                        if i == len(a)-1:
                            break
                        if i == len(a) -2:
                            lst_lastIterationTime.append(float(el))
                            
                        lst_iterationTime.append(float(el))
                        lst_iterationValue.append(float(lst_a[2+3*i]))
                        
                        i += 1
                    
                    # plt.scatter(lst_iterationTime, lst_iterationValue)
                    plt.axhline(y=score, color='r', linestyle='-')
                    # plt.ylabel("K value")
                    # plt.xlabel("Time (in seconds)")
                    # plt.show()
                    plt.scatter(list(range(len(lst_iterationValue))), lst_iterationValue)
                    lst_timeLimit = [18, 21, 22, 25, 26]
                    for el in lst_timeLimit:
                         plt.annotate("TL", (el, lst_iterationValue[el]+0.3), ha='center')
                        
                    plt.show()
                    print(lst_iterationTime)
                    
                    # consider all tested
                    # not if it passed

        except FileExistsError:
            print("Error: No LogFile.txt")
        except FileNotFoundError:
            print("Error: No LogFile.txt") 
    
    if notTesting:             
        df1 = pd.DataFrame.from_dict(mydic)
        df2 = pd.DataFrame.from_dict(mydic2)
        
        df = pd.concat([df1, df2], sort=False)
        if not df.empty:
            os.chdir(mat_files_directory)


            with pd.ExcelWriter('output.xlsx', mode='a', if_sheet_exists="replace", engine="openpyxl") as writer:  
                df.to_excel(writer, sheet_name=model)
        
    # Addition of the Details automatically
    # DEal with the fact that they'd be on top

# Trying to have all parameters on the same excel
# Distinguish between parameters - create different sheets
os.chdir("D:/Documents/Unif/MSP/BTR_MaCSBio/MaCSBio-GEM/General/Functions/Metabolic_tasks/Jelle/JuTest/Comparisons")

#plotTimeVsConvergence(lst_times, lst_gap)

#plotTaskVsConvergence(lst_tasks, lst_gap)
#plotIterationVsConvergence(lstNumberIteration, lst_gap)
#plotNumberReactionsVsConvergence(lstNumberReactions, lst_gap)
#plotNumberReactionsVsTime(lstNumberReactions, lst_times)
#plotTaskVsConvergence(lst_tasks, lstNumberReactions)
#plotNumberInitialReactionsVsReactionsiMAt(lstNumberReactionsInitial, lstNumberReactionsiMAT)
#plotTaskVsInitialReactions(lst_tasks, lstNumberReactionsInitial)
#plotTaskVsInitialndiMATReactions(lst_tasks, lstNumberReactionsInitial, lstNumberReactionsiMAT)
#plotTaskVsConvergence(lst_tasks, lst_gap)

#plotTaskVsInitialndiMATReactions(lst_tasks, lstNumberReactionsiMAT, lstNumberExpReactionsiMAT)
#plotTaskVsInitialndiMATReactions(lst_tasks, lstNumberReactionsFinal, lstNumberReactionsiMAT)

#plotTaskVsInitialndiMATReactions(lst_tasks, lstNumberReactionsFinal, lstNumberReactionsiMAT)
#plotTaskVsInitialndiMATReactions(lst_tasks, lstNumberReactionsFinal, lstNumberExpReactionsFinal)
plotTaskVsConvergence(lstTasks, lst_gap)

timeVsInitialReactions = False
if timeVsInitialReactions:
    lst_c = []
    for el in lst_gap:
        if el <= 0.1:
            lst_c.append("green")
        else:
            lst_c.append("red")
    plt.scatter(lst_times, lstNumberReactionsInitial, color=lst_c)
    i = 0
    for x,y in zip(lst_times,lstNumberReactionsInitial):
        plt.annotate(str(lstTasks[i]),
                    (x,y))
        i += 1
    plt.xlabel('Total time (in seconds)')
    plt.ylabel('Number of reactions in initial model')
    plt.show() 
ConvergenceVsNumberFinalReactions = False
if ConvergenceVsNumberFinalReactions:
    plt.scatter(lstNumberReactionsiMAT, lst_gap)
    plt.ylabel('Convergence')
    plt.xlabel('Number of reactions in final model')
    i = 0
    for x,y in zip(lstNumberReactionsiMAT,lst_gap):
            plt.annotate(str(int(lst_lastIterationTime[i])),
                        (x,y))
            i += 1
    plt.show()        



# result_list = [i for _, i in sorted(zip(lst_tasks, lst_times))]
# result_list2 = [i for _, i in sorted(zip(lst_tasks, lst_gap))]
# lst_tasks = sorted(lst_tasks)

# x = np.arange(len(result_list))
# width = 0.4

# fig, ax1 = plt.subplots()
# color = 'tab:blue'
# ax1.set_yscale("log")
# ax1.set_xlabel('Task Number')
# ax1.set_ylabel('Total time (in seconds)', color=color)
# ax1.bar(x, result_list, width=width, color=color)
# ax1.tick_params(axis='y', labelcolor=color)
# ax1.set_xticks(x + width/2,lst_tasks) 

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:orange'
# ax2.set_ylabel('Convergence', color=color)  # we already handled the x-label with ax1
# ax2.set_yscale("log")
# ax2.bar(x+width, result_list2, width=width, color=color)
# ax2.tick_params(axis='y', labelcolor=color)

# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()

