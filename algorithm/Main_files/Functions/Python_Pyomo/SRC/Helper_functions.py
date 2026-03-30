import pyomo.environ as pe
from cobra.util import remove_cons_vars_from_problem
import copy
import os
import glob
import re
import scipy
import cobra
from cobra import Model, Reaction, Metabolite
import pandas as pd
import warnings
import logging
import numpy as np
import pyomo.opt as po
import openpyxl
import time
import datetime
import pickle
import json
import hdf5storage as hdf5
from MetabolismModeler import file_management as fm
from Legacy import cModels as cM
import cProfile
import gurobipy
from openpyxl.utils import get_column_letter


### create function for reversible
# create function for taking different runs and summarizing them
# create function for reading a json test file
# create statistics for different runs
# create loop law for cobrapy
### adapt names to include digits 001
# start names with
### set filename to include testrun
#### make excel file width of colums correct
# make it so no need to constantly clcik save on excel
# create function for comparing reactions in 2 different solutions
# check in matlab file what expression was used > median of reversible that isn't 0
# calculate score
# fastcore + exchange reactions
# create filtered model from output and make json file
# send new json code to kayle
# download and set up base model for download, get newest version of taskStruct and run LT on it (matlab) then utilize this to create all reversible filtered models

#
#
# print("Fastcc")
# time_start = time.time()
# task_specific_model.solver = "gurobi"
# fast_cc_model = cobra.flux_analysis.fastcc(task_specific_model)
# cobra.io.save_matlab_model(fast_cc_model, main_folder + "\\" + f"Task_{task_number:03}" + f"\\fastcc_model_task{task_number}.mat")
# print(f"Time to run fastcc: {time.time() - time_start} seconds")
