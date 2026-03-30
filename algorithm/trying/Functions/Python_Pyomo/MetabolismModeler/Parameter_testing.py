###################################################### 3d parameter plot
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import os
import gurobi_logtools as glt_gurb
from io import StringIO
import json
import math
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import numpy as np
import re
from sklearn.datasets import make_classification
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from cobra.io import save_json_model, load_json_model
from cobra import Model
from os import path, makedirs, listdir
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import seaborn as sns


def check_SV_constraints_and_model_feasilibity(
        cobra_model: Model,
):
    solution = cobra_model.optimize()
    if solution.status == "optimal":
        print("Model is feasible")
        return True
    else:
        print("Model is infeasible")
        return False


def check_all_json_models_in_dir_for_feasibility(
        dir_with_models: str,
) -> Dict[str, bool]:
    is_feasible = {}
    for file in listdir(dir_with_models):
        if file.endswith(".json") and file.startswith("model_json_"):
            try:
                filename = file.split(".json")[0]
                filename = filename.split("model_json_")[1]
                model = load_json_model(path.join(dir_with_models, file))
                is_feasible[filename] = check_SV_constraints_and_model_feasilibity(model)
                print(f"Model {filename} checked")
            except:
                print(f"Model {file} could not be loaded")
                continue
    return is_feasible


from datetime import datetime


def find_closest_log_file(filename):
    files_dir = os.getcwd()
    log_file_name = f"{filename}_gurobi.log"
    filtered_time = os.path.getmtime(log_file_name)

    # Initialize variables to keep track of the closest log file and time difference
    closest_excel_file = None
    smallest_time_diff = float('inf')

    # Iterate over all .log files in the specified directory
    for excel_file in os.listdir(files_dir):
        if excel_file.endswith(".xlsx") and excel_file.startswith("filtered_"):
            excel_file_path = os.path.join(files_dir, excel_file)
            excel_file_time = os.path.getmtime(excel_file_path)

            # Calculate the time difference
            time_diff = abs(excel_file_time - filtered_time)

            # Update the closest log file if this one is closer in time
            if time_diff < smallest_time_diff:
                closest_excel_file = excel_file
                smallest_time_diff = time_diff

    print(closest_excel_file)
    return closest_excel_file


def read_xlsx_to_identify_integrality_violations(filename
                                                 ) -> bool:
    try:
        file = f"filtered_{filename}.xlsx"
        df = pd.read_excel(file)
        print(df.iloc[:, 7].unique())
        print(set(df.iloc[:, 7].unique()).issubset({0, 1}))
        if set(df.iloc[:, 7].unique()).issubset({0, 1}):
            return True
        else:
            return False
    except:
        print(f"could not find {file}")
        filename_2 = find_closest_log_file(filename)
        # .split("_gurobi.log")[0]
        # file_2 = f"filtered_{filename_2}.xlsx"
        try:
            df = pd.read_excel(filename_2)
            print(df.iloc[:, 7].unique())
            print(set(df.iloc[:, 7].unique()).issubset({0, 1}))
            if set(df.iloc[:, 7].unique()).issubset({0, 1}):
                return True
            else:
                return False
        except:
            print(f'Also could not find {filename_2}')
            return None
        return None


def read_xlsx_to_identify_number_of_reactions(filename
                                              ) -> bool:
    try:
        file = f"filtered_{filename}.xlsx"
        df = pd.read_excel(file)
        return len(df.index)
    except:
        print(f"could not find {file}")
        filename_2 = find_closest_log_file(filename)
        try:
            df = pd.read_excel(filename_2)
            return len(df.index)
        except:
            print(f'Also could not find {filename_2}')
            return None
        return None


def run_plotly_dash(path, z_is_column="Runtime", df_name="*"):
    os.chdir(path)
    # print(os.getcwd())
    summary, timelines = glt_gurb.get_dataframe([f"{df_name}.log"], timelines=True)
    df = summary
    # print(df)
    is_feasible_dict = check_all_json_models_in_dir_for_feasibility(path)

    df = df[["Runtime", 'FeasibilityTol (Parameter)', 'IntFeasTol (Parameter)', "Status", "LogFilePath", "ObjVal"]]
    df['isFeasible'] = df['LogFilePath'].apply(lambda x: is_feasible_dict.get(x.split("_gurobi.log")[0], False))
    df['notViolatedIntegrality'] = df['LogFilePath'].apply(lambda x: read_xlsx_to_identify_integrality_violations(x.split("_gurobi.log")[0]))
    df["numberReactions"] = df['LogFilePath'].apply(lambda x: read_xlsx_to_identify_number_of_reactions(x.split("_gurobi.log")[0]))
    print(df)
    # print(df['LogFilePath'])
    # print(is_feasible_dict.keys())
    # print(df.at[5, 'LogFilePath'].split("_gurobi.log")[0])
    # print(df["isFeasible"])
    gradient_column = "Epsilon"
    if z_is_column == "Epsilon":
        gradient_column = "Runtime"

    def extract_epsilon(text):
        # Find all the numbers that follow INIT_irrev
        matches = re.findall(r'(1e-?[0-9]?[0-9])', text)
        # print(matches)
        if len(matches) >= 3:
            # If there are 3 or more occurrences, take the first one
            return matches[0]
        elif len(matches) == 2:
            # If there are 2 occurrences, take the portion before the first 1e-
            match = re.search(r'INIT_irrev([0-9]+(?:\.[0-9]+)?)1e-', text)
            if match:
                return match.group(1)
        return None

    def is_correct(objv):
        if objv <= 600:
            return False
        return True

    def clip_runtime(time):
        if time >= 3000:
            return 3000
        return time

    def apply_log(text, base=10):
        if text is not None:
            text = float(text)
            return math.log(text, base)
        return None

    def get_value(value):
        return float(value)

    df["Epsilon"] = df["LogFilePath"].apply(extract_epsilon)
    df["Correct"] = df["ObjVal"].apply(is_correct)
    df["Runtime"] = df["Runtime"].apply(clip_runtime)
    df["FeasibilityTol (Parameter)"] = df["FeasibilityTol (Parameter)"].apply(apply_log)
    df["IntFeasTol (Parameter)"] = df["IntFeasTol (Parameter)"].apply(apply_log)
    df = df.dropna()
    epsilon = df["Epsilon"].apply(get_value)
    runtime = df["Runtime"].apply(get_value)

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     print(df)

    if gradient_column == "Epsilon":
        gradient_color_col = df["Epsilon"].apply(apply_log)
    else:
        gradient_color_col = df["Runtime"].apply(get_value)
        df["Epsilon"] = df["Epsilon"].apply(apply_log)
    # print(df["Epsilon"].apply(get_value) )
    # print(df["Runtime"].apply(get_value) )
    X = df[['FeasibilityTol (Parameter)', 'IntFeasTol (Parameter)', z_is_column]]
    Y = df['Correct']
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    svc = SVC(kernel='linear')
    svc.fit(X, Y)
    x_vals = np.linspace(df["FeasibilityTol (Parameter)"].min(), df["FeasibilityTol (Parameter)"].max(), 100)
    y_vals = np.linspace(df["IntFeasTol (Parameter)"].min(), df["IntFeasTol (Parameter)"].max(), 100)
    x, y = np.meshgrid(x_vals, y_vals)
    z = lambda x, y: (-svc.intercept_[0] - svc.coef_[0][0] * x - svc.coef_[0][1] * y) / svc.coef_[0][2]

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16, 7))

    # First subplot: Scatter plot of ObjVal by isFeasible
    for label, group in df.groupby('isFeasible'):
        color = 'red' if label == False else 'green'
        axes[0].scatter(group.index, group['ObjVal'], label=f'isFeasible={label}', color=color)
    axes[0].set_xlabel('Index')
    axes[0].set_ylabel('ObjVal')
    axes[0].set_title('ObjVal by isFeasible')
    axes[0].legend()

    # Second subplot: Scatter plot of ObjVal by notViolatedIntegrality
    for label, group in df.groupby('notViolatedIntegrality'):
        color = 'red' if label == False else 'green'
        axes[1].scatter(group.index, group['ObjVal'], label=f'notViolatedIntegrality={label}', color=color)
    axes[1].set_xlabel('Index')
    axes[1].set_ylabel('ObjVal')
    axes[1].set_title('ObjVal by notViolatedIntegrality')
    axes[1].legend()

    df['NoViolation_Feasible'] = df['notViolatedIntegrality'] & df['isFeasible']

    # Second subplot: Scatter plot of ObjVal by notViolatedIntegrality
    for label, group in df.groupby('NoViolation_Feasible'):
        color = 'red' if label == False else 'green'
        axes[2].scatter(group.index, group['ObjVal'], label=f'NoViolation_Feasible={label}', color=color)
    axes[2].set_xlabel('Index')
    axes[2].set_ylabel('ObjVal')
    axes[2].set_title('ObjVal by NoViolation_Feasible')
    axes[2].legend()
    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the combined plot
    plt.show()

    # Create subplots: 1 row and 2 columns
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 7))

    # Filter the DataFrame to only include rows where NoViolation_Feasible is True
    filtered_df = df[df['NoViolation_Feasible'] == True]
    filtered_df = filtered_df.reset_index(drop=True)
    print(filtered_df["Status"].unique())
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(filtered_df[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)", "ObjVal"]].sort_values("Runtime"))
    filtered_only_optimal = filtered_df[filtered_df['Status'] == "OPTIMAL"]
    print(f'Unique objective values = {filtered_only_optimal["ObjVal"].unique()}')

    # Further filtering based on ObjVal
    filtered_only_optimal = filtered_only_optimal[
        (filtered_only_optimal["ObjVal"] >= 606.707) &
        (filtered_only_optimal["ObjVal"] <= 606.80)
        ]
    print(filtered_only_optimal[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)"]])
    filtered_df[["Runtime", "Epsilon", "FeasibilityTol (Parameter)", "IntFeasTol (Parameter)", "ObjVal"]].sort_values("Runtime").to_excel(
        "filtered_data_frame.xlsx")
    # Ensure DataFrame has a reset index
    # Jitter for ObjVal and Runtime
    jitter_strength = 0.0  # Adjust this value to control the amount of jitter
    objval_jitter = filtered_df['ObjVal'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])
    runtime_jitter = filtered_df['Runtime'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])

    # Epsilon values and assign categorical colors
    filtered_df['Epsilon'] = pd.to_numeric(filtered_df['Epsilon'], errors='coerce')
    Epsilon_values = filtered_df['Epsilon'].values
    Epsilon_values_groups = filtered_df['Epsilon'].unique()

    # Generate a color palette for the unique Epsilon values
    palette = sns.color_palette("bright", len(Epsilon_values_groups))  # Use any color palette

    # Map each unique Epsilon to a color
    color_map = {eps: color for eps, color in zip(Epsilon_values_groups, palette)}
    colors = [color_map[eps] for eps in Epsilon_values]

    # Create the scatter plot for the first subplot (axes[0])
    sc = axes[0].scatter(objval_jitter, runtime_jitter, c=colors, label='NoViolation_Correct=True')

    # Create a custom legend for the Epsilon groups
    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', label=f'Epsilon={eps}')
               for eps, color in color_map.items()]
    axes[0].legend(handles=handles, title='Epsilon Groups')

    # Set labels and title for the first subplot (axes[0])
    axes[0].set_xlabel('ObjVal')
    axes[0].set_ylabel('Runtime')
    axes[0].set_title('Scatter Plot of Runtime vs ObjVal (NoViolation_Feasible=True)')

    # Annotate points with their index values and draw lines
    for i in range(len(filtered_df)):
        # Draw a line from the point to the annotation
        # i = 0
        axes[0].annotate(
            i,  # The index value
            (objval_jitter[i], runtime_jitter[i]),  # The coordinates of the point
            textcoords="offset points",  # Specify the position of the annotation text
            xytext=(5, 5),  # Offset from the point (so the text isn't on the point)
            arrowprops=dict(arrowstyle="->", lw=0.5),  # Arrow pointing to the point
            fontsize=8  # Font size for the annotation
        )

    # Jitter for ObjVal and numberReactions (for the second subplot)
    jitter_strength = 0.0  # Adjust this value to control the amount of jitter
    objval_jitter = filtered_df['ObjVal'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])
    number_reactions_jitter = filtered_df['numberReactions'] + np.random.uniform(-jitter_strength, jitter_strength, size=filtered_df.shape[0])

    # Create the scatter plot for the second subplot (axes[1])
    # axes[1].scatter(objval_jitter, number_reactions_jitter, color='green', label='NoViolation_Correct=True')

    # Create the scatter plot for the first subplot (axes[0])
    sc = axes[1].scatter(objval_jitter, number_reactions_jitter, c=colors, label='NoViolation_Correct=True')

    # Create a custom legend for the Epsilon groups
    handles = [plt.Line2D([1], [1], marker='o', color=color, linestyle='', label=f'Epsilon={eps}')
               for eps, color in color_map.items()]
    axes[1].legend(handles=handles, title='Epsilon Groups')

    # Annotate points with their index values and draw lines
    for i in range(len(filtered_df)):
        # Draw a line from the point to the annotation
        axes[1].annotate(
            i,  # The index value
            (objval_jitter[i], number_reactions_jitter[i]),  # The coordinates of the point
            textcoords="offset points",  # Specify the position of the annotation text
            xytext=(5, 5),  # Offset from the point (so the text isn't on the point)
            arrowprops=dict(arrowstyle="->", lw=0.5),  # Arrow pointing to the point
            fontsize=8  # Font size for the annotation
        )

    # Set labels and title for the second subplot (axes[1])
    axes[1].set_xlabel('ObjVal')
    axes[1].set_ylabel('numberReactions')
    axes[1].set_title('Scatter Plot of ObjVal vs numberReactions (NoViolation_Feasible=True)')
    axes[1].legend()

    # Show the plot
    plt.show()

    scatter = go.Scatter3d(
        x=df['FeasibilityTol (Parameter)'],
        y=df['IntFeasTol (Parameter)'],
        z=df[f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=10,  # Increase point size for better visibility
            color=gradient_color_col,
            colorscale=["green", "red"],
            opacity=0.8,
            colorbar=dict(
                title=f'{gradient_column}',
                # tickvals=[0, 0.5, 1],
                # ticktext=['Low', 'Medium', 'High']
            )
        ),
        text=df[f'{gradient_column}'],  # Tooltip shows the Runtime
    )
    highlight = go.Scatter3d(
        x=df.loc[df['Correct'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['Correct'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['Correct'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=18,  # Increase the size of the red circles
            color='red',
            symbol='circle-open',
            line=dict(width=3, color='red'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight2 = go.Scatter3d(
        x=df.loc[df['isFeasible'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['isFeasible'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['isFeasible'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=22,  # Increase the size of the red circles
            color='green',
            symbol='circle-open',
            line=dict(width=3, color='green'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight3 = go.Scatter3d(
        x=df.loc[df['notViolatedIntegrality'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['notViolatedIntegrality'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['notViolatedIntegrality'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=26,  # Increase the size of the red circles
            color='yellow',
            symbol='circle-open',
            line=dict(width=3, color='yellow'),  # Thicker red outline
        ),
        showlegend=False
    )
    highlight4 = go.Scatter3d(
        x=df.loc[df['NoViolation_Feasible'], 'FeasibilityTol (Parameter)'],
        y=df.loc[df['NoViolation_Feasible'], 'IntFeasTol (Parameter)'],
        z=df.loc[df['NoViolation_Feasible'], f'{z_is_column}'],
        mode='markers',
        marker=dict(
            size=30,  # Increase the size of the red circles
            color='Purple',
            symbol='circle-open',
            line=dict(width=4, color='Purple'),  # Thicker red outline
        ),
        showlegend=False
    )

    fig = go.Figure(data=[scatter,
                          # highlight,
                          # highlight2,
                          # highlight3,
                          highlight4,
                          ])
    fig.add_surface(x=x, y=y, z=z(x, y), colorscale='Greys', showscale=False)
    min_z = min(df[f'{z_is_column}'])
    max_z = max(df[f'{z_is_column}'])
    if gradient_column == "Runtime":
        fig.update_scenes(zaxis_range=[min_z, max_z])
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title='FeasibilityTol (Parameter)',
                # type='log',  # Log scale for x-axis
            ),
            yaxis=dict(
                title='IntFeasTol (Parameter)',
                # type='log',  # Log scale for y-axis
            ),
            zaxis=dict(
                title=f'{z_is_column}',
                # type='log',  # Log scale for z-axis
            ),
        ),
    )

    # Show the figure
    fig.show()
    fig.write_html(f"{z_is_column}_to_z.html")


if __name__ == '__main__':
    # main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_NT_LT_INIT_LINEAR_MIN_Tarray_3samples_8tasks"
    main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_no_threshold_fractional_options_test_enum"
    # if not os.getcwd().endswith("results"):
    #     os.chdir(os.path.join(os.getcwd(), "results"))
    run_plotly_dash(main_folder)

################################################################### example for troubleshooting some code

##### Some stuff that isn't working

# A = 5
# # B = 6

# def some_function(x):
#     print("this is running now")
#     return str(x) + "hello"

# C = some_function('d')


# if __name__ == "__main__":
#     D = 5
#

#### Should provide a list with integers
#    so we just test this part of the code by providing a list with integers

# list_with_integers = [2,43,5,12,7,4]

# def function_to_do_things_with_integers(integer_list):
#     for idx in range(len(integer_list)):
#         print(idx)

# function_to_do_things_with_integers(list_with_integers)

