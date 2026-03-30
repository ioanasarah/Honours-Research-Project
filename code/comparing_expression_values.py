import pandas as pd
import os
import glob

def compare_expression_values(
        folder_path: str,
        output_folder: str, 
        reaction_id: str):
        # Path to your folder with the Excel files
    expression_df = pd.DataFrame(columns=["Task", "Sample", "Expression"])

    for file_path in glob.glob(os.path.join(folder_path, "filtered_Task_*_Sample*.xlsx")):
        df = pd.read_excel(file_path)

        base_name = os.path.basename(file_path)
        sample_name = base_name.split("_")[3]
        task_number = base_name.split("_")[2]
        if reaction_id in df['ID'].values:
            expression_value = df.loc[df['ID'] == reaction_id, 'Expression'].values[0]
            expression_df = pd.concat([expression_df, pd.DataFrame({"Task": [task_number], "Sample": [sample_name],"Expression": [expression_value]})], ignore_index=True)
    print(expression_df)

    output_file = os.path.join(output_folder, f"reaction_expressions_{reaction_id}.csv")
    expression_df.to_csv(output_file, index=False)
    print(f"Saved: {output_file}")


folder_path = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\results\arachidonic_acid_metabolism"
output_folder = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\results\arachidonic_acid_metabolism\expression_values"
compare_expression_values(
    folder_path,
    output_folder, 
    "MAR00958")