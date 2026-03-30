import pandas as pd
import os
import glob

# change folder path
folder_path = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\results\starting with cytosolic APAP"

output_folder = os.path.join(folder_path, "csv_outputs")
os.makedirs(output_folder, exist_ok=True)

# Loop through all Excel files in folder starting with "filtered_Task_"
for file_path in glob.glob(os.path.join(folder_path, "filtered_Task_*_Sample*.xlsx")):
    df = pd.read_excel(file_path)
    # Select only the "ID" and "Expression" columns in the file
    if "ID" in df.columns and "Expression" in df.columns:
        subset = df[["ID", "Expression"]]
        # Reorder columns to make sure "ID" is before "Expression"
        subset = subset[["ID", "Expression"]] 

        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_file = os.path.join(output_folder, f"{base_name}_for_escher.csv")
        subset.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")
    else:
        print(f"Skipping {file_path}, missing required columns.")



