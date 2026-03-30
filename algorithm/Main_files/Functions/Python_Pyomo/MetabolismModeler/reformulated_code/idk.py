import pandas as pd

file  = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\expression_datasets\pfos_slp\data_pfos.csv"
#  turn csv into df
df = pd.read_csv(file)
#  if both second and third columns are 0 then remove the whole column from the df 
# df = df.loc[:, ~((df.iloc[1] == 0) & (df.iloc[2] == 0))]
df = df[~((df.iloc[:, 1] == 0) | (df.iloc[:, 2] == 0))]

print(df.head())
df.to_csv(r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\expression_datasets\pfos_slp\data_pfos_new.csv", index=False)