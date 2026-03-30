import os

main_folder = "C:\Ioana\_uni\honours\slp\algorithm\Main_files"
model_folder = os.path.join("models", "Human-Gem")
print("Looking in:", model_folder)
model_files = [file for file in os.listdir(model_folder) if file.lower().startswith("model")]
print("Model files found:", model_files)

if not os.path.exists(model_folder):
    print("Folder does NOT exist!")
else:
    files = os.listdir(model_folder)
    print("Files found:", files)