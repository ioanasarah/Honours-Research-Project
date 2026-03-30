import cobra
filename = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\models\xml_file\model.xml" # path to xml model
# cobra.io.sbml.validate_sbml_model(
#     filename, 
#     check_model=True)-> Tuple[cobra.Model, dict]

cobra.io.sbml.validate_sbml_model(
    filename)-> Tuple[Optional[cobra.core.Model], dict]