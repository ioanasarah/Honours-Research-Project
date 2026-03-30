data = load("C:\Ioana\_uni\honours\slp\algorithm\Main_files\models\consensus_2024\model_Human-GEM_COBRA_version_17.mat");
model = data.model;
writeCbModel(model, 'format', 'sbml', 'fileName', 'model.xml');
%saveas(model,"C:\Ioana\_uni\honours\slp\algorithm\Main_files\models\consensus_2024\model_matlab.")
