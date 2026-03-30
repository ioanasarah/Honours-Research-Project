initCobraToolbox(false)

inputFile = 'E:\GitClonesArchives\MaCSBio-GEM\General\Functions\Metabolic_tasks\Jelle\Main Functions\MetabolicTasks_MACSBIO_v0.5.xlsx';
taskStruct=generateTaskStructure_2022_Jelle(inputFile);

modelPath = 'SysBio_COBRA_v1.17_consensus.mat';
model = readCbModel(modelPath);

createIrrevModel = false;
if createIrrevModel
    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
    model = modelIrrev;
end

tic
param.taskStructure = taskStruct;
[taskReport, essentialRxns, taskStructure,usedRxns]=checkMetabolicTasksHumanGEM_2023_2_Jelle(model,param);
toc

expressionValuePath = 'geTMM_MAGNET_BC.csv';
expValues = readtable(expressionValuePath);

data.gene = expValues { :,1};

expValueCell = table2cell(expValues);
subsetCell =  expValueCell(:, 2:end);
expValueMatrix = cell2mat(subsetCell);

data.value = expValueMatrix;
data.parsedGPR = GPRparser(model);
tasksData.taskStructure = taskStruct;
tasksData.essentialRxns = essentialRxns;
param.doParallel =0;

[taskInfos]=createMATFiles(data,tasksData,model,param);
disp("Done")






