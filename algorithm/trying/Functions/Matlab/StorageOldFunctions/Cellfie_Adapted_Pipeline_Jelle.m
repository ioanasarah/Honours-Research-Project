%clearvars()
%initCobraToolbox

inputFile = 'E:\Git\MaCSBio-GEM\General\Data\Metabolic_tasks\MetabolicTasks_MACSBIO_v0.5.xlsx';
taskStruct=generateTaskStructure_2022_Jelle(inputFile);

modelPath = 'E:\Git\MaCSBio-GEM\General\Data\Models\Human-GEM_Cobra_v1.01.mat';
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

[reactionsToKeepPerTaskMILP,reactionsToKeepPerTaskLP,ScorebyTask_binaryMILP,ScorebyTask_binaryLP,scoresPerTaskLP, scoresPerTaskMILP,taskInfos, detailScoring]=CellFie_AverageGeneExpression_2023(data,tasksData,model,param);
%[score, score_binary ,taskInfos, detailScoring]=CellFie_2022(data,tasksData,model,param)







