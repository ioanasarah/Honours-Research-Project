%clearvars()

inputFile = 'E:\Git\MaCSBio-GEM\General\Data\Metabolic_tasks\MetabolicTasks_MACSBIO_v0.5.xlsx';
taskStruct=generateTaskStructure_2022_Jelle(inputFile);

modelPath = 'E:\Git\MaCSBio-GEM\General\Data\Models\Human-GEM_Cobra_v1.01.mat';
model = readCbModel(modelPath);
param.taskStructure = taskStruct;
[taskReport, essentialRxns, taskStructure,usedRxns]=checkMetabolicTasksHumanGEM_2023_2_Jelle(model,param);


expressionValuePath = 'geTMM_MAGNET_BC.csv';
expValues = readtable(expressionValuePath);

%data.gene = model.genes;
data.gene = expValues { :,1};

% %%% set expression data to be same size of model.genes, added 1 values for
% %%% 30 missing genes from expr data
% expColumn = expValues(:, 1).Variables;
% overlapIndices = ismember(expColumn, data.gene);
% filteredExpValues = expValues(overlapIndices, :);
% filteredExpColumn = filteredExpValues(:, 1).Variables;
% extraNames = setdiff(model.genes, filteredExpColumn);
% extraTable = array2table(ones(numel(extraNames), width(filteredExpValues)), 'VariableNames', filteredExpValues.Properties.VariableNames);
% extraTable.(filteredExpValues.Properties.VariableNames{1}) = extraNames;
% finalTable = [filteredExpValues; extraTable];

expValueCell = table2cell(expValues);
subsetCell =  expValueCell(:, 2:end);
expValueMatrix = cell2mat(subsetCell);

data.value = expValueMatrix;
data.parsedGPR = GPRparser(model);
tasksData.taskStructure = taskStruct;
tasksData.essentialRxns = essentialRxns;

[score, score_binary ,taskInfos, detailScoring]=CellFie_AverageGeneExpression_2023(data,tasksData,model);









