
load('Human-GEM_Cobra_v1.01.mat')
initCobraToolbox

param.taskStructure = generateTaskStructure_2022('MetabolicTasks_MACSBIO_v0.5.xlsx'); %path ? input excel file 
% generateTaskStructure_2022 model and path -> task list as matlab structure ; 
% subset 
% param.taskStructure = param.taskStructure([1, 246]); % put number of tasks here, inside square brackets
param.taskStructure = param.taskStructure([1 : 348]); % put number of tasks here, inside square brackets



%param.inputFile = (" "); % add path to excel file)
param.getEssential = 1; 
param.saveUsed = 1;

[taskReport, essentialRxns, taskStructure,usedRxns]=checkMetabolicTasksHumanGEM_2023_2(model,param);

% writecell(usedRxns{2,1},'usedRxns_Paps.txt');
% for task 273 Synthesis of ribose-5-phosphate:
% writecell(usedRxns{1,1},['GLUCOSE_ATP.txt']);
writecell(usedRxns{1,1},['Used_Rxns_100_200.txt']);

% writematrix is only for logical array 
%-> essential rxn list and used rxns for all the tasks 


