% Code from CellFie_AverageGeneExpression_2023_clean.m
% Adapted to work for Bachelor Thesis Research by Kayle Boessen

% .. Authors:
%     - Anne Richelle, January 2019
%	  - Modified by Eva (MACSBIO intern) to work on iHuman model, date unclear
%	  - Further modified by Bastien Nihant for flexibility and adapting to new cobraToolbox behaviors, November 2022
%     - Jelle Bonthuis
%     - Justin Cornélis
%     - Kayle Boessen

%% Clear workspace and cmd window
clear
clc

%% Initialize COBRA Toolbox
currentFolder = pwd;
initCobraToolbox(false);
cd(currentFolder);

%% Load expression reactions, significance, task data, and model

startTime = tic;

load("expressionRxnsCellfie")
load("significanceCellfie")
load("tasksDataCellfie")
load("modelCellfie")

expValues = expressionRxns;
TaskStructure = TasksData.taskStructure;

nSamples = size(expressionRxns,2);

%% Adjustment of expression data

% if no gene is associated with one of the reactions -
% remove the reaction from the count
noGene = find(sum(expValues,2)==-nSamples);
if ~isempty(noGene)
    expValues(noGene,:) = [];
end

noGene = find(isnan(sum(expValues,2)));
if ~isempty(noGene)
    expValues(noGene,:) = -1;
end

% include significance
expValuesAdjusted = expValues .* significance;
expValuesAdjusted(expValues == -1) = -1;

noGene = find(isnan(sum(expValuesAdjusted,2)));
if ~isempty(noGene)
    expValuesAdjusted(noGene,:) = -1;
end

%% Variable settings
epsilon = 0.001;
osenseStr = "min";
irrev = true;
orphanHandling = "median";
mainFolder = pwd;
addpath(mainFolder);


%% Loop

% Create an empty task model
EmptyTaskModel = createEmptyTaskModel(Model);

for iTask = 1:1 %size(TaskStructure,1)  % Choose which task(s) will run

    for jSample = 1:1 %nSamples         % For testing, adjust to 1 or whatever sample number you want

        % Separate expression values for current sample
        sampleExpValues = expValuesAdjusted(:, jSample);

        % Create a task specific model and checks if task passes 
        % (Utilizes code from CheckMetabolicTasks)
        [passTask, TaskSpecificModel] = openExchangeTaskReactionsAndCheckTask(Model,TaskStructure,iTask);

        % Adjust expression data and core model to correct size (based
        % on how many placeholder reactions were removed)
        diffInRxns = length(EmptyTaskModel.rxns) - length(TaskSpecificModel.rxns);
        TaskSpecificModelCore = adjustTaskModel(TaskSpecificModel, diffInRxns);

        sizeDiff = numel(TaskSpecificModel.rxns) - numel(sampleExpValues);
        sampleExpValuesAdjustedOld = [sampleExpValues; repmat(-1, sizeDiff, 1)];

        % % Make irreversible
        % [taskModel, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(TaskSpecificModel);
        % sampleExpValuesAdjustedTask = sampleExpValuesAdjustedOld;
        % 
        % % TODO
        % % is this for loop necessary? it seems like it is unused in further
        % % code... -> Ask Jelle
        % for g=1:size(rev2irrev,1)
        %     if size(rev2irrev{g},2) == 2
        %         sampleExpValuesAdjustedTask(rev2irrev{g}(2)) = sampleExpValuesAdjustedOld(rev2irrev{g}(1));
        %         %sampleExpValuesAdjustedTask(rev2irrev{g,1}(2)) = rev2irrev{g,1}(1);
        %     end
        % end

        if passTask
            
            cd(mainFolder)

            %%% OR adjustment code
            % TODO - Can this be moved to global variables?
            % -> ask Jelle
            alpha = 0.5; % even addition of the OR adjustment as regular 
            % objective funciton at 0.5
            % alpha gives the weights to the OR adjustment of the objective
            % function and to the original objective function by setting
            beta = 0;  % total addition to max score if all orphan 
            % transport reactions would be 0
            lbFromFVA = 1000; %lower bound from Flux Value Analysis

            currentDateTime = string(datetime('now'), 'yyyy-MM-dd_HH-mm-SS');
            folderNameRun = sprintf("Task%iSample%iAlpha%0.2fBeta%0.2f_Date%s",iTask,jSample,alpha,beta,currentDateTime);
            mkdir(folderNameRun)
            copyfile('convertToEscherModelWithInputs.py', folderNameRun)
            copyfile("getCode.py", folderNameRun)
            cd(folderNameRun)
            diary logFile.txt

            % FastCC filtering
            tic
            [ATask, ~, ~] = fastcc(TaskSpecificModelCore);
            fastCCTime = toc;

            fprintf("\nFastCC was done in %f seconds.\n", fastCCTime);
            elapsedTime = toc(startTime);
            fprintf("\nTime elapsed since start: %f seconds.\n", elapsedTime);

            ANew = zeros(size(TaskSpecificModelCore.rxns));
            indices = ismember(TaskSpecificModelCore.rxns, TaskSpecificModelCore.rxns(ATask));
            ANew(indices) = 1;
            reactionsNotPartOfTask = ANew == 0;
            NewReversibleModelForKApproximation = removeRxns(TaskSpecificModelCore, TaskSpecificModelCore.rxns(reactionsNotPartOfTask));
            sampleExpValuesAdjustedTaskReversibleKApproximation = sampleExpValuesAdjustedOld(ANew == 1);

            save('NewReversibleModelForKApproximation', 'NewReversibleModelForKApproximation')
            save('sampleExpValuesAdjustedTaskReversibleKApproximation', 'sampleExpValuesAdjustedTaskReversibleKApproximation')
            %load('sampleExpValuesAdjustedTaskReversibleKApproximation')

            [NewIrrevModelKapprox, matchRev, rev2irrevIrrevKapprox, irrev2rev] = convertToIrreversible(NewReversibleModelForKApproximation);
            expValuesMeanAdjustedTaskIrrevKapprox = sampleExpValuesAdjustedTaskReversibleKApproximation;

            for g=1:size(rev2irrevIrrevKapprox,1)
                if size(rev2irrevIrrevKapprox{g},2) == 2
                    expValuesMeanAdjustedTaskIrrevKapprox(rev2irrevIrrevKapprox{g}(2)) = sampleExpValuesAdjustedTaskReversibleKApproximation(rev2irrevIrrevKapprox{g}(1));
                end
            end

            currentDateTime = string(datetime('now'), 'yyyy-MM-dd_HH-mm-SS');
            %filenamenewIrrevModelKapprox = sprintf('LastCreatednewIrrevModelKapprox_Task%iSample%i_date_%s.mat',iTask,jSample,currentDateTime);
            save('NewIrrevModelKapprox', 'NewIrrevModelKapprox')

            %filenamerev2irrevIrrevKapprox = sprintf('LastCreatedrev2irrevIrrevKapprox_Task%iSample%i_date_%s.mat',iTask,jSample,currentDateTime);
            save('rev2irrevIrrevKapprox', 'rev2irrevIrrevKapprox')

            %filenameexpValuesMeanAdjustedTaskIrrevKapprox = sprintf('LastCreatedexpValuesMeanAdjustedTaskIrrevKapprox_Task%iSample%i_date_%s.mat',iTask,jSample,currentDateTime);
            save('expValuesMeanAdjustedTaskIrrevKapprox', 'expValuesMeanAdjustedTaskIrrevKapprox')


            [solutionIrrevFilterLoopLawKapproxTROR_Adjustment,scoreIrrevFilterLoopLawKapproxTROR_Adjustment,reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment, LogData] = runKApproxmitationIrrevTROR_AdjustmentV4(NewIrrevModelKapprox,expValuesMeanAdjustedTaskIrrevKapprox,epsilon,rev2irrevIrrevKapprox,alpha,beta,iTask,jSample,osenseStr,irrev,orphanHandling,lbFromFVA);
            ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndIrrevTROR_Adjustment(NewIrrevModelKapprox, solutionIrrevFilterLoopLawKapproxTROR_Adjustment,expValuesMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapproxTROR_Adjustment , iTask, osenseStr, irrev,jSample, alpha, beta, epsilon);
            writecell(ActiveReactionFluxesFilterIrrevKApprox,"ActiveReactionsTask_ConvergedSolution.xlsx",'Sheet','ActiveReactionsTask','Range','A1');


            % Function to use optmizeCB to remove any redundant reactions from model
            %prunedModel = removeRedundantReactionsFromSolution(reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment,expValueMeanAdjustedTaskIrrevKapprox, NewIrrevModelKapprox);

            turnOnLoopConstraints = true;
            [newSolution, prunedModel, originalModel, adjustedExpressionData, logDataIMAT] = runRxnsMinimalizationUsingIMAT(NewIrrevModelKapprox, solutionIrrevFilterLoopLawKapproxTROR_Adjustment, expValuesMeanAdjustedTaskIrrevKapprox, length(NewIrrevModelKapprox.rxns), epsilon, turnOnLoopConstraints);
            activeReactionFluxesFilterIrrevKApproxPruned = getActiveReactionsWithFluxesFromModelAndIrrevIMATMinimization(prunedModel, newSolution, adjustedExpressionData,scoreIrrevFilterLoopLawKapproxTROR_Adjustment, iTask, osenseStr, irrev, jSample, alpha, beta, epsilon);
            writecell(ActiveReactionFluxesFilterIrrevKApprox,"ActiveReactionsTask_IMAT_Min_ConvergedSol.xlsx",'Sheet','ActiveReactionsIMATMinimized','Range','A1');

            totalTime = toc;
            logDataOut = sprintf("\nTotal Time required is %f \n", totalTime);
            myLogData = strcat(LogData, logDataIMAT, logDataOut);
            writematrix(myLogData, "myLog.txt");

            diary off

            command = sprintf('python3 getCode.py');
            system(command);
            command = sprintf('python3 convertToEscherModelWithInputs.py "%s"', "");  % Replace 'python' with 'python3' if needed
            [status, cmdout] = system(command);
            disp(cmdout)

            diary logFile2.txt
        end
    end
end