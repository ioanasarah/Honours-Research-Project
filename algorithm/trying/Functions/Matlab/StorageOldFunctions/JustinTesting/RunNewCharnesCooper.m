load("expressionRxnsCellfie")
load("significanceCellfie")
load("tasksDataCellfie")
load("modelCellfie")

%%% Adjustment of expression data to include significance
expValue=expressionRxns;
SampleNumber = size(expressionRxns,2);

% if no gene is associated with one of the reaction -
% remove the reactions from the count
noGene = find(sum(expValue,2)==-SampleNumber);
if ~isempty(noGene)
    %signValue(noGene,:)=[];
    expValue(noGene,:)=[];
end
noGene = find(isnan(sum(expValue,2)));
if ~isempty(noGene)
    %signValue(noGene,:)=[];
    expValue(noGene,:)=-1;
end
expValueMean = mean(expValue, 2);
significanceMean = mean(significance, 2);

%expValueMeanAdjusted = expValueMean .* significanceMean;
%expValueMeanAdjusted(expValueMean == -1) = -1;
%noGene = find(isnan(sum(expValueMeanAdjusted,2)));
%if ~isempty(noGene)
%     %signValue(noGene,:)=[];
%     expValueMeanAdjusted(noGene,:)=-1;
% end




taskStructure = tasksData.taskStructure;

ExpValueAdjusted = expValue .*significance;
ExpValueAdjusted(expValue == -1) = -1;
noGene = find(isnan(sum(ExpValueAdjusted,2)));
if ~isempty(noGene)
    %signValue(noGene,:)=[];
    ExpValueAdjusted(noGene,:)=-1;
end

%%% Variable setting (might want to include this in the function later)
epsilon =1;
osenseStr = "min";
withMilpKapprox = false;
irrev = true;
orphanHandling = "median";
taskThreshold   = 5*log(2);

%%% Array setting for saving values later
iSize = size(taskStructure, 1);
jSize = SampleNumber;
ScorebyTask_binaryMILP = zeros(iSize, jSize+1);
reactionsToKeepPerTaskMILP = cell(iSize, jSize+1);
scoresPerTaskMILP = zeros(iSize, jSize+1);
ScorebyTask_binaryLP = zeros(iSize, jSize+1);
reactionsToKeepPerTaskLP = cell(iSize, jSize+1);
scoresPerTaskLP = zeros(iSize, jSize+1);

ScorebyTask_binaryMILPCooper = zeros(iSize, jSize+1);
reactionsToKeepPerTaskMILPCooper = cell(iSize, jSize+1);
scoresPerTaskMILPCooper = zeros(iSize, jSize+1);
ScorebyTask_binaryMILPCooperIrrev = zeros(iSize, jSize+1);
reactionsToKeepPerTaskMILPCooperIrrev = cell(iSize, jSize+1);
scoresPerTaskMILPCooperIrrev = zeros(iSize, jSize+1);

% For testing (when running the loops below, these are overtaken)
j = jSize;
i = 1;
j = 1;
MainFolder = pwd;
addpath(MainFolder)

%%% Create an empty task model that can be used in the "create core MILP +
%%% looplaw method" 
emptyTaskModel = createEmptyTaskModel(model);


%%% NOT YET WORKING (code should work but computation seems to get stuck on
%%% LoopLaw) 
% coreMILPProblemFractional = createMILPproblemFractionalCore(emptyTaskModel,epsilon);
% save('coreMILPProblemFractional', 'coreMILPProblemFractional')
% load('coreMILPProblemFractional')


%%% Code for creating coreMILP problems that already include the LoopLaw
%%% matrix, this code specifically uses an irrev model to do this and is
%%% thus only usuable for the K-approximation approach
[emptyTaskModelIrrev, matchRev, Emptyrev2irrev, Emptyirrev2rev] = convertToIrreversible(model);
%%MILPproblemCoreIrrev = createKMILPProblemCoreEmptydaptedNoYiMinus(emptyTaskModelIrrev, Emptyrev2irrev );
%%%save(" ","MILPproblemCoreIrrev")
%%load("coreMILPProblemAllOpenIrrev")
%coreMILPProblem = createCoreMILPProblem(emptyTaskModel);
%%save("coreMILPProblemAllOpenFinal","coreMILPProblem")
%load("coreMILPProblemAllOpenFinal") % only for testing
%save("ExpValueAdjusted")


for i=1:1%size(taskStructure,1)%, 1)%param.doParallel)

    for j=1:1%SampleNumber % For testing, adjust to 1 or whatever sample number you want
        expValueMeanAdjusted = ExpValueAdjusted (:,j );
        % expValueMeanAdjusted = mean(ExpValueAdjusted,2);

        %%% Creates a taskSpecific model and checks if task passes (Utilizes
        %%% code from CheckMetabolicTasks) 
        [PassTask, taskModelOld] = openExchangeTaskReactionsAndCheckTask(model,taskStructure,i);
        
        %%% Adjusts expression data and core model to correct size (based
        %%% on how many placeholder reactions were removed)
        sizeDiff =numel(emptyTaskModel.rxns) - numel(model.rxns);
        diffInReactions = length(emptyTaskModel.rxns) -length(taskModelOld.rxns);
        taskModelOldCore = adjustTaskModel(taskModelOld, diffInReactions);
        expValueMeanAdjustedTask = [expValueMeanAdjusted; repmat(-1, sizeDiff, 1)]; 

       
        % LP > convert > fastcc > 

        %%% Make Irrev for subsequent tasks
        sizeDiff = numel(taskModelOld.rxns) - numel(expValueMeanAdjusted);
        expValueMeanAdjustedTaskOld = [expValueMeanAdjusted; repmat(-1, sizeDiff, 1)];
        [taskModel, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(taskModelOld);
        expValueMeanAdjustedTask = expValueMeanAdjustedTaskOld;
        for g=1:size(rev2irrev,1)
            if size(rev2irrev{g},2) == 2
                expValueMeanAdjustedTask(rev2irrev{g}(2)) = expValueMeanAdjustedTaskOld(rev2irrev{g}(1));
                %expValueMeanAdjustedTask(rev2irrev{g,1}(2)) = rev2irrev{g,1}(1);
            end
        end


        if PassTask
            %%
            cd(MainFolder)
            %%% OR adjustment code
            alpha = 0.5; % even addition of the OR adjustment as regular obj funciton at 0.5
            % alpha gives the weights to the OR adjustment of the obj
            % function and to the original objective function by setting
            beta = 0;  % total addition to max score if all oprhan transport reactions would be 0


            currentDateTimeFolder = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            FolderNameRun = sprintf("Task%iSample%iAlpha%0.2fBeta%0.2f_Date%s",i,j,alpha, beta,currentDateTimeFolder);
            mkdir(FolderNameRun)
            copyfile('convertToEscherModelWithInputs.py', FolderNameRun)
            cd(FolderNameRun)
            diary logFile.txt

            %%% Code snippet for fastcc filtering >turn Irrervisible > LoopLaw > K approxitation
            %%% LoopLaw on irrerversibible models seems to take enormous
            %%% amounts of time, also for fastCC filtered ones, but maybe
            %%% the combination of fastCC filtering and SURFsara can do the
            %%% trick. NOW ADDED TROR ADJUSTMENT INTO MILP             
            tic
            [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
            toc        
            Anew = zeros(size(taskModelOldCore.rxns));
            indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
            Anew(indices) =1;
            reactionsNotPartOfTask = Anew==0;   
            newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
            expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTaskOld(Anew ==1);

            save('newReversibleModelForKApproxmiation', 'newReversibleModelForKApproxmiation')
            save('expValueMeanAdjustedTaskReversibleKApproximation', 'expValueMeanAdjustedTaskReversibleKApproximation')

            [newIrrevModelKapprox, matchRev, rev2irrevIrrevKapprox, irrev2rev] = convertToIrreversible(newReversibleModelForKApproxmiation);
            expValueMeanAdjustedTaskIrrevKapprox = expValueMeanAdjustedTaskReversibleKApproximation;

            for g=1:size(rev2irrevIrrevKapprox,1)
                if size(rev2irrevIrrevKapprox{g},2) == 2
                    expValueMeanAdjustedTaskIrrevKapprox(rev2irrevIrrevKapprox{g}(2)) = expValueMeanAdjustedTaskReversibleKApproximation(rev2irrevIrrevKapprox{g}(1));
                end
            end

            currentDateTime = datestr(now, 'yyyy-mm-dd___HH-MM-SS');
            %filenamenewIrrevModelKapprox = sprintf('LastCreatednewIrrevModelKapprox_Task%iSample%i_date_%s.mat', i,j,currentDateTime);
            save('newIrrevModelKapprox', 'newIrrevModelKapprox')

            %filenamerev2irrevIrrevKapprox = sprintf('LastCreatedrev2irrevIrrevKapprox_Task%iSample%i_date_%s.mat', i,j,currentDateTime);
            save('rev2irrevIrrevKapprox', 'rev2irrevIrrevKapprox')    

            %filenameexpValueMeanAdjustedTaskIrrevKapprox = sprintf('LastCreatedexpValueMeanAdjustedTaskIrrevKapprox_Task%iSample%i_date_%s.mat', i,j,currentDateTime);
            save('expValueMeanAdjustedTaskIrrevKapprox', 'expValueMeanAdjustedTaskIrrevKapprox')

            %% Charnes cooper using already filtered reversible model
            OrphanReactionValue = median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0)); % set to 0 if orphan reactions should not be included
            %median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0))
            CharnesCooperMILPproblem = createMILPproblemFractional(newReversibleModelForKApproxmiation, expValueMeanAdjustedTaskReversibleKApproximation, epsilon, true,OrphanReactionValue); 
            save('CharnesCooperMILPproblem', 'CharnesCooperMILPproblem')
            %
            varargin = parseSolverParameters('MILP');
            varargin.Heuristics = 0.05; %default 0.05
            varargin.printLevel  = 3;
            varargin.Presolve = 2;
            varargin.Method = 4;
            varargin.Cuts = 3;
            varargin.NodeFileStart = 1;
            varargin.MIPFocus = 2;
            varargin.logFile = 'CharnesCooperLog.txt';
            varargin.NoRelHeurTime = 600;
            BoundsValue = 70;
            OrphanReactionValue = median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0)); % set to 0 if orphan reactions should not be included
            BoundsAdjust = true;
                        

            %
            orphanReactionsIds = find(expValueMeanAdjustedTaskReversibleKApproximation ==-1);
            c_z = CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns));
            c_z(orphanReactionsIds) = OrphanReactionValue;
            CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns)) = c_z;
            CharnesCooperMILPproblem.c(5*length(newReversibleModelForKApproxmiation.rxns)+1:6*length(newReversibleModelForKApproxmiation.rxns)) = c_z;  

            %
            CharnesCooperMILPproblemEdited = CharnesCooperMILPproblem;
            if BoundsAdjust
                tempUB = CharnesCooperMILPproblem.ub(1:length(newReversibleModelForKApproxmiation.rxns));
                tempLB = CharnesCooperMILPproblem.lb(1:length(newReversibleModelForKApproxmiation.rxns));
                tempUB(tempUB==1000) = BoundsValue;
                tempLB(tempLB== -1000) = -BoundsValue;
    
                CharnesCooperMILPproblemEdited.ub(1:length(newReversibleModelForKApproxmiation.rxns)) = tempUB;
                CharnesCooperMILPproblemEdited.lb(1:length(newReversibleModelForKApproxmiation.rxns)) = tempLB;            
            end
           %

           % Define parameter names and values
           parameterNames = {
               'Heuristics', 'printLevel', 'Presolve', 'Method', 'Cuts', 'NodeFileStart', 'MIPFocus', 'logFile', 'NoRelHeurTime', 'BoundsValue', 'BoundsAdjust','OrphanReactionValue'
               };

           parameterValues = {
               varargin.Heuristics, varargin.printLevel, varargin.Presolve, varargin.Method, varargin.Cuts, varargin.NodeFileStart, varargin.MIPFocus, varargin.logFile, varargin.NoRelHeurTime,BoundsValue,BoundsAdjust,OrphanReactionValue
               };



            % Open a text file for writing (create or overwrite)
            fileID = fopen('parameter_values_CharnesCooper.txt', 'w');

            % Check if the file was successfully opened
            if fileID == -1
                error('Error opening the file for writing.');
            end

            % Loop through parameter names and values and write them to the file
            for i = 1:numel(parameterNames)
                fprintf(fileID, '%s = %s\n', parameterNames{i}, mat2str(parameterValues{i}));
            end

            % Close the file
            fclose(fileID);

            
            tic
            rng(1)
            CharnesCooperSolution = solveCobraMILP_Adjusted(CharnesCooperMILPproblemEdited, varargin);
            save('CharnesCooperSolution', 'CharnesCooperSolution')
            sumAverageExpressionCharnesCooper = sum(expValueMeanAdjustedTaskReversibleKApproximation(abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-0.1 &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
            scoreAverageCharnesCooper = sumAverageExpressionCharnesCooper/length(find( abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-0.1 &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
            toc
            ActiveReactionFluxesCharnesCoopers = getActiveReactionsWithFluxesFromModelAndSaveModelCooper(newReversibleModelForKApproxmiation, CharnesCooperSolution,expValueMeanAdjustedTaskReversibleKApproximation,scoreAverageCharnesCooper , i, osenseStr , irrev,j);
            writecell(ActiveReactionFluxesCharnesCoopers,"ActiveReactionsTask_CharnesCooper.xlsx",'Sheet','ActReactionsTaskCharnes','Range','A1');
            %

            % % doesn't work, for some reason requires to turn on some
            % % expression reactions!!!!!!
            % load('CharnesCooperSolution')
            % load('expValueMeanAdjustedTaskReversibleKApproximation')
            % load('newReversibleModelForKApproxmiation')

            % [newSolutionCharnes, prunedModelCharnes, originalModelCharnes, adjustedExpressionDataCharnes] = runRxnsMinimalizationUsingIMATForCharnesCooper(newReversibleModelForKApproxmiation, CharnesCooperSolution, expValueMeanAdjustedTaskReversibleKApproximation,length(newReversibleModelForKApproxmiation.rxns), epsilon);
            % ActiveReactionFluxesFilterIrrevKApproxPrunedCharnes = getActiveReactionsModelAndIMATMinimizationCharnes(prunedModelCharnes, newSolutionCharnes,adjustedExpressionDataCharnes,scoreAverageCharnesCooper , i, osenseStr, irrev,j, alpha, beta, epsilon);
            % writecell(ActiveReactionFluxesFilterIrrevKApproxPrunedCharnes,"ActiveReactionsTask_IMAT_Charnes_ConvergedSol.xlsx",'Sheet','ActiveReactionsIMATMinCharnes','Range','A1');
            % 
            % % 
            
            command = sprintf('python convertToEscherModelWithInputs.py "%s"', "");  % Replace 'python' with 'python3' if needed
            [status, cmdout] = system(command);

            diary off
            nowCurrentFolder = pwd;
            if string(nowCurrentFolder) ~= string(MainFolder)
                cd(MainFolder)
            end
        end
    end
end


