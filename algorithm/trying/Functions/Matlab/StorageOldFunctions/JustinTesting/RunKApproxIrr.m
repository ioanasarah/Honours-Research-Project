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

expValueMeanAdjusted = expValueMean .* significanceMean;
expValueMeanAdjusted(expValueMean == -1) = -1;
noGene = find(isnan(sum(expValueMeanAdjusted,2)));
if ~isempty(noGene)
    %signValue(noGene,:)=[];
    expValueMeanAdjusted(noGene,:)=-1;
end
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
load("coreMILPProblemAllOpenFinal") % only for testing


for i=1:size(taskStructure,1)%, 1)%param.doParallel)

    for j=322:322%SampleNumber % For testing, adjust to 1 or whatever sample number you want
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

        % %%% Code for attempting to use fastcc filtering on coreMILP problem
        % %%% (reversible) this SHOULD work, however at some point (at a
        % %%% specific K) the computation gets stuck (12+ hours).
        % tic
        % [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
        % toc        
        % %%load("Atask322") %Used for testing so that fastCC doesn't have to
        % %%be run multiple times
        % 
        % Anew = zeros(size(taskModelOldCore.rxns));
        % indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
        % Anew(indices) =1;
        % reactionsNotPartOfTask = Anew==0;    
        % taskSpecificMILP = createTaskSpecificMILP(coreMILPProblem, taskModelOldCore,expValueMeanAdjustedTask ,epsilon);
        % [taskSpecificMILPFiltered,expValueMeanAdjustedTaskFiltered] = adjustTaskSpecificMILPForCCFiltering(taskSpecificMILP,reactionsNotPartOfTask,taskModelOldCore,expValueMeanAdjustedTask);
        % % this gets stuck at some point
        % [solution,score,reactionsToKeep] = runKApproxmitation (taskSpecificMILPFiltered,taskModelOldCore, expValueMeanAdjustedTaskFiltered, epsilon);
        % ActiveReactionFluxesKApproxFilteredReversible = getActiveReactionsWithFluxesFromModelAndSaveModelLoopFunction(taskModelOldCore, solution,expValueMeanAdjustedTask,score,taskSpecificMILPFiltered , i, osenseStr , irrev);


        % %%% Code for a fastCC charnes cooper approach (this would
        % %%% eventually become like above (a taskSpecificMILP based on a
        % %%% core MILP with a single LoopLaw run), however was never fully
        % %%% developed.
        % tic
        % [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
        % toc        
        % Anew = zeros(size(taskModelOldCore.rxns));
        % indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
        % Anew(indices) =1;
        % reactionsNotPartOfTask = Anew==0;   
        % newReversibleModelForCharnesCooper = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
        % expValueMeanAdjustedTaskReversibleForCharnesCooper = expValueMeanAdjustedTask(Anew ==1);
        % MILPproblemFractionalReversibleFilteredCharnesCooper = createMILPproblemFractional(newReversibleModelForCharnesCooper, expValueMeanAdjustedTaskReversibleForCharnesCooper, epsilon, true);
        % solutionFractionalReversibleFilteredCharnesCooper = solveCobraMILP(MILPproblemFractionalReversibleFilteredCharnesCooper);
        % save("solutionFractionalReversibleFilteredCharnesCooper", 'solutionFractionalReversibleFilteredCharnesCooper')

        
        % 
        % %%% Code for filtering, followed by subsequent K approximation with
        % %%% looplaw on a reversible model (should be like the Core MILP
        % %%% approach, but since there is some experimental filtering and
        % %%% re-adding of LoopLaw matrix, this can be used for testing)
        % %%% Added in a version that doesn't have any objective function 
        % %%% (no maximization) to make performance better, unsure how to use
        % %%% the following data.
        % tic
        % [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
        % toc        
        % Anew = zeros(size(taskModelOldCore.rxns));
        % indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
        % Anew(indices) =1;
        % reactionsNotPartOfTask = Anew==0;   
        % newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
        % expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTask(Anew ==1);
        % [solutionFilteredKapprox,scoreFilteredKapprox,reactionsToKeepFilteredKapprox] = runKApproxmitationLoopsForTestingCAdjusted(newReversibleModelForKApproxmiation, expValueMeanAdjustedTaskReversibleKApproximation, epsilon);
        % ActiveReactionFluxesFilteredKapprox = getActiveReactionsWithFluxesFromModelAndSaveModelRevesibleK(newReversibleModelForKApproxmiation, solutionFilteredKapprox,expValueMeanAdjustedTaskReversibleKApproximation,scoreFilteredKapprox , i, osenseStr, irrev,j);



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

            %%% Code snippet for turn Irrervisible > fastcc filtering > LoopLaw > K approxitation
            %%% LoopLaw on irrerversibible models seems to take enormous
            %%% amounts of time, also for fastCC filtered ones, but maybe
            %%% the combination of fastCC filtering and SURFsara can do the
            %%% trick:
            tic
            [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
            toc        
            Anew = zeros(size(taskModelOldCore.rxns));
            indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
            Anew(indices) =1;
            reactionsNotPartOfTask = Anew==0;   
            newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
            expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTask(Anew ==1);
            
            [newIrrevModelKapprox, matchRev, rev2irrevIrrevKapprox, irrev2rev] = convertToIrreversible(newReversibleModelForKApproxmiation);
            expValueMeanAdjustedTaskIrrevKapprox = expValueMeanAdjustedTaskReversibleKApproximation;
            for g=1:size(rev2irrevIrrevKapprox,1)
                if size(rev2irrevIrrevKapprox{g},2) == 2
                    expValueMeanAdjustedTaskIrrevKapprox(rev2irrevIrrevKapprox{g}(2)) = expValueMeanAdjustedTaskReversibleKApproximation(rev2irrevIrrevKapprox{g}(1));
                end
            end            
            [solutionIrrevFilterLoopLawKapprox,scoreIrrevFilterLoopLawKapprox,reactionsToKeepIrrevFilterLoopLawKapprox] = runKApproxmitationAdaptedNoYMinusWithLoopLaw (newIrrevModelKapprox, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox);
            ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndSaveModelMinYKapprox(newIrrevModelKapprox, solutionIrrevFilterLoopLawKapprox,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapprox , i, osenseStr, irrev,j);
            save("solutionIrrevFilterLoopLawKapprox", "solutionIrrevFilterLoopLawKapprox")
            save("scoreIrrevFilterLoopLawKapprox", "scoreIrrevFilterLoopLawKapprox")

        end
    end
end


