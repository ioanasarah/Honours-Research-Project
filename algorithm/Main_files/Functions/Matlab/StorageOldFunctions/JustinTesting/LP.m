
load("/home/jcornelis/Documents/MaCSBio-10.10/General/Functions/Metabolic_tasks/Jelle/Main Functions/expressionRxnsCellfie")
load("/home/jcornelis/Documents/MaCSBio-10.10/General/Functions/Metabolic_tasks/Jelle/Main Functions/significanceCellfie")
load("/home/jcornelis/Documents/MaCSBio-10.10/General/Functions/Metabolic_tasks/Jelle/Main Functions/tasksDataCellfie")
load("/home/jcornelis/Documents/MaCSBio-10.10/General/Functions/Metabolic_tasks/Jelle/Main Functions/modelCellfie")

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


for i=108:108%size(taskStructure,1)%, 1)%param.doParallel)

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


           

            %%% Calculate the weighted minimum reactions (see word
            %%% document, minimum weighted flux, based on expression data
            %%% using 'median' for orphan reactions
            disp( "Weighted minimum LP")
            %taskModelAdjusted = changeRxnBounds(taskModel, { 'MAR08971'}, 0);
            %orphanHandling = "zero";
            orphanHandling = "median";        
                         
            % rxnsToEdit = {'MAR03925_f','MAR04523_b', 'MAR03925_b','MAR04523_f', 'MAR00479','MAR00483'};
            % a =  findRxnIDs(taskModel,rxnsToEdit);
            % expValueMeanAdjustedTask(a) = 0.001;



            [ActiveReactionFluxesLP, reactionsToKeepLP,scoreLP] = runLPProblemOneOverGeneExpression(taskModel, expValueMeanAdjustedTask, i,orphanHandling,j,irrev,osenseStr);
            reactionsToKeepLP = taskModel.rxns(reactionsToKeepLP);
        end
    end
end
