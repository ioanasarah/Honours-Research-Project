function[reactionsToKeepPerTaskMILP,reactionsToKeepPerTaskLP,ScorebyTask_binaryMILP,ScorebyTask_binaryLP,scoresPerTaskLP, scoresPerTaskMILP, taskInfos, detailScoring]=CellFie_AverageGeneExpression_2023(data,tasksData,model,param)
% Compute the score associated to each metabolic task listed in taskstructure based on transcriptomic data
%
% USAGE:
%    [score, score_binary ,taskInfos, detailScoring]=CellFie_2022(data,tasksData,model,param)
%
% INPUTS:
%	    data
%       .gene                   cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  mRNA expression data structure (genes x samples)associated to each gene metioned in data.gene
%       .parsedGPR          ParsedGPR for the used model, with the gene
%                                   names in same format as model.genes
%   tasksData
%       .taskStructure      Task structure as used in checkMetabolicTasks
%       .essentialRxns      Cell array of essential rxns for each task as
%                                   returned by checkMetabolicTasks
%
%
%   model                         Reference model used to compute the
%                               metabolic task scores (as a matlab stucture)
% OPTIONAL INPUTS:
%   param.ThreshType            Type of thresholding approach used
%                               (i.e.,'global' or 'local') (default - local)
% related to the use of a GLOBAL thresholding approach - the threshold value is the same for all the genes
%   param.percentile_or_value   the threshold can be defined using a value introduced by the user ('value')
%                               or based on a percentile of the distribution of expression value for all the
%                               genes and across all samples of your
%                               dataset ('percentile')
%   param.percentile            percentile from the distribution of
%                               expression values for all the genes and across all samples that will be
%                               used to define the threshold value
%   param.value                 expression value for which a gene is
%                               considered as active or not (e.g., 5)
%
% related to the use of a LOCAL thresholding approach - the threshold value is different for all the genes
%   param.percentile_or_value   the threshold can be defined using a value introduced by the user ('value')
%                               or based on a percentile of the distribution of expression value of a
%                               specific gene across all samples of your
%                               dataset ('percentile'-default)
%   param.LocalThresholdType    option to define the type of local thresholding approach to use
%                               - 'minmaxmean' (default options )- the threshold for a gene is determined by the mean of expression
%                                  values observed for that gene among all the samples, tissues, or conditions BUT
%                                   the threshold :(i) must be higher or equal to a lower bound and (ii) must be lower
%                                   or equal to an upper bound.
%                               - 'mean' -the threshold for a gene is defined as the mean expression value
%                                  of this gene across all the samples, tissues, or conditions
%   param.percentile_low        lower percentile used to define which gene
%                               are always inactive in the case of use 'MinMaxMean' local thresholding
%                               approach (default = 25)
%   param.percentile_high       upper percentile used to define which gene
%                               are always active in the case of use 'MinMaxMean' local thresholding
%                               approach (default= 75)
%   param.value_low             lower expression value used to define which gene
%                               are always inactive in the case of use 'MinMaxMean' local thresholding
%                               approach (e.g., 5)
%   param.value_high            upper expression value used to define which gene
%                               are always active in the case of use 'MinMaxMean' local thresholding
%                               approach (e.g., 5)
%
% related to the gene mapping approach used
%   param.minSum:               instead of using min and max, use min for AND and Sum
%                               for OR (default: false, i.e. use min)
% Related to parallel computation
%   param.doParallel             Controls the use of parpool for the GPR
%                               portion of the code (time consuming with
%                               high number of subjects.
%                               0: No parallel computation (default)
%                               1: Automatically creates a pool based on
%                               defaults or uses an existing pool.
%                               >1: The number of workers which will be
%                               attempted to be used. Limited by the actual
%                               CPU capacity or by the existing pool's
%                               worker numbers if there is one.
% Related to the required outputs
%   param.getDetails            Boolean, do we want the detailed scoring
%                                (extremely memory intensive for large models and large numers of
%                                samples.
% OUTPUTS:
%   score                       relative quantification of the activity of a metabolic task in a specific condition
%                               based on the availability of data for multiple conditions
%   score_binary                binary version of the metabolic task score
%                               to determine whether a task is active or inactive in specific
%                               conditions
%   taskInfos                   Description of the metabolic task assessed
%   detailScoring               Matrix detailing the scoring
%       1st column = sample ID
%       2nd column = task ID
%       3th column = task score for this sample
%       4th column = task score in binary version for this sample
%       5th column = essential reaction associated to this task
%       6th column = expression score associated  to the reaction listed in the 5th column
%       7th column = gene used to determine the expression of the reaction listed in the 5th column
%       8th column = original expression value of the gene listed in the 7th column
%
% .. Authors:
%       - Anne Richelle, January 2019
%	  - Modified by Eva (MACSBIO intern) to work on iHuman model, date unclear
%	  - Further modified by Bastien Nihant for flexibility and adapting to new cobraToolbox behaviors, November 2022

if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end
if ~exist('model','var')
    error('The reference model has not been defined - please choose a reference model')
end
if ~isfield(data,'parsedGPR')
    data.parsedGPR = GPRparser(model);
end
if ~exist('ref','var')
    ref = model.modelID;
end

if ~exist('param','var') || isempty(param)
    param = struct();
end
if ~isfield(param,'ThreshType')
    param.ThreshType='local';
end
if ~isfield(param,'percentile_or_value')
    param.percentile_or_value='percentile';
end
if ~isfield(param,'LocalThresholdType')
    param.LocalThresholdType='minmaxmean';
end
if ~isfield(param,'percentile_low')
    param.percentile_low=25;
end
if ~isfield(param,'percentile_high')
    param.percentile_high=75;
end
if ~isfield(param,'doParallel')
    param.doParallel = 0;
end
if ~isfield(param,'getDetails')
    param.getDetails = true;
end
%load the info about the task structure
taskInfos=struct2cell(tasksData.taskStructure);
taskInfos=taskInfos';
taskInfos(:,5:end)=[];

%% Depending on the model, load the list of reactions associated with the task
% All these files have the following format
% essentialRxnsbyTask_name_of_model and are located in essentialRxns folder
%load(strcat('essentialRxns/essentialRxnsbyTask_',ref));
essentialRxns = tasksData.essentialRxns;
% check that at least part of list of gene provided are in the model loaded
ID_model=[];
gene_notInModel=[];
for i=1:length(data.gene)
    if isempty(find(strcmp(data.gene{i},model.genes)))
        gene_notInModel(end+1)=i;
    else
        tmpid=find(strcmp(data.gene{i},model.genes)==1);
        ID_model(end+1)=tmpid(1);
    end
end

% introduce a warning in the webtool about how many of the genes provided are
% actually mapped to the model
if isempty(gene_notInModel)
    disp('All genes provided in data are included in the reference model')
else
    disp([num2str(length(gene_notInModel)),' genes provided are not included in the reference model:'])
    if length(gene_notInModel)>10
        disp(data.gene(gene_notInModel(1:10)))
        ntemp = length(gene_notInModel)-10;
        disp(['Plus ',num2str(ntemp),' other genes'])
    else
        disp(data.gene(gene_notInModel))
    end

end

% remove the gene and associated value provided by the user in data that are not in the model
data.gene(gene_notInModel)=[];
SampleNumber = size(data.value,2);
if SampleNumber==1
    data.value(gene_notInModel)=[];
else
    data.value(gene_notInModel,:)=[];
end

% get the threshold value and the histogram for the complete dataset and
% print a figure
if SampleNumber>1
    linData = reshape(data.value,numel(data.value),1);
else
    linData=data.value;
end
linData(linData==0)=[];

% definition of the thresholds
if strcmp(param.ThreshType,'global') && strcmp(param.percentile_or_value,'percentile')
    disp('RUN - global: percentile')
    l_global = (prctile(log10(linData),param.percentile));
    data.ths=10^l_global;
elseif strcmp(param.ThreshType,'global') && strcmp(param.percentile_or_value,'value')
    disp('RUN - global: value')
    data.ths=param.value;
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'mean')
    disp('RUN - local: mean')
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'minmaxmean')&& strcmp(param.percentile_or_value,'percentile')
    disp('RUN - local: minmaxmean: percentile')
    l_high = (prctile(log10(linData),param.percentile_high));
    data.ths_high=10^l_high;
    l_low = (prctile(log10(linData),param.percentile_low));
    data.ths_low=10^l_low;
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'minmaxmean')&& strcmp(param.percentile_or_value,'value')
    disp('RUN - local: minmaxmean: value')
    data.ths_high=param.value_high;
    data.ths_low=param.value_low;
else
    error('No analysis triggered')
end

%% Compute the threshold(s) depending on the approach used
Gene_score=[];
switch param.ThreshType
    case 'local'
        if strcmp(param.LocalThresholdType,'mean')
            %the threshold for each gene is equal to its mean value over
            %all the samples
            threshold=mean(data.value,2)';
        else
            threshold=[];
            for i=1:length(data.gene)
                expressionValue=data.value(i,:);
                if  mean(expressionValue)>=data.ths_high
                    threshold(i)=data.ths_high;
                else
                    threshold(i)=max(mean(expressionValue),data.ths_low);
                end
            end
        end
        % every single gene is associated to an expression score
        for i=1:SampleNumber
            Gene_score(:,i)=5.*log(1+(data.value(:,i)./threshold'));
        end
    case 'global'
        Gene_score=5.*log(1+(data.value./data.ths));
end

% Mapping of the expression data to the model
expression.gene=data.gene;
expression.Rxns=nan(length(model.rxns),SampleNumber);
expression.gene_used=cell(length(model.rxns),SampleNumber);
expression.count=nan(length(model.rxns),SampleNumber);
minSum = false;

disp('Load GPR parse');
%% load parsedGPR for each model
tic
%deal with parallelisation

if isnumeric(param.doParallel) && param.doParallel > 1 && size(gcp("nocreate"),1) == 0
    ppool = parpool(param.doParallel);
elseif  isnumeric(param.doParallel) && param.doParallel > 1 && size(gcp("nocreate"),1) ~= 0
    ppool.NumWorkers = param.doParallel;
elseif isnumeric(param.doParallel) && param.doParallel == 0
    ppool.NumWorkers = param.doParallel;
elseif isnumeric(param.doParallel) && param.doParallel == 1
    ppool = gcp();
end

%load(strcat('parsedGPR/parsedGPR_',ref));
%%parsedGPR = GPRparser(model,minSum);%code to compute the parsed GPR using
%%cobratoolbox
parsedGPR = data.parsedGPR;
disp('Mapping of the expression data to the model');
tmp_rxns = expression.Rxns;
tmp_geneUsed  =expression.gene_used;
tmp_count  =expression.count;
env = getEnvironment();
parfor (i=1:SampleNumber,ppool.NumWorkers)
    restoreEnvironment(env)
    %disp(num2str(i))
    expr = expression;
    if SampleNumber==1
        expr.value=Gene_score;
    else
        expr.value=Gene_score(:,i);
    end
    % Find wich genes in expression data are used in the model
    [gene_id, gene_expr] = findUsedGenesLevels(model,expr);
    % Link the gene to the model reactions
    [expressionRxns, gene_used] = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, minSum);


    if ~strcmp('MT_iHuman.mat', ref) && ~strcmp('MT_Recon3D_entrez.mat', ref) && ~strcmp('ihuman', ref) && ~isempty(str2num(gene_used{1}{1}))
        gene_all=[];
        for j=1:length(gene_used)
            if ~isempty(gene_used{j})
                gene_all(end+1)=str2num(gene_used{j}{1});
            end
        end
    else
        gene_all={};
        for j = 1:length(gene_used)
            if ~isempty(gene_used{j})
                gene_all{end+1} = gene_used{j}{1};
            end
        end
    end
    countGene = tabulate(gene_all);
    if ~strcmp('MT_iHuman.mat', ref) && ~strcmp('MT_Recon3D_entrez.mat', ref) && ~strcmp('ihuman', ref) && ~isempty(str2num(gene_used{1}{1}))
        count=[];
        for k=1:length(gene_used)
            if ~isempty(gene_used{k})
                tmp=countGene(str2num(gene_used{k}{1}),2);
                count(k)=tmp;
            else
                count(k)=0;
            end
        end
    else
        count=nan(1,length(model.rxns));
        keys = (countGene(:, 1));
        %         if  strcmp('MT_iHuman.mat', ref)
        val = cell2mat(countGene(:, 2));
        %         else
        %             val = countGene(:, 2);
        %         end
        geneMap = containers.Map(keys,val);
        for k=1:length(gene_used)
            if ~isempty(gene_used{k})
                tmp= geneMap(gene_used{k}{1});
                count(k)=tmp;
            end
        end
    end
    tmp_rxns(:,i) = expressionRxns;
    %expression.Rxns =expressionRxns;
    tmp_geneUsed(:,i) = gene_used';
    % expression.gene_used(:,i)=gene_used' ;
    tmp_count(:,i) = count';
    %expression.count(:,i)=count';
end
expression.Rxns = tmp_rxns;
expression.gene_used = tmp_geneUsed;
expression.count = tmp_count;
%% Compute the score
expressionRxns=expression.Rxns;
assignin('base','expression',expression)
significance=1./expression.count;
significance(isinf(significance))=0;
ScorebyTask=[];
ScorebyTask_binary=[];
toc
display('Compute the task activity score');


assignin('base','taskInfos', taskInfos)
%
% save("expressionRxnsCellfie","expressionRxns" )
% save("significanceCellfie","significance")
% save("tasksDataCellfie","tasksData")
% save("modelCellfie","model")
%
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
epsilon = 0.001;
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
i = 18;
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
changeCobraSolverParams( 'LP' ,'feasTol',1e-8)

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
i =1;
j = 1;

for i=77:size(taskStructure,1)%, 1)%param.doParallel)

    for j=1:1%SampleNumber % For testing, adjust to 1 or whatever sample number you want
        expValueMeanAdjusted = ExpValueAdjusted (:,j );
        % expValueMeanAdjusted = mean(ExpValueAdjusted,2);

        % TODO
        median_value = median(expValueMeanAdjusted(expValueMeanAdjusted>0));


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
            %
            % [ATask1, fastCCModel1, V] = fastcc(taskModel,1);
            % Anew1 = zeros(size(taskModel.rxns));
            % indices1 = ismember(taskModel.rxns, taskModel.rxns(ATask));
            % Anew1(indices1) =1;
            % reactionsNotPartOfTask1 = Anew1 ==0;
            % newReversibleModelForKApproxmiation1 = removeRxns(taskModel, taskModel.rxns(reactionsNotPartOfTask));
            % exBLpValueMeanAdjustedTaskReversibleKApproximation1 = expValueMeanAdjustedTask(Anew1 ==1);

            %
            % [TEMP1, TEMP4, TEMP2, TEMP3] = convertToIrreversible(taskModelOldCore);
            % revVector = zeros(length(taskModelOldCore.rxns),1);
            % for ghgh = 1:length(TEMP2)
            %     if size(TEMP2{ghgh}, 2) == 2
            %         revVector(ghgh) = 1;
            %     end
            % end
            %
            % consrevVectoristent = swiftcc(taskModelOldCore.S,revVector, 'gurobi');
            % coreInd = length(model.rxns)+1:length(model.rxns)+sizeDiff;
            % weights = ones(length(taskModelOldCore.rxns),1);
            % tic
            % [reconstruction, reconInd, LP] = swiftcore(taskModelOldCore, coreInd, weights, epsilon, false, 'gurobi');
            % toc




            %%
            cd(MainFolder)
            %%% OR adjustment code
            alpha = 0.5; % even addition of the OR adjustment as regular obj funciton at 0.5
            % alpha gives the weights to the OR adjustment of the obj
            % function and to the original objective function by setting
            beta = 0;  % total addition to max score if all oprhan transport reactions would be 0
            lbFromFVA = 1000;

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
            [ATask, fastCCModel, V] = fastcc(taskModelOldCore);
            toc
            Anew = zeros(size(taskModelOldCore.rxns));
            indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
            Anew(indices) =1;
            reactionsNotPartOfTask = Anew==0;
            newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
            expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTaskOld(Anew ==1);

            save('newReversibleModelForKApproxmiation', 'newReversibleModelForKApproxmiation')
            save('expValueMeanAdjustedTaskReversibleKApproximation', 'expValueMeanAdjustedTaskReversibleKApproximation')
            %load('expValueMeanAdjustedTaskReversibleKApproximation')
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
            %%
            
            vararginFVA.osenseStr = 'min';
            vararginFVA.optPercentage = 100;
            vararginFVA.allowLoops = 'original';
            tic
            %[minFlux, maxFlux, Vmin, Vmax] = fluxVariability(newIrrevModelKapprox, vararginFVA);
            %lbFromFVA = maxFlux;
            toc

            %%% Testing
            %load("rev2irrevIrrevKapprox")
            %load("newIrrevModelKapprox")
            %load("expValueMeanAdjustedTaskIrrevKapprox")


            %%%

            % %%% FOR MARIAN CHECK ONLY TEMP
            % newIrrevModelKapproxAdjustedMarian =  changeRxnBounds(newIrrevModelKapprox, newIrrevModelKapprox.rxns(402), 30); %MAR06916
            % newIrrevModelKapproxAdjustedMarian =  changeRxnBounds(newIrrevModelKapproxAdjustedMarian, newIrrevModelKapproxAdjustedMarian.rxns(3489), 30); %temporary_exchange_MAM02751[c]
            % [solutionIrrevFilterLoopLawKapproxTROR_Adjustment,scoreIrrevFilterLoopLawKapproxTROR_Adjustment,reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment] = runKApproxmitationIrrevTROR_Adjustment (newIrrevModelKapproxAdjustedMarian, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox, alpha, beta);
            % ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndIrrevTROR_Adjustment(newIrrevModelKapproxAdjustedMarian, solutionIrrevFilterLoopLawKapproxTROR_Adjustment,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapproxTROR_Adjustment , i, osenseStr, irrev,j, alpha, beta, epsilon);

            % %%%% FOR MARIAN CHECK 2 EDIT TASK
            % listTOMake31 = { 'temporary_exchange_MAM02751[c]','temporary_exchange_MAM01285[c]','temporary_exchange_MAM02039[c]','temporary_exchange_MAM01371[c]'};
            % IDsToMake31 = findRxnIDs(newIrrevModelKapprox, listTOMake31 );
            % newIrrevModelKapproxAdjustedMarian2 =  changeRxnBounds(newIrrevModelKapprox, newIrrevModelKapprox.rxns(IDsToMake31), 31.5);
            % ID37 = findRxnIDs(newIrrevModelKapproxAdjustedMarian2 , {'temporary_exchange_MAM02040[c]'});
            % newIrrevModelKapproxAdjustedMarian2 =  changeRxnBounds(newIrrevModelKapproxAdjustedMarian2, newIrrevModelKapproxAdjustedMarian2.rxns(ID37), 37.5);
            % [solutionIrrevFilterLoopLawKapproxTROR_Adjustment2,scoreIrrevFilterLoopLawKapproxTROR_Adjustment2,reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment] = runKApproxmitationIrrevTROR_Adjustment (newIrrevModelKapproxAdjustedMarian2, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox, alpha, beta);
            % ActiveReactionFluxesFilterIrrevKApprox2 = getActiveReactionsWithFluxesFromModelAndIrrevTROR_Adjustment(newIrrevModelKapproxAdjustedMarian2, solutionIrrevFilterLoopLawKapproxTROR_Adjustment2,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapproxTROR_Adjustment2 , i, osenseStr, irrev,j, alpha, beta, epsilon);

            % % taskModel1WITHOUTLooplaaw= load('MILPProblemTask1WITHOUTLoopLaw').MILPproblem ;
            % % taskModel1Loops = load('MILPProblemTask1LoopLaw.mat').MILPproblem;
            median(expValueMeanAdjustedTaskIrrevKapprox(expValueMeanAdjustedTaskIrrevKapprox>0));
            tic




            %[solutionMILP, scoreMILP, ActiveReactionFluxesFilterIrrevKApprox] = LP7_based_MILP_K_iteration(newIrrevModelKapprox,expValueMeanAdjustedTaskIrrevKapprox, rev2irrevIrrevKapprox, i, j, epsilon);
            




            init_like_solution_irrev = init_like_MILP(newIrrevModelKapprox, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox);
            selectedRows = round(abs(init_like_solution_irrev.cont(1:length(newIrrevModelKapprox.rxns))),6)>= (epsilon-(epsilon/100));
            validValues = expValueMeanAdjustedTaskIrrevKapprox(selectedRows) ~= -1;
            selectedExpression = expValueMeanAdjustedTaskIrrevKapprox(selectedRows);
            averageScore = sum(selectedExpression(validValues)) / length(selectedExpression(validValues)) ;
            %filename = sprintf('lastValidSolutionTask%iSample%iMILP___Scored_%.3f_Average%.3f_date_%s.mat',i,j, kfeasibletest,averageScore,currentDateTime);
            lastValidSolution = init_like_solution_irrev;
            %save(filename, 'lastValidSolution')
            getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(newIrrevModelKapprox, ...
                lastValidSolution,expValueMeanAdjustedTaskIrrevKapprox, ...
                5 , i, osenseStr, irrev,j, alpha, beta, epsilon);
            
            
            
            
            % [solutionIrrevFilterLoopLawKapproxTROR_Adjustment,scoreIrrevFilterLoopLawKapproxTROR_Adjustment,reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment, LogData] = runKApproxmitationIrrevTROR_AdjustmentV4 (newIrrevModelKapprox, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox, alpha, beta,i,j,osenseStr,irrev,orphanHandling,lbFromFVA);
            % ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndIrrevTROR_Adjustment(newIrrevModelKapprox, solutionIrrevFilterLoopLawKapproxTROR_Adjustment,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapproxTROR_Adjustment , i, osenseStr, irrev,j, alpha, beta, epsilon);
            % writecell(ActiveReactionFluxesFilterIrrevKApprox,"ActiveReactionsTask_ConvergedSolution.xlsx",'Sheet','ActiveReactionsTask','Range','A1');
            % 
            % 
            % 
            % % Function to use optmizeCB to remove any redundant reactions from model
            % %prunedModel = removeRedundantReactionsFromSolution(reactionsToKeepIrrevFilterLoopLawKapproxTROR_Adjustment,expValueMeanAdjustedTaskIrrevKapprox, newIrrevModelKapprox);
            % 
            % TurnOnLoopConstraints = true;
            % [newSolution, prunedModel, originalModel, adjustedExpressionData, logDataIMAT] = runRxnsMinimalizationUsingIMAT(newIrrevModelKapprox, solutionIrrevFilterLoopLawKapproxTROR_Adjustment, expValueMeanAdjustedTaskIrrevKapprox,length(newIrrevModelKapprox.rxns), epsilon,TurnOnLoopConstraints);
            % ActiveReactionFluxesFilterIrrevKApproxPruned = getActiveReactionsWithFluxesFromModelAndIrrevIMATMinimization(prunedModel, newSolution,adjustedExpressionData,scoreIrrevFilterLoopLawKapproxTROR_Adjustment , i, osenseStr, irrev,j, alpha, beta, epsilon);
            % writecell(ActiveReactionFluxesFilterIrrevKApprox,"ActiveReactionsTask_IMAT_Min_ConvergedSol.xlsx",'Sheet','ActiveReactionsIMATMinimized','Range','A1');

            Totaltime = toc;
            % LogDataOut = sprintf("\nTotal Time required is %f \n", Totaltime);
            % myLogData = strcat(LogData, logDataIMAT, LogDataOut);
            % writematrix(myLogData, "myLog.txt");

            diary off
            command = sprintf('python convertToEscherModelWithInputs.py "%s"', "");  % Replace 'python' with 'python3' if needed
            [status, cmdout] = system(command)
            disp(cmdout)




            % diary logFile2.txt
            % %% Charnes cooper using already filtered reversible model
            % %OrphanReactionValue = median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0)); % set to 0 if orphan reactions should not be included
            % %median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0))
            % OrphanReactionValue = median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0)); % THis one doesn't actually work, just necessary to send to charnes cooper, to be fully save, adjust both
            % CharnesCooperMILPproblem = createMILPproblemFractional(newReversibleModelForKApproxmiation, expValueMeanAdjustedTaskReversibleKApproximation, epsilon, true,OrphanReactionValue);
            % save('CharnesCooperMILPproblem', 'CharnesCooperMILPproblem')
            % %
            % varargin = parseSolverParameters('MILP');
            % varargin.Heuristics = 0.05; %default 0.05
            % varargin.printLevel  = 3;
            % varargin.Presolve = 2;
            % varargin.Method = 4;
            % varargin.Cuts = 3;
            % varargin.NodeFileStart = 1;
            % varargin.MIPFocus = 2;
            % varargin.logFile = 'CharnesCooperLog.txt';
            % varargin.NoRelHeurTime = 600;
            % BoundsValue = 70;
            % BoundsAdjust = true;
            % 
            % 
            % %
            % orphanReactionsIds = find(expValueMeanAdjustedTaskReversibleKApproximation ==-1);
            % c_z = CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns));
            % c_z(orphanReactionsIds) = OrphanReactionValue;
            % CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns)) = c_z;
            % CharnesCooperMILPproblem.c(5*length(newReversibleModelForKApproxmiation.rxns)+1:6*length(newReversibleModelForKApproxmiation.rxns)) = c_z;
            % 
            % %
            % CharnesCooperMILPproblemEdited = CharnesCooperMILPproblem;
            % if BoundsAdjust
            %     tempUB = CharnesCooperMILPproblem.ub(1:length(newReversibleModelForKApproxmiation.rxns));
            %     tempLB = CharnesCooperMILPproblem.lb(1:length(newReversibleModelForKApproxmiation.rxns));
            %     tempUB(tempUB==1000) = BoundsValue;
            %     tempLB(tempLB== -1000) = -BoundsValue;
            % 
            %     CharnesCooperMILPproblemEdited.ub(1:length(newReversibleModelForKApproxmiation.rxns)) = tempUB;
            %     CharnesCooperMILPproblemEdited.lb(1:length(newReversibleModelForKApproxmiation.rxns)) = tempLB;
            % end
            % %
            % 
            % % Define parameter names and values
            % parameterNames = {
            %     'Heuristics', 'printLevel', 'Presolve', 'Method', 'Cuts', 'NodeFileStart', 'MIPFocus', 'logFile', 'NoRelHeurTime', 'BoundsValue', 'BoundsAdjust','OrphanReactionValue'
            %     };
            % 
            % parameterValues = {
            %     varargin.Heuristics, varargin.printLevel, varargin.Presolve, varargin.Method, varargin.Cuts, varargin.NodeFileStart, varargin.MIPFocus, varargin.logFile, varargin.NoRelHeurTime,BoundsValue,BoundsAdjust,OrphanReactionValue
            %     };
            % 
            % 
            % 
            % % Open a text file for writing (create or overwrite)
            % fileID = fopen('parameter_values_CharnesCooper.txt', 'w');
            % 
            % % Check if the file was successfully opened
            % if fileID == -1
            %     error('Error opening the file for writing.');
            % end
            % 
            % % Loop through parameter names and values and write them to the file
            % for i = 1:numel(parameterNames)
            %     fprintf(fileID, '%s = %s\n', parameterNames{i}, mat2str(parameterValues{i}));
            % end
            % 
            % % Close the file
            % fclose(fileID);
            % 
            % 
            % tic
            % rng(1)
            % CharnesCooperSolution = solveCobraMILP_Adjusted(CharnesCooperMILPproblemEdited, varargin);
            % save('CharnesCooperSolution', 'CharnesCooperSolution')
            % sumAverageExpressionCharnesCooper = sum(expValueMeanAdjustedTaskReversibleKApproximation(abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-(espilon/100) &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
            % scoreAverageCharnesCooper = sumAverageExpressionCharnesCooper/length(find( abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-(espilon/100) &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
            % toc
            % ActiveReactionFluxesCharnesCoopers = getActiveReactionsWithFluxesFromModelAndSaveModelCooper(newReversibleModelForKApproxmiation, CharnesCooperSolution,expValueMeanAdjustedTaskReversibleKApproximation,scoreAverageCharnesCooper , i, osenseStr , irrev,j);
            % writecell(ActiveReactionFluxesCharnesCoopers,"ActiveReactionsTask_CharnesCooper.xlsx",'Sheet','ActReactionsTaskCharnes','Range','A1');
            % %
            % 
            % % % doesn't work, for some reason requires to turn on some
            % % % expression reactions!!!!!!
            % % load('CharnesCooperSolution')
            % % load('expValueMeanAdjustedTaskReversibleKApproximation')
            % % load('newReversibleModelForKApproxmiation')
            % 
            % % [newSolutionCharnes, prunedModelCharnes, originalModelCharnes, adjustedExpressionDataCharnes] = runRxnsMinimalizationUsingIMATForCharnesCooper(newReversibleModelForKApproxmiation, CharnesCooperSolution, expValueMeanAdjustedTaskReversibleKApproximation,length(newReversibleModelForKApproxmiation.rxns), epsilon);
            % % ActiveReactionFluxesFilterIrrevKApproxPrunedCharnes = getActiveReactionsModelAndIMATMinimizationCharnes(prunedModelCharnes, newSolutionCharnes,adjustedExpressionDataCharnes,scoreAverageCharnesCooper , i, osenseStr, irrev,j, alpha, beta, epsilon);
            % % writecell(ActiveReactionFluxesFilterIrrevKApproxPrunedCharnes,"ActiveReactionsTask_IMAT_Charnes_ConvergedSol.xlsx",'Sheet','ActiveReactionsIMATMinCharnes','Range','A1');
            % %
            % % %
            % 
            % command = sprintf('python convertToEscherModelWithInputs.py "%s"', "");  % Replace 'python' with 'python3' if needed
            % [status, cmdout] = system(command)
            % disp(cmdout)
            % 
            % diary off
            % nowCurrentFolder = pwd;
            % if string(nowCurrentFolder) ~= string(MainFolder)
            %     cd(MainFolder)
            % end


            %%

            %%




            %%




            %%





            %
            % [solutionIrrevFilterLoopLawKapprox,scoreIrrevFilterLoopLawKapprox,reactionsToKeepIrrevFilterLoopLawKapprox] = runKApproxmitationAdaptedNoYMinusWithLoopLaw (newIrrevModelKapprox, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox);
            % ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndSaveModelMinYKapprox(newIrrevModelKapprox, solutionIrrevFilterLoopLawKapprox,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapprox , i, osenseStr, irrev,j);
            % % % solutionIrrevFilterLoopLawKapprox = load("solutionIrrevFilterLoopLawKapprox").lastValidSolution;
            % scoreIrrevFilterLoopLawKapprox = load("scoreIrrevFilterLoopLawKapprox").lastValidK;
            %
            % save("solutionIrrevFilterLoopLawKapprox", "solutionIrrevFilterLoopLawKapprox")
            % save("scoreIrrevFilterLoopLawKapprox", "scoreIrrevFilterLoopLawKapprox")
            %



            % tic
            % [ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
            % toc
            % Anew = zeros(size(taskModelOldCore.rxns));
            % indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
            % Anew(indices) =1;
            % reactionsNotPartOfTask = Anew==0;
            % newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
            % expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTask(Anew ==1);
            %
            % [newIrrevModelKapprox, matchRev, rev2irrevIrrevKapprox, irrev2rev] = convertToIrreversible(newReversibleModelForKApproxmiation);
            % expValueMeanAdjustedTaskIrrevKapprox = expValueMeanAdjustedTaskReversibleKApproximation;
            % for g=1:size(rev2irrevIrrevKapprox,1)
            %     if size(rev2irrevIrrevKapprox{g},2) == 2
            %         expValueMeanAdjustedTaskIrrevKapprox(rev2irrevIrrevKapprox{g}(2)) = expValueMeanAdjustedTaskReversibleKApproximation(rev2irrevIrrevKapprox{g}(1));
            %     end
            % end
            % [solutionIrrevFilterLoopLawKapprox,scoreIrrevFilterLoopLawKapprox,reactionsToKeepIrrevFilterLoopLawKapprox] = runKApproxmitationAdaptedNoYMinusWithLoopLaw (newIrrevModelKapprox, expValueMeanAdjustedTaskIrrevKapprox, epsilon, rev2irrevIrrevKapprox);
            % % solutionIrrevFilterLoopLawKapprox = load("solutionIrrevFilterLoopLawKapprox").lastValidSolution;
            % % scoreIrrevFilterLoopLawKapprox = load("scoreIrrevFilterLoopLawKapprox").lastValidK;
            % ActiveReactionFluxesFilterIrrevKApprox = getActiveReactionsWithFluxesFromModelAndSaveModelMinYKapprox(newIrrevModelKapprox, solutionIrrevFilterLoopLawKapprox,expValueMeanAdjustedTaskIrrevKapprox,scoreIrrevFilterLoopLawKapprox , i, osenseStr, irrev,j);
            % save("solutionIrrevFilterLoopLawKapprox", "solutionIrrevFilterLoopLawKapprox")
            % save("scoreIrrevFilterLoopLawKapprox", "scoreIrrevFilterLoopLawKapprox")
            %


            %%% Original Charnes Cooper formulation with (true) or without
            %%% (false) loopLaw, solving time = 20+ hours
            % MILPproblemFractional = createMILPproblemFractional(taskModel, expValueMeanAdjustedTask, epsilon, true);
            % solutionFractional = solveCobraMILP(MILPproblemFractional);


            %%% Attempt to perform Charnes Cooper fractional problem with
            %%% fastCC filtered model (still doesn't work)
            % tic
            % [A, fastCCModel, V] = fastcc(taskModelOld,epsilon);
            % toc
            % Anew = zeros(size(taskModelOld.rxns));
            % indices = ismember(taskModelOld.rxns, taskModelOld.rxns(A));
            % Anew(indices) =1;
            % Anew2 = Anew==0;
            % fastCCModel = removeRxns(taskModelOld,taskModelOld.rxns(Anew2) );
            % expValueMeanAdjustedTaskCCOld = expValueMeanAdjustedTaskOld(Anew ==1);
            % % MILPproblemFractional = createMILPproblemFractional(fastCCModel, expValueMeanAdjustedTaskCCOld, epsilon, true);
            % % solutionFractional = solveCobraMILP(MILPproblemFractional);


            %%% FastCC filtering and using the filtered (requires above
            %%% code fastCCmodel and expValueMeanAjustedCCOld) to turn into
            %%% irrev model and solve using K approx
            % [fastCCModelIrrev, matchRev, fastCCrev2irrev, irrev2rev] = convertToIrreversible(fastCCModel);
            % expValueMeanAdjustedTaskCC = expValueMeanAdjustedTaskCCOld;
            % for g=1:size(fastCCrev2irrev,1)
            %     if size(fastCCrev2irrev{g},2) == 2
            %         expValueMeanAdjustedTaskCC(fastCCrev2irrev{g}(2)) = expValueMeanAdjustedTaskCCOld(fastCCrev2irrev{g}(1));
            %         %expValueMeanAdjustedTask(rev2irrev{g,1}(2)) = rev2irrev{g,1}(1);
            %     end
            % end
            % [solutionKApproxCCIrrev,scoreMILPKApproxCCIrrev,reactionsToKeepMILPKApproxCC] = runKApproxmitationAdaptedNoYMinus (fastCCModelIrrev, expValueMeanAdjustedTaskCC, epsilon, fastCCrev2irrev);
            % ActiveReactionFluxesMILPCCIrrev = getActiveReactionsWithFluxesFromModelAndSaveModelAdaptedMinY(fastCCModelIrrev, solutionKApproxCCIrrev,expValueMeanAdjustedTask,scoreMILPKApproxCCIrrev , i, osenseStr, irrev,j);


            %%% EssentialReactionBased Fractional MILP solution: this code
            %%% takes the output of essentialReactoins and filters the
            %%% model based on this, then performs both charnes Cooper and
            %%% a regular K approximation (doesn't work for all tasks,
            %%% which is a bit strange, not super useful as of right now)
            %%% necessary for creation of excel tables as of now
            % disp( "EssentialRxns reactions")
            % [essentialRxns,usedRxns] = essentialRxnsTasks2023_1(taskModelOld,true,1e-06);
            % charValues = usedRxns;  % Array of char values
            % 
            % names = taskModelOld.rxns(startsWith(taskModelOld.rxns, 'temporary_exchange_'));
            % indicesExact = [];  % Initialize indices for exact matches
            % indicesFB = [];  % Initialize indices for "_f" and "_b" matches
            % tempExchangeReactions = findRxnIDs(taskModelOld,names);
            % for b = 1:numel(charValues)
            %     charValue = charValues(b);
            %     indicesExact = [indicesExact; find(ismember(taskModelOld.rxns, charValue))];
            %     indicesFB = [indicesFB; findRxnIDs(taskModelOld, {strjoin([charValue '_f'], "")})];
            %     indicesFB = [indicesFB; findRxnIDs(taskModelOld, {strjoin([charValue '_b'], "")})];
            % end
            % indicesFB = [indicesFB;indicesExact;tempExchangeReactions];
            % indicesFB = unique(indicesFB(indicesFB ~= 0));
            % idsCharnes = ones(size(expValueMeanAdjustedTaskOld,1),1);
            % idsCharnes(indicesFB) = 0;
            % taskModelForCharnesCooper = removeRxns(taskModelOld,taskModelOld.rxns(idsCharnes==1));
            % expValueMeanAdjustedTaskCharnesCooper =expValueMeanAdjustedTaskOld(idsCharnes==0);
            % saveUsedReactionModel (taskModelForCharnesCooper,i,j)
            % 
            % disp("MILP Cooper reversible")
            % %MILPproblemFractional = createMILPproblemFractional(taskModelForCharnesCooper, expValueMeanAdjustedTaskCharnesCooper, epsilon, false);
            % %solutionFractional = solveCobraMILP(MILPproblemFractional);
            % solutionFractional.stat =0;
            % if solutionFractional.stat ==1
            %     ActiveReactionFluxesMILPCooper = getActiveReactionsWithFluxesFromModelAndSaveModelCooper(taskModelForCharnesCooper, solutionFractional,expValueMeanAdjustedTaskCharnesCooper,solutionFractional.obj , i, osenseStr, irrev,j);
            % else
            %     nanString = ' ';
            %     ActiveReactionFluxesMILPCooper = repmat({nanString}, 1, 6);
            % end
            % %convert to irrev model for k approx
            % [taskModelCooperIrrev, matchRevCooper, rev2irrevCooper, irrev2revCooper] = convertToIrreversible(taskModelForCharnesCooper);
            % expValueMeanAdjustedTaskCharnesCooperIrrev = expValueMeanAdjustedTaskCharnesCooper;
            % for g=1:size(rev2irrevCooper,1)
            %     if size(rev2irrevCooper{g},2) == 2
            %         expValueMeanAdjustedTaskCharnesCooperIrrev(rev2irrevCooper{g}(2)) = expValueMeanAdjustedTaskCharnesCooper(rev2irrevCooper{g}(1));
            %         %expValueMeanAdjustedTask(rev2irrev{g,1}(2)) = rev2irrev{g,1}(1);
            %     end
            % end
            % rxnFormulas = printRxnFormula(taskModelCooperIrrev,taskModelCooperIrrev.rxns,false,true,true);
            % ActiveReactionUsedReactions = [num2cell(expValueMeanAdjustedTaskCharnesCooperIrrev),taskModelCooperIrrev.rxns ,rxnFormulas];
            % %[solutionUsedRxnsKApprox,scoreUsedRxnsKApprox,reactionsToKeepUsedRxnsKApprox] = runKApproxmitationAdaptedNoYMinusErrorDebugging (taskModelCooperIrrev, expValueMeanAdjustedTaskCharnesCooperIrrev, epsilon, rev2irrevCooper);
            % % ActiveReactionFluxesUsedRxnsKApprox = getActiveReactionsWithFluxesFromModelAndSaveModelUsedRxns(taskModelCooperIrrev, solutionUsedRxnsKApprox,expValueMeanAdjustedTaskCharnesCooperIrrev,scoreUsedRxnsKApprox , i, osenseStr, irrev,j);
            % disp("MILP Cooper Irrev ")
            % %MILPproblemFractionalIrrev = createMILPproblemFractional(taskModelCooperIrrev, expValueMeanAdjustedTaskCharnesCooperIrrev, epsilon,false);
            % %solutionFractionalIrrev = solveCobraMILP(MILPproblemFractionalIrrev);
            % solutionFractionalIrrev.stat =0;
            % if solutionFractionalIrrev.stat ==1
            %     ActiveReactionFluxesMILPCooperIrrev = getActiveReactionsWithFluxesFromModelAndSaveModelCooperIrrev(taskModelCooperIrrev, solutionFractionalIrrev,expValueMeanAdjustedTaskCharnesCooperIrrev,solutionFractionalIrrev.obj , i, osenseStr, irrev,j);
            % else
            %     nanString = ' ';
            %     ActiveReactionFluxesMILPCooperIrrev = repmat({nanString}, 1, 6);
            % end
            % 
            % 
            % %%% Performs a simple minimum reactions LP to get a set of
            % %%% minimum reactions (different from those of Ella and Bastien
            % %%% who uses pFBA to create this minimum reaction set if I
            % %%% recall correctly)
            % disp( "Minimimumreactions")
            % taskModelAdjusted = changeRxnBounds(taskModel,{'MAR04898_f'}, 6 );
            % taskModelAdjusted = changeRxnBounds(taskModelAdjusted,{'MAR07638'}, 30 );
            % taskModelAdjusted = changeRxnBounds(taskModelAdjusted,{'MAR04922_b'}, 6 );
            % taskModelAdjusted = changeRxnBounds(taskModelAdjusted,{'MAR04888_b'}, 36 );
            % 
            % 
            % 
            % 
            % [ActiveReactionFluxesMinimumReactions, reactionsToKeepUsedReactions,scoreUsedReactions] = runLPMinimumReactions(taskModel, expValueMeanAdjustedTask, i,orphanHandling,j,irrev,osenseStr);
            % 
            % 
            % %%% Calculate the weighted minimum reactions (see word
            % %%% document, minimum weighted flux, based on expression data
            % %%% using 'median' for orphan reactions
            % disp( "Weighted minimum LP")
            % %taskModelAdjusted = changeRxnBounds(taskModel, { 'MAR08971'}, 0);
            % %orphanHandling = "zero";
            % orphanHandling = "median";
            % 
            % % rxnsToEdit = {'MAR03925_f','MAR04523_b', 'MAR03925_b','MAR04523_f', 'MAR00479','MAR00483'};
            % % a =  findRxnIDs(taskModel,rxnsToEdit);
            % % expValueMeanAdjustedTask(a) = 0.001;
            % 
            % 
            % 
            % [ActiveReactionFluxesLP, reactionsToKeepLP,scoreLP] = runLPProblemOneOverGeneExpression(taskModel, expValueMeanAdjustedTask, i,orphanHandling,j,irrev,osenseStr);
            % reactionsToKeepLP = taskModel.rxns(reactionsToKeepLP);
            % % ActiveReactionFluxesLP1 = ActiveReactionFluxesLP;
            % % ActiveReactionFluxesLP2 = ActiveReactionFluxesLP;
            % % ActiveReactionFluxesLP3 = ActiveReactionFluxesLP;
            % % ActiveReactionFluxesLP4 = ActiveReactionFluxesLP;
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % %%% If turned on (see variable setting before the loops) then
            % %%% will run the K approximation of this task (can take some
            % %%% time (20 min) for a task).
            % if withMilpKapprox
            %     disp("MILP K approx")
            %     [solution,scoreMILP,reactionsToKeepMILP] = runKApproxmitationAdaptedNoYMinus (taskModel, expValueMeanAdjustedTask, epsilon, rev2irrev);
            %     ActiveReactionFluxesMILP = getActiveReactionsWithFluxesFromModelAndSaveModelAdaptedMinY(taskModel, solution,expValueMeanAdjustedTask,scoreMILP , i, osenseStr, irrev,j);
            % end
            % 
            % %%% Creates excel tables
            % if withMilpKapprox
            %     createExcelOfReactions(i,j,withMilpKapprox,ActiveReactionFluxesLP,ActiveReactionFluxesMinimumReactions,ActiveReactionUsedReactions,ActiveReactionFluxesMILPCooper,ActiveReactionFluxesMILPCooperIrrev,ActiveReactionFluxesMILP);
            %     %reactionsToKeepMILP = taskModel.rxns(reactionsToKeepMILP);
            % else
            %     createExcelOfReactions(i,j,withMilpKapprox,ActiveReactionFluxesLP,ActiveReactionFluxesMinimumReactions,ActiveReactionUsedReactions,ActiveReactionFluxesMILPCooper,ActiveReactionFluxesMILPCooperIrrev);
            %     reactionsToKeepMILP = 0;
            % end
            % 
            % %%% Save the data gotten from all the different approaches
            % if solutionFractional.stat ==1
            %     ScorebyTask_binaryMILPCooper(i,j) = 1;
            %     reactionsToKeepPerTaskMILPCooper{i}(:,j) = {0};
            %     scoresPerTaskMILPCooper(i,j) = solutionFractional.obj;
            % else
            %     ScorebyTask_binaryMILPCooper(i,j) = 0;
            %     reactionsToKeepPerTaskMILPCooper{i}(:,j) = {0};
            %     scoresPerTaskMILPCooper(i,j) = 0;
            % end
            % 
            % if solutionFractionalIrrev.stat ==1
            %     ScorebyTask_binaryMILPCooperIrrev(i,j) = 1;
            %     reactionsToKeepPerTaskMILPCooperIrrev{i}(:,j) = {0};
            %     scoresPerTaskMILPCooperIrrev(i,j) = solutionFractionalIrrev.obj;
            % else
            %     ScorebyTask_binaryMILPCooperIrrev(i,j) = 0;
            %     reactionsToKeepPerTaskMILPCooperIrrev{i}(:,j) = {0};
            %     scoresPerTaskMILPCooperIrrev(i,j) = 0;
            % end
            % 
            % 
            % if withMilpKapprox
            %     if scoreMILP > taskThreshold
            %         ScorebyTask_binaryMILP(i, j) = 1;
            %         reactionsToKeepPerTaskMILP{i}(:, j) = {reactionsToKeepMILP};
            %         scoresPerTaskMILP(i,j) = scoreMILP;
            %     else
            %         ScorebyTask_binaryMILP(i, j) = 0;
            %         reactionsToKeepPerTaskMILP{i}(:, j) = {0};
            %         scoresPerTaskMILP(i,j) = scoreMILP;
            %     end
            % end
            % if scoreLP > taskThreshold
            %     ScorebyTask_binaryLP(i, j) = 1;
            %     reactionsToKeepPerTaskLP{i}(:, j) = {reactionsToKeepLP};
            %     scoresPerTaskLP(i,j) = scoreLP;
            % else
            %     ScorebyTask_binaryLP(i, j) = 0;
            %     reactionsToKeepPerTaskLP{i}(:, j) = 0;
            %     scoresPerTaskLP(i,j)= scoreLP;
            % end

        else
            ScorebyTask_binaryLP(i,j) = 0;
            reactionsToKeepPerTaskLP{i}(:, j) = 0;
            scoresPerTaskLP(i,j) = 0;
            ScorebyTask_binaryMILP(i,j) = 0;
            reactionsToKeepPerTaskMILP{i}(:,j) = {0};
            scoresPerTaskMILP(i,j) = 0;

            ScorebyTask_binaryMILPCooper(i,j) = 0;
            reactionsToKeepPerTaskMILPCooper{i}(:,j) = {0};
            scoresPerTaskMILPCooper(i,j) = 0;
            ScorebyTask_binaryMILPCooperIrrev(i,j) = 0;
            reactionsToKeepPerTaskMILPCooperIrrev{i}(:,j) = {0};
            scoresPerTaskMILPCooperIrrev(i,j) = 0;
        end

    end

end

% Save specific variables to a .mat file
save('ScoresAndReactions_workspace.mat', 'ScorebyTask_binaryLP', 'reactionsToKeepPerTaskLP', 'scoresPerTaskLP', ...
    'ScorebyTask_binaryMILP', 'reactionsToKeepPerTaskMILP', 'scoresPerTaskMILP', ...
    'ScorebyTask_binaryMILPCooper', 'reactionsToKeepPerTaskMILPCooper', 'scoresPerTaskMILPCooper', ...
    'ScorebyTask_binaryMILPCooperIrrev', 'reactionsToKeepPerTaskMILPCooperIrrev', 'scoresPerTaskMILPCooperIrrev');

load('ScoresAndReactions_workspace.mat')

end