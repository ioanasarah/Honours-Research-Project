function[score, score_binary ,taskInfos, detailScoring]=CellFie_2022(data,tasksData,model,param)
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
significance=1./expression.count;
significance(isinf(significance))=0;
ScorebyTask=[];
ScorebyTask_binary=[];
display('Compute the task activity score');
for i=1:size(taskInfos,1)
	if ~isempty(essentialRxns{i})
    	rxns=essentialRxns{i};
        rxnID=findRxnIDs(model,rxns);
    	rxnID(rxnID==0)=[];
      	if ~isempty(rxnID)
        	expValue=expressionRxns(rxnID,:);
           	signValue=significance(rxnID,:);
            % if no gene is associated with one of the reaction -
          	% remove the reactions from the count
            noGene = find(sum(expValue,2)==-SampleNumber);
            if ~isempty(noGene)
                signValue(noGene,:)=[];
                expValue(noGene,:)=[];
            end
            noGene = find(isnan(sum(expValue,2)));
            if ~isempty(noGene)
                signValue(noGene,:)=[];
                expValue(noGene,:)=[];
            end
           	if ~isempty(expValue)
            	if size(expValue,1)>1
                    ScorebyTask(i,:)=sum(expValue.*signValue)./size(expValue,1);
                    Val=sum(expValue)./size(expValue,1);
                    ID_up=find(Val>=5*log(2));
                    ScorebyTask_binary(i,:)=zeros(1,SampleNumber);
                    ScorebyTask_binary(i,ID_up)=1;
                else
                    ScorebyTask(i,:)=expValue.*signValue;
                    ID_up=find(expValue>=5*log(2));
                    ScorebyTask_binary(i,:)=zeros(1,SampleNumber);
                    ScorebyTask_binary(i,ID_up)=1;
                end
            else
            	ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
                ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
           	end
        else
        	ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
            ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
    	end
    else
        ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
        ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
   end
end

detailScoring={};
if param.getDetails
    disp('Format the score')
    for j=1:SampleNumber
    
        incR=1;
        for i=1:size(taskInfos,1)
    
            if ~isempty(essentialRxns{i})
                rxns=essentialRxns{i};
                rxnID=findRxnIDs(model,rxns);
                rxnID(rxnID==0)=[];
                if ~isempty(rxnID)
                    for k=1:length(rxnID)
                        %1st column = sample ID
                        detailScoring{incR,((j-1)*8)+1}=j;
                        %2nd column = task ID
                        detailScoring{incR,((j-1)*8)+2}=i;
                        %3th column = task score for this sample
                        detailScoring{incR,((j-1)*8)+3}=ScorebyTask(i,j);
                        %4th column = task score in binary version for this sample
                        detailScoring{incR,((j-1)*8)+4}=ScorebyTask_binary(i,j);
                        %5th column = essential reaction associated to this
                        %task
                        detailScoring{incR,((j-1)*8)+5}=rxns(k);
                        %6th column = expression score associated  to the
                        %reaction listed in the 5th column
                        detailScoring{incR,((j-1)*8)+6}=expression.Rxns(rxnID(k),SampleNumber);
                        %7th column = gene used to determine the expression of the
                        %reaction listed in the 5th column
                        geneName=expression.gene_used(rxnID(k),SampleNumber);
                        detailScoring{incR,((j-1)*8)+7}=geneName;
                        %8th column = original expression value of the gene
                        %listed in the 7th column
                        detailScoring{incR,((j-1)*8)+8}=data.value(strcmp(data.gene,geneName),SampleNumber);
                        incR=incR+1;
                    end
               end
            end
        end

    end
else
    disp('Detailed scores are not computed')
end
score=ScorebyTask;
score_binary=ScorebyTask_binary;
end