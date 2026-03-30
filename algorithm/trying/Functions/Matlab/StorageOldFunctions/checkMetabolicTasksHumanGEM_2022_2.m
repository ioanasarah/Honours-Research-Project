function [taskReport, essentialRxns, taskStructure,usedRxns]=checkMetabolicTasksHumanGEM_2022_2(model,inputFile,printOutput,printOnlyFailed,printDetails,getEssential,taskStructure,saveUsed)
% Modified to get usedRxns, fixed some errors in reaction naming, and fixed
% to use essentialRxnsTasks2022_2 and generateTaskStructure_2022
% 
% Performs a set of simulations as defined in a task file
% to check if the model is able to pass a list of metabolic tasks.
% A metabolic task is defined as the capacity of producing a list of
% output products when ONLY a defined list of inputs substrates are
% available. In other words, a model successfuly pass a metabolic task if it is
% still a solvable LP problem when ONLY exchange reactions related to
% the inputs substrates and outputs products are allowed to carry fluxes
%
% USAGE:
%
%    [taskReport, essentialRxns, taskStructure,usedRxns]=checkMetabolicTasksHumanGEM_2022_2(model,inputFile,printOutput,printOnlyFailed,printDetails,getEssential,taskStructure,saveUsed)
%
% INPUTS:
%    model:              a model structure
%    inputFile:          a task list in Excel format. See the function
%                        parseTaskList for details (opt if taskStructure is
%                        supplied)
%
% OPTIONAL INPUTS:
%    printOutput:        true if the results of the test should be displayed
%                        (default - true)
%    printOnlyFailed:    true if only tasks that failed should be displayed
%                        (default - false)
%    getEssential:       true if the minimal number of reactions that need to be
%                        active to pass the task need to be computed
%                        (default - false)
%    taskStructure:      structure with the tasks, as from `generateTaskStructure`. If
%                        this is supplied then inputFile is ignored
%    saveUsed:           boolean (false by default), should the function
%                       return usedRxns
%
% OUTPUTS:
%    taskReport:         structure with the results:
%
%                          * firstcolumn - id - cell array with the id of the task
%                          * secondcolumn - description - cell array with the description of the task
%                          * thirdcolumn - ok - boolean array with true if the task was successful
%
%    essentialRxns:      cell array containing the essential reactions required
%                        to pass a task
%
%    taskStructure:      structure with the tasks, as from `generateTaskStructure`.
%   NOTE: ALLMETS can be used to freely allow import/excretion of all
%   metabolites (which have an exchange reaction in the model!). Bounds
%   specified for other metabolites in the task will overwrite this. In
%   addition, non-reversible exchanges will not be made reversible.
% .. Authors:
%       - Originally written for RAVEN toolbox by Rasmus Agren, 2013-11-17
%       - Adapted for cobratoolbox and modified to rely only on flux constraints by Richelle Anne, 2017-05-18
%       - Modified by Sonia Balan then by Bastien Nihant in 2022.
if nargin < 8 || ~exist("saveUsed",'var')
    saveUsed = false;
end
if nargin < 3 || isempty(printOutput)
    printOutput=true;
end
if nargin < 4 || isempty(printOnlyFailed)
    printOnlyFailed=false;
end
if nargin < 5 || isempty(printDetails)
    printDetails=false;
end
if nargin < 6 || isempty(getEssential)
    getEssential=false;
end
% Generate a task structure from a list of task in excell format
if nargin < 7 || isempty(taskStructure)
  taskStructure=generateTaskStructure_2022(inputFile); 
end

% Check the format of the model
if size(model.rxns,2)>size(model.rxns,1)
    model.rxns=model.rxns';
end

%Find all exchange/demand/sink reactions
Exchange = {};

for k=1:length(model.rxns)
    if sum(abs(model.S(:,k))) == 1  
        Exchange(end+1) = model.rxns(k);
    end
end
Exchange=unique(Exchange);
ExchangeInd = findRxnIDs(model,Exchange);
nonExchangesInd = setdiff(1:length(model.lb),ExchangeInd) ;

%Checking for lower bounds above 0 in exchanges.
if sum(model.lb(ExchangeInd) > 0) > 0
    msg = 'At least one exchange reaction has a lower bound higher than 0. This will make tasks with "ALLMETS" fail';
    warning(msg);
end
%Checking for lower bounds above 0 in non-exchanges.
if sum(model.lb(nonExchangesInd) > 0) > 0
    msg = 'At least one non-exchange reaction has a lower bound higher than 0. This could make ANY task fail';
    warning(msg);
end
%Close all exchange reactions
previousLB = model.lb;
previousUB = model.ub;
model.lb(findRxnIDs(model,Exchange))=0;
model.ub(findRxnIDs(model,Exchange))=0;

score=0;
totalTask=0;
notPresent=0;
taskReport = cell(numel(taskStructure),3);

essentialRxns=cell(numel(taskStructure),1);
usedRxns = cell(numel(taskStructure),1);
metabolites={};
%Checking for lower bounds above 0.

for i=1:numel(taskStructure)
    
    tModel=model;
    modelMets=upper(tModel.mets);
    
    %%SETUP of the input model
    %suppress objective function if any
    tModel.c(tModel.c==1)=0;

    if isfield(model,'csense')
        if size(tModel.csense,2)>size(tModel.csense,1)
            tModel.csense=tModel.csense(:);
        end
        tModel.csense(length(model.b),1) = 'E';
    end

    taskReport{i,1}=taskStructure(i).id;
    taskReport{i,2}=taskStructure(i).system;
    taskReport{i,3}=taskStructure(i).subsystem;
    taskReport{i,4}=taskStructure(i).description;
    
    %Set the inputs
    if ~isempty(taskStructure(i).inputs)
       
        rxn_Subs={};
        %Adding code to resolve ALLMETS in inputs prior to any other input
        if sum(strcmp(taskStructure(i).inputs,'ALLMETS')) >=1
            tModel.lb(findRxnIDs(tModel,Exchange)) = max(previousLB(findRxnIDs(model,Exchange)),-1000);% by convention in Human1, -1000 means unrestricted input for all metabolites.
            %tModel.ub(findRxnIDs(tModel,Exchange)) = 0;
            inputIDsToUse = 1:length(taskStructure(i).inputs);
            inputIDsToUse = inputIDsToUse(~strcmp(taskStructure(i).inputs,'ALLMETS'));
        else
            inputIDsToUse = 1:length(taskStructure(i).inputs);
        end
        
        for n=inputIDsToUse
            INPUT=taskStructure(i).inputs(n);

            metabolites(end+1)=INPUT;
            INPUT=INPUT{1};
            if strcmp(INPUT(end),']')
                match_INPUTS = strncmpi(INPUT,modelMets,length(INPUT(1:end-3)));
            else
                match_INPUTS = strncmpi(INPUT,modelMets,length(INPUT(1:end-1)));
            end
            
            match_INPUTS = modelMets(match_INPUTS==1);
                        
            compSymbol=cell(1,length(match_INPUTS));
            
            for k=1:length(match_INPUTS)
                if strcmp(match_INPUTS{k}(end),']')
                    tokens = match_INPUTS{k}(end-1); 
                else
                    tokens = match_INPUTS{k}(end); 
                end
                compSymbol{k} = tokens;
            end
        
            % Definition of the compartment for the exchange reaction
%             if strcmp(INPUT(end),']')
%                 comp_used=INPUT(end-1);
%             else
%                 comp_used=INPUT(end);
%             end
            % Set the exchange reactions for the inputs
            AddExchange=0;

            if ismember(upper(INPUT(end-1)),compSymbol) || ismember(upper(INPUT(end)),compSymbol) 

                Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,INPUT));
                Tsp_rxn = full(tModel.S(:,Tsp_ID));
                Nb_React=sum(abs(Tsp_rxn),1);
                
                % If an exchange reaction already exist
                if ~isempty(find(Nb_React==1))
                    ID_exc=find(Nb_React==1);
                   
                    % If the input is also member of the outputs, let the exchange reversible
                    if ismember(INPUT,taskStructure(i).outputs)==1
                        tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -1000, 'l');
                        tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), 1000, 'u');
                        rxn_Subs(end+1) = tModel.rxns(Tsp_ID(ID_exc));
                    else
                        tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -taskStructure(i).UBin(n), 'l');
                        tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -taskStructure(i).LBin(n), 'u');
                        rxn_Subs(end+1) = tModel.rxns(Tsp_ID(ID_exc));
                    end
                    
                else 
                    AddExchange=1;
                end
            else
            	AddExchange=1;
            end
            

            % Add a temporary exchange reaction that allows the import of
            % the metabolite 
            if AddExchange==1
    
                % If the input is also member of the outputs, let the exchange reversible
                if ismember(INPUT,taskStructure(i).outputs)==1
                    [tModel]=addReaction(tModel,['temporary_exchange_',INPUT],[' <=> ',INPUT],[],[],-1000,1000);
                    taskStructure(i).inputs(n)={[INPUT]};
                else
                    [tModel]=addReaction(tModel,['temporary_exchange_',INPUT],[' => ',INPUT],[],[],taskStructure(i).LBin(n),taskStructure(i).UBin(n));
                    taskStructure(i).inputs(n)={[INPUT]};
                end
                rxn_Subs(end+1) = {['temporary_exchange_',INPUT]};
            end           
        end
    end

    modelMets=upper(tModel.mets);
    [I J]=ismember(upper(taskStructure(i).inputs(inputIDsToUse)),modelMets);
    J=J(I);
    
    %Check that all metabolites exist and are defined only once
    if ~all(I)
        fprintf(['ERROR: Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
    	taskReport{i,5}='Could not find all inputs';
    	notPresent=notPresent+1;    
    end
    if numel(J)~=numel(unique(J))
    	display(['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);  
    end

    %Set the outputs
    if ~isempty(taskStructure(i).outputs)

        rxn_Prod={};

        %Adding code to resolve ALLMETS in outputs prior to any other
        %output
        if sum(strcmp(taskStructure(i).outputs,'ALLMETS')) >=1
            tModel.ub(findRxnIDs(tModel,Exchange)) = min(previousUB(findRxnIDs(model,Exchange)),1000);% by convention in Human1, 1000 means unrestricted output for all metabolites.
            %tModel.lb(findRxnIDs(tModel,Exchange)) = 0;
            outputIDsToUse = 1:length(taskStructure(i).outputs);
            outputIDsToUse = outputIDsToUse(~strcmp(taskStructure(i).outputs,'ALLMETS'));
        else
            outputIDsToUse = 1:length(taskStructure(i).outputs);
        end
        for n=outputIDsToUse
            OUTPUT=taskStructure(i).outputs(n);
            metabolites(end+1)=OUTPUT;
        	OUTPUT=OUTPUT{1};
            %%%%%%%%%%%% This chuck added from checkMetabolicTasks_BVD from
            %%%%%%%%%%%% iCardio
            if strcmp(OUTPUT,'ALLMETS')
                % exchange = findRxnIDs(tModel, Exchange);
                for p = 1:length(Exchange)
                    tModel.ub(findRxnIDs(tModel,Exchange(p)))=1000;
                end
            end
            %%%%%%%%%%%%
        	%skip the setup if output is also input as it has already been
        	%setup
            if ismember(upper(OUTPUT),upper(taskStructure(i).inputs))==1
            	continue
            end
            if strcmp(OUTPUT(end),']')
                match_OUTPUTS = strncmpi(OUTPUT,modelMets,length(OUTPUT(1:end-3)));
            else
                match_OUTPUTS = strncmpi(OUTPUT,modelMets,length(OUTPUT(1:end-1)));
            end
            
            match_OUTPUTS = modelMets(match_OUTPUTS==1);
            compSymbol=cell(1,length(match_OUTPUTS));
            for k=1:length(match_OUTPUTS)
                if strcmp(match_OUTPUTS{k}(end),']')
                    tokens = match_OUTPUTS{k}(end-1);                   
                else
                    tokens = match_OUTPUTS{k}(end); 
                end
                compSymbol{k} = tokens;
            end
        

            % Definition of the compartment for the exchange reaction
%            comp_used=OUTPUT(end);
            % Set the exchange reactions for the outputs
            AddExchange=0;
            if ismember(upper(OUTPUT(end-1)),compSymbol)==1 || ismember(upper(OUTPUT(end)),compSymbol)==1
            	Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,OUTPUT));
                Tsp_rxn = full(tModel.S(:,Tsp_ID));
                Nb_React=sum(abs(Tsp_rxn),1);
                
                % If an exchange reaction already exist
                if ~isempty(find(Nb_React==1))
                	ID_exc=find(Nb_React==1);
                     tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), taskStructure(i).LBout(n), 'l');
                     tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), taskStructure(i).UBout(n), 'u');
                     rxn_Prod(end+1)=tModel.rxns(Tsp_ID(ID_exc));
                else
                    AddExchange=1;
                end
            else

                AddExchange=1;
            end
            
            % Add a temporary exchange reaction that allows the export of
            % the metabolite 
            if AddExchange==1
            	[tModel]=addReaction(tModel,['temporary_exchange_',OUTPUT],[OUTPUT,' => '],[],[],taskStructure(i).LBout(n),taskStructure(i).UBout(n));
                taskStructure(i).outputs(n)={[OUTPUT]};
                rxn_Prod(end+1) = {['temporary_exchange_',OUTPUT]};
            end
        end
    end
    if isfield(taskStructure,'EQU')
        if ~isempty(taskStructure(i).EQU)
        currentTask = taskStructure(i).EQU;
        for p = 1:length(taskStructure(i).EQU)
            equation = currentTask{p,1};
            UB = taskStructure(i).EQUUB(p,1);
            LB = taskStructure(i).EQULB(p,1);
%             if UB == 0 
%                 UB = 1000;
%             end
%             if LB == 0
%                 LB = -1000;
%             end
            % add reaction to the model
            [tModel, rxnIDexists] = addReaction(tModel, strcat('TEMPORARY_',taskStructure(i).id,'_',num2str(p)), 'reactionFormula', equation, 'lowerBound', LB, 'upperBound', UB);
            
            % Check to see if the reaction is already in the model
            % temp = 
        end
        end
    end
    

    modelMets=upper(tModel.mets);
    [I J]=ismember(upper(taskStructure(i).outputs(outputIDsToUse)),modelMets);
    J=J(I);

    %Check that all metabolites exist and are defined only once
    if ~all(I)
        fprintf(['ERROR: Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
        taskReport{i,5}='Could not find all outputs';
        notPresent=notPresent+1;
    end
    if numel(J)~=numel(unique(J))
        display(['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);  
    end

	%Solve the constrained problem

    if isfield(model,'csense')
        if size(tModel.csense,2)>size(tModel.csense,1)
            tModel.csense=tModel.csense(:);
        end
        tModel.csense(length(tModel.mets),1) = 'E';
    end
	tModel.osense = -1;
	tModel.A=tModel.S;
	sol=solveCobraLP(tModel);


    if printDetails==true
        if sol.stat~=0
               SUBS=rxn_Subs;
               PROD=rxn_Prod;
            display 'Reactions associated with substrates'
            printRxnFormula(tModel,rxn_Subs);
            display 'Bounds of reactions associated with substrates'
            [tModel.lb(findRxnIDs(tModel,rxn_Subs)) tModel.ub(findRxnIDs(tModel,rxn_Subs))]
            display 'Flux values of reactions associated with substrates'
            sol.full(findRxnIDs(tModel,rxn_Subs))
            display 'Reactions associated with products'
            printRxnFormula(tModel,rxn_Prod);
            display 'Bounds of reactions associated with products'
            [tModel.lb(findRxnIDs(tModel,PROD)) tModel.ub(findRxnIDs(tModel,rxn_Prod))]
            display 'Flux values of reactions associated with products'
            sol.full(findRxnIDs(tModel,rxn_Prod))

        end 
    end
    
    if ~isempty(sol.full)
        if sum(abs(sol.full))~=0
            if taskStructure(i).shouldFail==0
                taskReport{i,5}='true';
                if printOnlyFailed==false && printOutput==true
                    fprintf(['PASS: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
                    score=score+1;
                end
                %Calculate the minimal number of reactions to pass the task
                if getEssential==true 
                    if saveUsed
                        [Rxns_taskEssential,Rxns_Used]=essentialRxnsTasks2022_2(tModel,true);
                        usedRxns{i} = Rxns_Used;
                        essentialRxns{i}= Rxns_taskEssential';
                    else
                        [Rxns_taskEssential]=essentialRxnsTasks2022_2(tModel);
                        essentialRxns{i}= Rxns_taskEssential';
                    end
                    
                end  
            else
                taskReport{i,5}='PASS (should fail)';
                if printOutput==true
                    fprintf(['PASS (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
                end
            end
        else
            if taskStructure(i).shouldFail==0
                taskReport{i,5}='FAIL (should NOT fail)';
                if printOutput==true
                    fprintf(['FAIL: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
                end
            else
                taskReport{i,5}='FAIL (should fail)';
                if printOnlyFailed==false && printOutput==true
                    fprintf(['FAIL (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
                    score=score+1;
                end
            end
        end            
    else
        if taskStructure(i).shouldFail==0
            taskReport{i,5}='FAIL (should NOT fail)';
            if printOutput==true
                fprintf(['FAIL: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        else
            taskReport{i,5}='FAIL (should fail)';
            if printOnlyFailed==false && printOutput==true
                fprintf(['FAIL (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
                score=score+1;
            end
        end
    end
totalTask=totalTask+1;
end

fprintf(['Pass ',num2str(score),' over ',num2str(totalTask),' metabolic tasks tested','\n']);
fprintf([num2str(notPresent),' failed metabolic task are due to absence of metabolites in the model','\n']);
taskReport{end+1,1}='Final Score';
taskReport{end,5}=[num2str(score),'/',num2str(i)];
taskReport{end+1,1}='Not present in the model';
taskReport{end,5}=num2str(notPresent);
