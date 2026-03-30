function [PassTask, model] = openExchangeTaskReactionsAndCheckTask(model,taskStructure,i);
model = closeExchangeReactions(model);


%Set the inputs
if ~isempty(taskStructure(i).inputs)

    rxn_Subs={};
    inputIDsToUse = 1:length(taskStructure(i).inputs);
    inputReactions = taskStructure(i).inputs(inputIDsToUse);
    %inputReactions = inputReactions';
    charList = char(cell2mat(inputReactions));
    inputReactions = findRxnIDs(model,findRxnsFromMets(model,charList));
    
    findRxnsFromMets(taskStructure(i).inputs(inputIDsToUse))
    
    for n=inputIDsToUse
        INPUT=taskStructure(i).inputs(n);

        %metabolites(end+1)=INPUT;
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
        %metabolites(end+1)=OUTPUT;
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






















%
% model = closeExchangeReactions(model);
%     tModel=model;
%     modelMets=upper(tModel.mets);
%
%     %%SETUP of the input model
%     %suppress objective function if any
%     tModel.c(tModel.c==1)=0;
%
%     if isfield(model,'csense')
%         if size(tModel.csense,2)>size(tModel.csense,1)
%             tModel.csense=tModel.csense(:);
%         end
%         tModel.csense(length(model.b),1) = 'E';
%     end
%
%     taskReport{i,1}=taskStructure(i).id;
%     taskReport{i,2}=taskStructure(i).system;
%     taskReport{i,3}=taskStructure(i).subsystem;
%     taskReport{i,4}=taskStructure(i).description;
%
%     %Set the inputs
%     if ~isempty(taskStructure(i).inputs)
%
%         rxn_Subs={};
%         %Adding code to resolve ALLMETS in inputs prior to any other input
%         if sum(strcmp(taskStructure(i).inputs,'ALLMETS')) >=1
%             tModel.lb(findRxnIDs(tModel,Exchange)) = max(previousLB(findRxnIDs(model,Exchange)),-1000);% by convention in Human1, -1000 means unrestricted input for all metabolites.
%             %tModel.ub(findRxnIDs(tModel,Exchange)) = 0;
%             inputIDsToUse = 1:length(taskStructure(i).inputs);
%             inputIDsToUse = inputIDsToUse(~strcmp(taskStructure(i).inputs,'ALLMETS'));
%         else
%             inputIDsToUse = 1:length(taskStructure(i).inputs);
%         end
%
%         for n=inputIDsToUse
%             INPUT=taskStructure(i).inputs(n);
%
%             metabolites(end+1)=INPUT;
%             INPUT=INPUT{1};
%             if strcmp(INPUT(end),']')
%                 match_INPUTS = strncmpi(INPUT,modelMets,length(INPUT(1:end-3)));
%             else
%                 match_INPUTS = strncmpi(INPUT,modelMets,length(INPUT(1:end-1)));
%             end
%
%             match_INPUTS = modelMets(match_INPUTS==1);
%
%             compSymbol=cell(1,length(match_INPUTS));
%
%             for k=1:length(match_INPUTS)
%                 if strcmp(match_INPUTS{k}(end),']')
%                     tokens = match_INPUTS{k}(end-1);
%                 else
%                     tokens = match_INPUTS{k}(end);
%                 end
%                 compSymbol{k} = tokens;
%             end
%
%             % Definition of the compartment for the exchange reaction
%             %             if strcmp(INPUT(end),']')
%             %                 comp_used=INPUT(end-1);
%             %             else
%             %                 comp_used=INPUT(end);
%             %             end
%             % Set the exchange reactions for the inputs
%             AddExchange=0;
%
%             if ismember(upper(INPUT(end-1)),compSymbol) || ismember(upper(INPUT(end)),compSymbol)
%
%                 Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,INPUT));
%                 Tsp_rxn = full(tModel.S(:,Tsp_ID));
%                 Nb_React=sum(abs(Tsp_rxn),1);
%
%                 % If an exchange reaction already exist
%                 if ~isempty(find(Nb_React==1))
%                     ID_exc=find(Nb_React==1);
%
%                     % If the input is also member of the outputs, let the exchange reversible
%                     if ismember(INPUT,taskStructure(i).outputs)==1
%                         tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -1000, 'l');
%                         tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), 1000, 'u');
%                         rxn_Subs(end+1) = tModel.rxns(Tsp_ID(ID_exc));
%                     else
%                         tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -taskStructure(i).UBin(n), 'l');
%                         tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), -taskStructure(i).LBin(n), 'u');
%                         rxn_Subs(end+1) = tModel.rxns(Tsp_ID(ID_exc));
%                     end
%
%                 else
%                     AddExchange=1;
%                 end
%             else
%                 AddExchange=1;
%             end
%
%
%             % Add a temporary exchange reaction that allows the import of
%             % the metabolite
%             if AddExchange==1
%
%                 % If the input is also member of the outputs, let the exchange reversible
%                 if ismember(INPUT,taskStructure(i).outputs)==1
%                     [tModel]=addReaction(tModel,['temporary_exchange_',INPUT],[' <=> ',INPUT],[],[],-1000,1000);
%                     taskStructure(i).inputs(n)={[INPUT]};
%                 else
%                     [tModel]=addReaction(tModel,['temporary_exchange_',INPUT],[' => ',INPUT],[],[],taskStructure(i).LBin(n),taskStructure(i).UBin(n));
%                     taskStructure(i).inputs(n)={[INPUT]};
%                 end
%                 rxn_Subs(end+1) = {['temporary_exchange_',INPUT]};
%             end
%         end
%     end
%
%     modelMets=upper(tModel.mets);
%     [I J]=ismember(upper(taskStructure(i).inputs(inputIDsToUse)),modelMets);
%     J=J(I);
%
%     %Check that all metabolites exist and are defined only once
%     if ~all(I)
%         fprintf(['ERROR: Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
%         taskReport{i,5}='Could not find all inputs';
%         notPresent=notPresent+1;
%     end
%     if numel(J)~=numel(unique(J))
%         display(['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);
%     end
%
%     %Set the outputs
%     if ~isempty(taskStructure(i).outputs)
%
%         rxn_Prod={};
%
%         %Adding code to resolve ALLMETS in outputs prior to any other
%         %output
%         if sum(strcmp(taskStructure(i).outputs,'ALLMETS')) >=1
%             tModel.ub(findRxnIDs(tModel,Exchange)) = min(previousUB(findRxnIDs(model,Exchange)),1000);% by convention in Human1, 1000 means unrestricted output for all metabolites.
%             %tModel.lb(findRxnIDs(tModel,Exchange)) = 0;
%             outputIDsToUse = 1:length(taskStructure(i).outputs);
%             outputIDsToUse = outputIDsToUse(~strcmp(taskStructure(i).outputs,'ALLMETS'));
%         else
%             outputIDsToUse = 1:length(taskStructure(i).outputs);
%         end
%         for n=outputIDsToUse
%             OUTPUT=taskStructure(i).outputs(n);
%             metabolites(end+1)=OUTPUT;
%             OUTPUT=OUTPUT{1};
%             %%%%%%%%%%%% This chuck added from checkMetabolicTasks_BVD from
%             %%%%%%%%%%%% iCardio
%             if strcmp(OUTPUT,'ALLMETS')
%                 % exchange = findRxnIDs(tModel, Exchange);
%                 for p = 1:length(Exchange)
%                     tModel.ub(findRxnIDs(tModel,Exchange(p)))=1000;
%                 end
%             end
%             %%%%%%%%%%%%
%             %skip the setup if output is also input as it has already been
%             %setup
%             if ismember(upper(OUTPUT),upper(taskStructure(i).inputs))==1
%                 continue
%             end
%             if strcmp(OUTPUT(end),']')
%                 match_OUTPUTS = strncmpi(OUTPUT,modelMets,length(OUTPUT(1:end-3)));
%             else
%                 match_OUTPUTS = strncmpi(OUTPUT,modelMets,length(OUTPUT(1:end-1)));
%             end
%
%             match_OUTPUTS = modelMets(match_OUTPUTS==1);
%             compSymbol=cell(1,length(match_OUTPUTS));
%             for k=1:length(match_OUTPUTS)
%                 if strcmp(match_OUTPUTS{k}(end),']')
%                     tokens = match_OUTPUTS{k}(end-1);
%                 else
%                     tokens = match_OUTPUTS{k}(end);
%                 end
%                 compSymbol{k} = tokens;
%             end
%
%
%             % Definition of the compartment for the exchange reaction
%             %            comp_used=OUTPUT(end);
%             % Set the exchange reactions for the outputs
%             AddExchange=0;
%             if ismember(upper(OUTPUT(end-1)),compSymbol)==1 || ismember(upper(OUTPUT(end)),compSymbol)==1
%                 Tsp_ID=findRxnIDs(tModel,findRxnsFromMets(tModel,OUTPUT));
%                 Tsp_rxn = full(tModel.S(:,Tsp_ID));
%                 Nb_React=sum(abs(Tsp_rxn),1);
%
%                 % If an exchange reaction already exist
%                 if ~isempty(find(Nb_React==1))
%                     ID_exc=find(Nb_React==1);
%                     tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), taskStructure(i).LBout(n), 'l');
%                     tModel = changeRxnBounds(tModel,tModel.rxns(Tsp_ID(ID_exc)), taskStructure(i).UBout(n), 'u');
%                     rxn_Prod(end+1)=tModel.rxns(Tsp_ID(ID_exc));
%                 else
%                     AddExchange=1;
%                 end
%             else
%
%                 AddExchange=1;
%             end
%
%             % Add a temporary exchange reaction that allows the export of
%             % the metabolite
%             if AddExchange==1
%                 [tModel]=addReaction(tModel,['temporary_exchange_',OUTPUT],[OUTPUT,' => '],[],[],taskStructure(i).LBout(n),taskStructure(i).UBout(n));
%                 taskStructure(i).outputs(n)={[OUTPUT]};
%                 rxn_Prod(end+1) = {['temporary_exchange_',OUTPUT]};
%             end
%         end
%     end
%     if isfield(taskStructure,'EQU')
%         if ~isempty(taskStructure(i).EQU)
%             currentTask = taskStructure(i).EQU;
%             for p = 1:length(taskStructure(i).EQU)
%                 equation = currentTask{p,1};
%                 UB = taskStructure(i).EQUUB(p,1);
%                 LB = taskStructure(i).EQULB(p,1);
%                 %             if UB == 0
%                 %                 UB = 1000;
%                 %             end
%                 %             if LB == 0
%                 %                 LB = -1000;
%                 %             end
%                 % add reaction to the model
%                 [tModel, rxnIDexists] = addReaction(tModel, strcat('TEMPORARY_',taskStructure(i).id,'_',num2str(p)), 'reactionFormula', equation, 'lowerBound', LB, 'upperBound', UB);
%
%                 % Check to see if the reaction is already in the model
%                 % temp =
%             end
%         end
%     end
%
%
%     modelMets=upper(tModel.mets);
%     [I J]=ismember(upper(taskStructure(i).outputs(outputIDsToUse)),modelMets);
%     J=J(I);
%
%     %Check that all metabolites exist and are defined only once
%     if ~all(I)
%         fprintf(['ERROR: Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
%         taskReport{i,5}='Could not find all outputs';
%         notPresent=notPresent+1;
%     end
%     if numel(J)~=numel(unique(J))
%         display(['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time']);
%     end
%
%     %Solve the constrained problem
%
%     if isfield(model,'csense')
%         if size(tModel.csense,2)>size(tModel.csense,1)
%             tModel.csense=tModel.csense(:);
%         end
%         tModel.csense(length(tModel.mets),1) = 'E';
%     end
%     tModel.osense = -1;
%     tModel.A=tModel.S;
%     sol=solveCobraLP(tModel);

end