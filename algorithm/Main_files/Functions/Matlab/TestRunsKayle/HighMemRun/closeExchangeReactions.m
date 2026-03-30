function Model = closeExchangeReactions(Model)

OldModel = Model;

if size(OldModel.rxns,2)>size(OldModel.rxns,1)
    OldModel.rxns=OldModel.rxns';
end

% Find all exchange/demand/sink reactions
Exchange = {};

for k = 1:length(OldModel.rxns)

    if strlength(string(OldModel.rxns(k)))>2
        if  extractBetween(string(OldModel.rxns(k)),1,3) == 'EX_'
            Exchange(end+1) = OldModel.rxns(k);
            
        end
    end
    if sum(abs(OldModel.S(:,k))) == 1
        Exchange(end+1) = OldModel.rxns(k);
        %disp(k)
    end
end
Exchange=unique(Exchange);

ExchangeInd = findRxnIDs(OldModel,Exchange);
nonExchangesInd = setdiff(1:length(OldModel.lb),ExchangeInd) ;
%Oldmodel.lb(16) = 0;

% Checking for lower bounds above 0 in exchanges.
if sum(OldModel.lb(ExchangeInd) > 0) > 0
    msg = 'At least one exchange reaction has a lower bound higher than 0. This will make tasks with "ALLMETS" fail';
    warning(msg);
end

% Checking for lower bounds above 0 in non-exchanges.
if sum(OldModel.lb(nonExchangesInd) > 0) > 0
    msg = 'At least one non-exchange reaction has a lower bound higher than 0. This could make ANY task fail';
    warning(msg);
end


% Close all exchange reactions
previousLB = OldModel.lb;
previousUB = OldModel.ub;
OldModel.lb(findRxnIDs(OldModel,Exchange))=0;
OldModel.ub(findRxnIDs(OldModel,Exchange))=0;



Model = OldModel;
end
