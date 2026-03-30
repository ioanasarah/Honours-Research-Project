function model = closeExchangeReactions(model);

Oldmodel = model;

if size(Oldmodel.rxns,2)>size(Oldmodel.rxns,1)
    Oldmodel.rxns=Oldmodel.rxns';
end

%Find all exchange/demand/sink reactions
Exchange = {};

for k=1:length(Oldmodel.rxns)

    if strlength(string(Oldmodel.rxns(k)))>2
        if  extractBetween(string(Oldmodel.rxns(k)),1,3) == 'EX_'
            Exchange(end+1) = Oldmodel.rxns(k);
            
        end
    end
    if sum(abs(Oldmodel.S(:,k))) == 1
        Exchange(end+1) = Oldmodel.rxns(k);
        %disp(k)
    end
end
Exchange=unique(Exchange);

ExchangeInd = findRxnIDs(Oldmodel,Exchange);
nonExchangesInd = setdiff(1:length(Oldmodel.lb),ExchangeInd) ;
%Oldmodel.lb(16) = 0;

%Checking for lower bounds above 0 in exchanges.
if sum(Oldmodel.lb(ExchangeInd) > 0) > 0
    msg = 'At least one exchange reaction has a lower bound higher than 0. This will make tasks with "ALLMETS" fail';
    warning(msg);
end

%Checking for lower bounds above 0 in non-exchanges.
if sum(Oldmodel.lb(nonExchangesInd) > 0) > 0
    msg = 'At least one non-exchange reaction has a lower bound higher than 0. This could make ANY task fail';
    warning(msg);
end


%Close all exchange reactions
previousLB = Oldmodel.lb;
previousUB = Oldmodel.ub;
Oldmodel.lb(findRxnIDs(Oldmodel,Exchange))=0;
Oldmodel.ub(findRxnIDs(Oldmodel,Exchange))=0;



model = Oldmodel;
end
