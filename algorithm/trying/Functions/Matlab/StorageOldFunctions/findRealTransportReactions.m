function transportRxnIDs = findRealTransportReactions(taskModel)

prunedModel = removeRxns(taskModel, taskModel.rxns(taskModel.grRules~="")); 
a = findTransRxns(prunedModel);
b = printRxnFormula(prunedModel, a, false, true, false);

pattern = 'MAM\d+\[[a-z]\]';
TransportReactionBinary = zeros(length(b),1);

for bIterator = 1:length(b)
    c = string(b(bIterator));
    metabolites = regexp(c, pattern, 'match');
    splitString = split(c, '->');
    metabolitesBefore = regexp(splitString{1}, pattern, 'match');
    metabolitesBefore = cellfun(@(x) regexprep(x, '\[[a-z]\]', ''), metabolitesBefore, 'UniformOutput', false);
    metabolitesAfter = regexp(splitString{2}, pattern, 'match');
    metabolitesAfter = cellfun(@(x) regexprep(x, '\[[a-z]\]', ''), metabolitesAfter, 'UniformOutput', false);

    if numel(metabolitesBefore) == numel(metabolitesAfter)
        metabolitesBeforeSet = unique(metabolitesBefore);
        metabolitesAfterSet = unique(metabolitesAfter);

        if isequal(metabolitesBeforeSet, metabolitesAfterSet)
            TransportReactionBinary(bIterator) = 1;
        end
    end
end
transportRxns = a(TransportReactionBinary ==1);
transportRxnIDs = findRxnIDs(taskModel,transportRxns);
end