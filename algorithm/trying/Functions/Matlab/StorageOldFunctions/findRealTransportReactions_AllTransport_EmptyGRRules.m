function transportRxnIDs = findRealTransportReactions_AllTransport_EmptyGRRules(taskModel)

prunedModel = removeRxns(taskModel, taskModel.rxns(taskModel.grRules~="")); 
a = findTransRxns(prunedModel);
transportRxnIDs = findRxnIDs(taskModel,a);

end




% % Find namesbefore/after
% 
% rxnFormulas = cell(length(b),1)
% for bIterator = 1:length(b)
%     fullFormula = "";
%     pattern = 'MAM\d+\[[a-z]\]';
%     metabolitesBefore = regexp(splitString{1}, pattern, 'match');
%     pattern2 = '\[(.*?)\]';
%     beforeMetsArray = cell(length(metabolitesBefore),1)
%     for cIterator = 1: length(metabolitesBefore)
%         idx = findMetIDs(taskModel,metabolitesBefore(cIterator));
%         metName = taskModel.metNames(idx);
%         extension = regexp(metabolitesBefore{cIterator}, pattern2, 'tokens');
%         newName = sprintf("%s[%s]",string(metName), string(extension));
%         beforeMetsArray(cIterator) = cellstr(newName);
%     end
% 
%     metabolitesAfter = regexp(splitString{2}, pattern, 'match');
%     afterMetsArray = cell(length(metabolitesAfter),1)
%     for cIterator = 1: length(metabolitesAfter)
%         idx = findMetIDs(taskModel,metabolitesAfter(cIterator));
%         metName = taskModel.metNames(idx);
%         extension = regexp(metabolitesAfter{cIterator}, pattern2, 'tokens');
%         newName = sprintf("%s[%s]",string(metName), string(extension));
%         afterMetsArray(cIterator) = cellstr(newName);
%     end
% 
%     separator1 = ' + ';
%     separator2 = ' -> ';
%     beforePart = strjoin(beforeMetsArray, separator1);
%     afterPart = strjoin(afterMetsArray, separator1);
%     fullFormula = [beforePart, separator2, afterPart];
%     fullFormula = cellstr(fullFormula);
% 
% end

% find IDs before/after
% get Metname, get extension
% Put together