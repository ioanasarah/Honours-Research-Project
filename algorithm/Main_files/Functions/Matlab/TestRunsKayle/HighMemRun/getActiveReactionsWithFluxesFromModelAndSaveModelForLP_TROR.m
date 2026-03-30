function [ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelForLP_TROR(taskModel, solution,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j)


if irrev 
    irrevStr = "Irrev";
else
    irrevStr = "NoIrrev";
end


if osenseStr == "min"
    osenseStr = "Min";
elseif osenseStr == "max"
    osenseStr = "Max";
end




reactionLength = length(taskModel.rxns);


% yiFlagsPlus =  solution.int(1:reactionLength);
% yiFlagsMinus = solution.int(reactionLength+1:2*reactionLength);
% yiFlagsCombined = zeros(reactionLength,1);
% 
% for j = 1:reactionLength
%     if yiFlagsPlus(j) == 1 && yiFlagsMinus(j) == 1
%         fprintf("rxn %s has both forward and backward reactions active\n", j)
%     end
% 
%     if yiFlagsPlus(j) == 1
%         yiFlagsCombined(j) = 1;
%     elseif yiFlagsMinus(j) == 1
%         yiFlagsCombined(j) = -1;
%     end
% end

newVector = zeros(size(solution.v)); 
newVector(abs(solution.v) > 0.001 )= 1;
selectedRows = newVector == 1;
printFlag = false;
tmodel = taskModel;
%tmodel.mets = tmodel.metNames;
b = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );


pattern = 'MAM\d+\[[a-z]\]';
pattern2 = '\[(.*?)\]';

rxnFormulas = cell(length(b),1);
for bIterator = 1:length(b)
    c = string(b(bIterator));
    fullFormula = "";
    splitString = split(c, '->');
    metabolitesBefore = regexp(splitString{1}, pattern, 'match');
    beforeMetsArray = cell(length(metabolitesBefore),1);
    for cIterator = 1: length(metabolitesBefore)
        idx = findMetIDs(taskModel,metabolitesBefore(cIterator));
        metName = taskModel.metNames(idx);
        extension = regexp(metabolitesBefore{cIterator}, pattern2, 'tokens');
        newName = sprintf("%s[%s]",string(metName), string(extension));
        beforeMetsArray(cIterator) = cellstr(newName);
    end

    metabolitesAfter = regexp(splitString{2}, pattern, 'match');
    afterMetsArray = cell(length(metabolitesAfter),1);
    for cIterator = 1: length(metabolitesAfter)
        idx = findMetIDs(taskModel,metabolitesAfter(cIterator));
        metName = taskModel.metNames(idx);
        extension = regexp(metabolitesAfter{cIterator}, pattern2, 'tokens');
        newName = sprintf("%s[%s]",string(metName), string(extension));
        afterMetsArray(cIterator) = cellstr(newName);
    end

    separator1 = ' + ';
    separator2 = ' -> ';
    beforePart = strjoin(beforeMetsArray, separator1);
    afterPart = strjoin(afterMetsArray, separator1);
    fullFormula = [beforePart, separator2, afterPart];
    fullFormula = cellstr(fullFormula);
    rxnFormulas(bIterator) =fullFormula; 
end


%disp(rxnFormulas)





newCellArray = [num2cell(solution.v(selectedRows)), num2cell(expValueMeanAdjustedTask(selectedRows)), num2cell(taskModel.c(selectedRows)), taskModel.rxns(selectedRows), taskModel.rxnNames(selectedRows), rxnFormulas];

validValues = expValueMeanAdjustedTask(selectedRows) ~= -1;
selectedExpression = expValueMeanAdjustedTask(selectedRows);
selectedFluxes = solution.v(selectedRows);

score = sum(selectedExpression(validValues)) / length(selectedFluxes(validValues));
dispName =  sprintf("LPProblemTask%iSample%iScored%.3f%s%s%s",i,j, score, osenseStr,irrevStr);

removeRows = newVector == 0;

newModel = removeRxns(taskModel, taskModel.rxns(removeRows));

solutionName = string("lastValidSolution"+dispName+".mat");
newModelName = string("newModel"+dispName+".mat");
fluxSolutionName = string("fluxTable"+dispName+".mat");
% 
% save(solutionName, "solution");
% save(newModelName, "newModel");
% save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;
scoreLP = score;

end