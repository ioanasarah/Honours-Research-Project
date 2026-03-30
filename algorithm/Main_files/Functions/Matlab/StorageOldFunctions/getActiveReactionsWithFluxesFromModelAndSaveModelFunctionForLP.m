function [ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelFunctionForLP(taskModel, solution,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j)


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
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

b = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

pattern = '(\d+\s)?MAM\d+\[[a-z]\]';
pattern3 = 'MAM\d+\[[a-z]\]';
pattern2 = '\[(.*?)\]';

rxnFormulas = cell(length(b),1);
for bIterator = 1:length(b)
    c = string(b(bIterator));
    fullFormula = "";
    splitString = split(c, '->');
    metabolitesBeforeForIdx = regexp(splitString{1}, pattern3, 'match');
    metabolitesBefore = regexp(splitString{1}, pattern, 'match');
    beforeMetsArray = cell(length(metabolitesBefore),1);
    for cIterator = 1: length(metabolitesBefore)
        idx = findMetIDs(taskModel,metabolitesBeforeForIdx(cIterator));
        metName = taskModel.metNames(idx);
        numberStr = regexprep(metabolitesBefore{cIterator}, '(\d+\s)?MAM\d+\[[a-z]\]', '$1');
        numberOfMetabolites = str2double(numberStr);
        extension = regexp(metabolitesBefore{cIterator}, pattern2, 'tokens');
        if  isnan(numberOfMetabolites)  
             newName = sprintf("%s[%s]", string(metName), string(extension));
        elseif mod(numberOfMetabolites, 1) == 0
            newName = sprintf("%.0f %s[%s]",numberOfMetabolites , string(metName), string(extension));

        else
           newName = sprintf("%.2f %s[%s]",numberOfMetabolites , string(metName), string(extension));
        end
       
        beforeMetsArray(cIterator) = cellstr(newName);
    end

    metabolitesAfter = regexp(splitString{2}, pattern, 'match');
    metabolitesAfterForIDX = regexp(splitString{2}, pattern3, 'match');
    afterMetsArray = cell(length(metabolitesAfter),1);
    for cIterator = 1: length(metabolitesAfter)
        idx = findMetIDs(taskModel,metabolitesAfterForIDX(cIterator));
        metName = taskModel.metNames(idx);
        numberStr = regexprep(metabolitesAfter{cIterator}, '(\d+\s)?MAM\d+\[[a-z]\]', '$1');
        numberOfMetabolites = str2double(numberStr);
        extension = regexp(metabolitesAfter{cIterator}, pattern2, 'tokens');
        if  isnan(numberOfMetabolites)  
             newName = sprintf("%s[%s]", string(metName), string(extension));
        elseif mod(numberOfMetabolites, 1) == 0
            newName = sprintf("%.0f %s[%s]",numberOfMetabolites , string(metName), string(extension));
           
        else
            newName = sprintf("%.2f %s[%s]",numberOfMetabolites , string(metName), string(extension));
        end      
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

save(solutionName, "solution");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;
scoreLP = score;

end