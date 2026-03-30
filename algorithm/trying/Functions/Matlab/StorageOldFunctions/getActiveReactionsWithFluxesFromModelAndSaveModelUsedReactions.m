function [ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelUsedReactions(taskModel, solution,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j)

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
newVector(abs(solution.v) > 0.1 )= 1;
selectedRows = newVector == 1;
printFlag = false;
tmodel = taskModel;
tmodel.mets = tmodel.metNames;
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );


newCellArray = [num2cell(solution.v(selectedRows)), num2cell(expValueMeanAdjustedTask(selectedRows)), num2cell(taskModel.c(selectedRows)), taskModel.rxns(selectedRows), taskModel.rxnNames(selectedRows), rxnFormulas];

validValues = expValueMeanAdjustedTask(selectedRows) ~= -1;
selectedExpression = expValueMeanAdjustedTask(selectedRows);
selectedFluxes = solution.v(selectedRows);

score = sum(selectedExpression(validValues)) / length(selectedFluxes(validValues));
dispName =  sprintf("UsedReactionsTask%iSample%iScored%.3f%s%s%s",i,j, score, osenseStr,irrevStr);

removeRows = newVector == 0;

newModel = removeRxns(taskModel, taskModel.rxns(removeRows));

solutionName = string("lastValidSolutionUsedReactions"+dispName+".mat");
newModelName = string("newModelUsedReactions"+dispName+".mat");
fluxSolutionName = string("fluxTableUsedReactions"+dispName+".mat");

save(solutionName, "solution");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;
scoreLP = score;

end