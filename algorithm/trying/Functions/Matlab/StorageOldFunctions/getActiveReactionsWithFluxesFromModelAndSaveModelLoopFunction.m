function ActiveReactionFluxes = getActiveReactionsWithFluxesFromModelAndSaveModelLoopFunction(taskModel, solution,expValueMeanAdjustedTask,score ,taskSpecificMILPFiltered, i, osenseStr , irrev)


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

dispName =  sprintf("Task%iScored%.3f",i, score);

reactionLength = size(taskSpecificMILPFiltered.S,2);


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

newVector = zeros(size(solution.cont(1:reactionLength))); 
newVector(abs(solution.cont(1:reactionLength)) > 0.99)  = 1;
selectedRows = newVector == 1;

sizeDiff = length(taskModel.rxns)-reactionLength;
selectedRows1 = [selectedRows;zeros(sizeDiff,1)];
selectedRows1 = selectedRows1==1;
printFlag = false;


tmodel = taskModel;
tmodel.mets = tmodel.metNames;
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows1),printFlag );

yiFlagsCombined = ones(size(reactionLength,1),1);
%newCellArray = [num2cell(yiFlagsCombined(selectedRows)), num2cell(solution.cont(selectedRows)), num2cell(expValueMeanAdjustedTask(selectedRows)), taskModel.rxns(selectedRows1), taskModel.rxnNames(selectedRows1),rxnFormulas];

 num2cell(solution.cont(selectedRows))
taskModel.rxns(selectedRows1)
 taskModel.rxnNames(selectedRows1)
 
newCellArray = [num2cell(solution.cont(selectedRows)), taskModel.rxns(selectedRows1), taskModel.rxnNames(selectedRows1),rxnFormulas];


removeRows = newVector == 0;

newModel = removeRxns(taskModel, taskModel.rxns(removeRows));

solutionName = string("lastValidSolutionLoopFunction"+dispName+".mat");
newModelName = string("newModelLoopFunction"+dispName+".mat");
fluxSolutionName = string("fluxTableLoopFunction"+dispName+".mat");

save(solutionName, "solution");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;

end