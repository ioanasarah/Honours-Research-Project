function ActiveReactionFluxes = getActiveReactionsWithFluxesFromModelAndSaveModelMinYKapprox(taskModel, solution,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j)


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

dispName =  sprintf("Task%iSample%iScored%.3f%s%s%s",i,j, score, osenseStr,irrevStr);

reactionLength = length(taskModel.rxns);
disp(reactionLength)
disp(length(solution.cont))

if ~isstruct(solution)
     ActiveReactionFluxes=0;
else
yiFlagsPlus =  solution.int(1:reactionLength);
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
newVector(abs(solution.cont(1:reactionLength)) > 0.99 | yiFlagsPlus ~= 0) = 1;
selectedRows = newVector == 1;

printFlag = false;


tmodel = taskModel;
tmodel.mets = tmodel.metNames;
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

newCellArray = [num2cell(yiFlagsPlus(selectedRows)), num2cell(solution.cont(selectedRows)), num2cell(expValueMeanAdjustedTask(selectedRows)), taskModel.rxns(selectedRows), taskModel.rxnNames(selectedRows),rxnFormulas];


removeRows = newVector == 0;

newModel = removeRxns(taskModel, taskModel.rxns(removeRows));

solutionName = string("lastValidSolutionMILPMinYKApprox"+dispName+".mat");
newModelName = string("newModelMILPMinYKApprox"+dispName+".mat");
fluxSolutionName = string("fluxTableMILPMinYKApprox"+dispName+".mat");

save(solutionName, "solution");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;


end