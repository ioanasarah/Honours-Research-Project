function ActiveReactionFluxes = getActiveReactionsWithFluxesFromModelAndSaveModelCooper(taskModelCooperIrrev, solutionFractionalIrrev,expValueMeanAdjustedTaskCharnesCooperIrrev,score , i, osenseStr , irrev,j)


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

dispName =  sprintf("Task%iSample%iScored%.3f%s",i,j, score);

reactionLength = length(taskModelCooperIrrev.rxns);

if ~isstruct(solutionFractionalIrrev)
     ActiveReactionFluxes=0;
else
yiFlagsPlus =  solutionFractionalIrrev.int(1:reactionLength);
yiFlagsMinus = solutionFractionalIrrev.int(reactionLength+1:2*reactionLength);
yiFlagsCombined = zeros(reactionLength,1);

for j = 1:reactionLength
    if yiFlagsPlus(j) == 1 && yiFlagsMinus(j) == 1
        fprintf("rxn %s has both forward and backward reactions active\n", j)
    end

    if yiFlagsPlus(j) == 1
        yiFlagsCombined(j) = 1;
    elseif yiFlagsMinus(j) == 1
        yiFlagsCombined(j) = -1;
    end
end

newVector = zeros(reactionLength,1); 
newVector(abs(solutionFractionalIrrev.cont(1:reactionLength)) > 0.99  | yiFlagsCombined ~= 0| yiFlagsMinus ~=0) = 1;
selectedRows = newVector == 1;

printFlag = false;


tmodel = taskModelCooperIrrev;
tmodel.mets = tmodel.metNames;
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

newCellArray = [num2cell(yiFlagsCombined(selectedRows)), num2cell(solutionFractionalIrrev.cont(selectedRows)), num2cell(expValueMeanAdjustedTaskCharnesCooperIrrev(selectedRows)), taskModelCooperIrrev.rxns(selectedRows), taskModelCooperIrrev.rxnNames(selectedRows),rxnFormulas];


removeRows = newVector == 0;

newModel = removeRxns(taskModelCooperIrrev, taskModelCooperIrrev.rxns(removeRows));

solutionName = string("lastValidSolutionCooper"+dispName+".mat");
newModelName = string("newModelCooper"+dispName+".mat");
fluxSolutionName = string("fluxTableCooper"+dispName+".mat");

save(solutionName, "solutionFractionalIrrev");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;


end