function [taskSpecificMILPFiltered,expValueMeanAdjustedTaskFiltered] = adjustTaskSpecificMILPForCCFiltering(taskSpecificMILPFunction,reactionsNotPartOfTask,taskModelOldCore,expValueMeanAdjustedTask)

%taskSpecificMILPFunction = taskSpecificMILP;
%spy(taskSpecificMILP.A)
tic
expValueMeanAdjustedTaskFiltered = expValueMeanAdjustedTask;
expValueMeanAdjustedTaskFiltered(reactionsNotPartOfTask==1) = -1; 

SizeA = zeros(size(taskModelOldCore.S,1),1);
extendedReactionsNotPartOfTaskRows = [SizeA;reactionsNotPartOfTask;reactionsNotPartOfTask;reactionsNotPartOfTask;reactionsNotPartOfTask;reactionsNotPartOfTask];

taskSpecificMILPFunction.b(extendedReactionsNotPartOfTaskRows==1) = 0;

extendedReactionsNotPartOfTaskColumns = [reactionsNotPartOfTask;reactionsNotPartOfTask;reactionsNotPartOfTask;reactionsNotPartOfTask];

taskSpecificMILPFunction.lb(extendedReactionsNotPartOfTaskColumns==1) = 0;
taskSpecificMILPFunction.ub(extendedReactionsNotPartOfTaskColumns==1) = 0;
taskSpecificMILPFunction.c(extendedReactionsNotPartOfTaskColumns==1)  = 0;

taskSpecificMILPFunction.A(size(taskModelOldCore.S,1)+1:length(extendedReactionsNotPartOfTaskRows), extendedReactionsNotPartOfTaskColumns ==1) =0;
taskSpecificMILPFunction.A(length(extendedReactionsNotPartOfTaskRows)+1,extendedReactionsNotPartOfTaskColumns==1) = 0;


% nonEmptyRows = any(taskSpecificMILPFunction.A(1:length(extendedReactionsNotPartOfTaskRows),1:length(extendedReactionsNotPartOfTaskColumns)), 2);
% nonEmptyCols =  any(taskSpecificMILPFunction.A(1:length(extendedReactionsNotPartOfTaskRows),1:length(extendedReactionsNotPartOfTaskColumns)), 1);

nonEmptyRows = any(taskSpecificMILPFunction.A, 2);
nonEmptyCols =  any(taskSpecificMILPFunction.A, 1);

expValueMeanAdjustedTaskFiltered = expValueMeanAdjustedTaskFiltered(reactionsNotPartOfTask==0);

A_cleaned = taskSpecificMILPFunction.A(nonEmptyRows,nonEmptyCols);

taskSpecificMILPFunction.lb = taskSpecificMILPFunction.lb(nonEmptyCols);
taskSpecificMILPFunction.ub = taskSpecificMILPFunction.ub(nonEmptyCols);
taskSpecificMILPFunction.c = taskSpecificMILPFunction.c(nonEmptyCols);
taskSpecificMILPFunction.vartype = taskSpecificMILPFunction.vartype(nonEmptyCols);

taskSpecificMILPFunction.b = taskSpecificMILPFunction.b(nonEmptyRows);
taskSpecificMILPFunction.csense = taskSpecificMILPFunction.csense(nonEmptyRows);

emptyRows = find(~nonEmptyRows);
emptyCols = find(~nonEmptyCols);

taskSpecificMILPFunction.A = A_cleaned;
sizeOldSRow = size(taskModelOldCore.S,1);
sizeOldSColumn = size(taskModelOldCore.S,2);
nonEmptyRowsS = nonEmptyRows(1:sizeOldSRow);
nonEmptyColsS = nonEmptyCols(1:sizeOldSColumn);

taskSpecificMILPFunction.S = taskModelOldCore.S(nonEmptyRowsS,nonEmptyColsS);
taskSpecificMILPFunction.rxns = taskModelOldCore.rxns(nonEmptyColsS);


taskSpecificMILPFiltered = taskSpecificMILPFunction;
spy(taskSpecificMILPFiltered.A)
toc


end