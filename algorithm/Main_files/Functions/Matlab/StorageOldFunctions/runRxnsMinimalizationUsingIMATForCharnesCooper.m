function [newSolution, prunedModel, originalModel, adjustedExpressionData] = runRxnsMinimalizationUsingIMATForCharnesCooper(taskModel, solution, expValueMeanAdjustedTask,reactionLength, epsilon)


% Original model
x = solution.cont(1:reactionLength);
rxnRemList = taskModel.rxns(abs(x) < epsilon-0.1);
originalModel = removeRxns(taskModel,rxnRemList);

% Find indices for expression
expressionIndices1 = expValueMeanAdjustedTask ~= -1;
expressionIndices = [expressionIndices1;expressionIndices1];
logicalIndicesExpression1 = find((solution.int(1:2*reactionLength) == 1) & (expressionIndices == 1));
logicalIndicesExpression1(logicalIndicesExpression1>reactionLength) = logicalIndicesExpression1(logicalIndicesExpression1>reactionLength)-reactionLength;
logicalIndicesExpressionNotInSolution =  find((solution.int(1:2*reactionLength) == 0) & (expressionIndices == 1));
logicalIndicesExpressionNotInSolution(logicalIndicesExpressionNotInSolution>reactionLength) = (logicalIndicesExpressionNotInSolution(logicalIndicesExpressionNotInSolution>reactionLength)-reactionLength);
logicalIndicesExpressionNotInSolution = unique(logicalIndicesExpressionNotInSolution);

% TESTING PURPOSES
% a = solution.cont(logicalIndicesExpression1);
% exprRxns = taskModel.rxns(logicalIndicesExpression1);
% exprReacitonNames = taskModel.rxnNames(logicalIndicesExpression1);
% exprExpr = expValueMeanAdjustedTask(logicalIndicesExpression1);
% ubs = taskModel.ub(logicalIndicesExpression1);
% lbs =  taskModel.lb(logicalIndicesExpression1);

% Adapt all expressionRxns other than those in solution to have LB/UB 0
%set all expressionsRxns in solution to have LB of 1
taskModel1 = taskModel;
negativeLBIdxExpression = taskModel.lb(logicalIndicesExpression1)<0;

taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpression1(negativeLBIdxExpression)), -epsilon, 'u');
taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpression1(negativeLBIdxExpression)), -70, 'l');

taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpression1(negativeLBIdxExpression==0)), epsilon, 'l');
taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpression1(negativeLBIdxExpression==0)), 70, 'u');

taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpressionNotInSolution), 6, 'u');
taskModel1 = changeRxnBounds(taskModel1, taskModel1.rxns(logicalIndicesExpressionNotInSolution), -6, 'l');



% Create thresholds so that all values are part of RLindex in IMAT
arbitrarilyHighValue = 100000;
threshold_lb = arbitrarilyHighValue;
threshold_ub = arbitrarilyHighValue+1;

% Modified IMAT version that includes ORs in the optimization
% Minimizes ORs, please note: the prunedModel contains modified LBs for
% expression carrying reactions
[prunedModel,newSolution] = iMAT_AdaptedForORRemoval(taskModel1, expValueMeanAdjustedTask, threshold_lb, threshold_ub);
findRxnIDs(taskModel1, prunedModel.rxns);
adjustedExpressionData = expValueMeanAdjustedTask (findRxnIDs(taskModel1, prunedModel.rxns));

if length(logicalIndicesExpression1) ~= length(find((newSolution.int(1:reactionLength) == 0) & (expressionIndices1 == 1)))
    disp("in IMAT Optimization to minimize ORs: Length of expression reactions in originalModel vs prunedModel are not the same, did not manage to protect reactions")
elseif length(originalModel.rxns)< length(prunedModel.rxns)
    disp("in IMAT Optimization to minimize ORs: Original model contains less reactions after IMAT pruning, this should not happen!")
elseif newSolution.stat ~= 1
    disp("in IMAT Optimization to minimize ORs: Did not get a feasible solution! This has to do with the bounds, unsure how to fix")
else
    fprintf("\nAmount of Expression carrying reactions in original and pruned model: %i -- %i \n", length(logicalIndicesExpression1),length(find((newSolution.int(1:reactionLength) == 0) & (expressionIndices == 1))) )
    % disp(length(logicalIndicesExpression1))
    % disp(length(find((newSolution.int(1:reactionLength) == 0) & (expressionIndices == 1))));
end

end