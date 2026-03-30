function [ActiveReactionFluxes,reactionsToKeepLP,scoreLP] = runLPProblemOneOverGeneExpression_LPiterationsTest(taskModel,expValueMeanAdjustedTask,orphanHandling,osenseStr,TRORAdjust,iterationTest)

taskModel.osense = 1;
taskModel.osenseStr = "min";
taskModel.c = zeros(size(expValueMeanAdjustedTask)); 
if orphanHandling == "one"
    taskModel.c(expValueMeanAdjustedTask == -1) = 1 ./ iterationTest;
elseif orphanHandling == "median"
    orphanHandleParameter = median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0));
    taskModel.c(expValueMeanAdjustedTask == -1) = 1 ./ orphanHandleParameter;
elseif orphanHandling == "zero"
    orphanHandleParameter = 0.0000001;
    taskModel.c(expValueMeanAdjustedTask == -1) = 1 .* orphanHandleParameter;
else
    disp("OrphanHandling should be 'one', 'median', or 'zero': Stopping function")
    return
end

taskModel.c(expValueMeanAdjustedTask ~= -1) = 1 ./ expValueMeanAdjustedTask(expValueMeanAdjustedTask ~= -1); % Already existed

if TRORAdjust
    transportRxnIDs = findRealTransportReactions_AllTransport_EmptyGRRules(taskModel);
    taskModel.c(transportRxnIDs) = 0.0000001;
end

solutionLP = optimizeCbModel(taskModel,osenseStr);

newVector = zeros(size(solutionLP.v)); 
newVector(abs(solutionLP.v) > 0.0001) = 1;
reactionsToKeepLP = newVector == 1;

[ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelForLP_TROR(taskModel,solutionLP,expValueMeanAdjustedTask);

end