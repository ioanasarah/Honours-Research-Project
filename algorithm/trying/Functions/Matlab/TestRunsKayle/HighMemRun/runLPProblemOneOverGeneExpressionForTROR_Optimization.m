function [ActiveReactionFluxes, reactionsToKeepLP,scoreLP] = runLPProblemOneOverGeneExpressionForTROR_Optimization(taskModel, expValueMeanAdjustedTask, i, orphanHandling,j,irrev,osenseStr,TRORAdjust)


taskModel.osense = 1;
taskModel.osenseStr = "min";
taskModel.c = zeros(size(expValueMeanAdjustedTask)); 
if orphanHandling == "one"
    orphanHandleParameter = 1;
    taskModel.c(expValueMeanAdjustedTask == -1) = 1;
elseif orphanHandling == "median"
    orphanHandleParameter = median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0));
    %disp(orphanHandleParameter)
    taskModel.c(expValueMeanAdjustedTask == -1) = 1 ./ orphanHandleParameter;
elseif orphanHandling == "zero"
    orphanHandleParameter = 0.0000001;
    %disp(median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0)))
    taskModel.c(expValueMeanAdjustedTask == -1) = 1.*orphanHandleParameter;
else
    disp("OrphanHandling should be 'one', 'median', or 'zero': Stopping function")
    return
end


aboveZero = taskModel.S > 0;
sumValuesAboveZero = transpose(sum(taskModel.S .* aboveZero));
belowZero = taskModel.S < 0;
sumValuesBelowZero = abs(transpose(sum(taskModel.S .* belowZero)));

perReactionModifierApproach1 = sumValuesAboveZero ./sumValuesBelowZero;
perReactionModifierApproach2 = sumValuesBelowZero ./sumValuesAboveZero;

perReactionModifierApproach1(perReactionModifierApproach1 == Inf) = 1;
perReactionModifierApproach2(perReactionModifierApproach2 == Inf) = 1;

perReactionModifierApproach3 = max(perReactionModifierApproach1, 1);
perReactionModifierApproach4 = max(perReactionModifierApproach2, 1);


taskModel.c(expValueMeanAdjustedTask ~= -1) = 1 ./ expValueMeanAdjustedTask(expValueMeanAdjustedTask ~= -1); % Already existed
%taskModel.c = full(taskModel.c .*perReactionModifierApproach4); 
if TRORAdjust
    transportRxnIDs = findRealTransportReactions_AllTransport_EmptyGRRules(taskModel);
    taskModel.c(transportRxnIDs) =0.0000001;
end
score = 0;
solutionLP = optimizeCbModel(taskModel,osenseStr);

newVector = zeros(size(solutionLP.v)); 
newVector(abs(solutionLP.v) > 0.0001 )= 1;
reactionsToKeepLP = newVector == 1;

[ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelForLP_TROR(taskModel, solutionLP,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j);
end