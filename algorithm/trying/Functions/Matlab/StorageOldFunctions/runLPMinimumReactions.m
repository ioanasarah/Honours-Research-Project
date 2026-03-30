function [ActiveReactionFluxes, reactionsToKeepLP,scoreLP] = runLPMinimumReactions(taskModel, expValueMeanAdjustedTask, i, orphanHandling,j,irrev,osenseStr)

taskModel.osense = 1;
taskModel.osenseStr = "min";
taskModel.c = ones(size(expValueMeanAdjustedTask)); 
% if orphanHandling == "one"
%     orphanHandleParameter = 1;
%     taskModel.c(expValueMeanAdjustedTask == -1) = 1;
% elseif orphanHandling == "median"
%     orphanHandleParameter = median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0));
%     taskModel.c(expValueMeanAdjustedTask == -1) = 1 ./ orphanHandleParameter;
% elseif orphanHandling == "zero"
%     orphanHandleParameter = 0;
%     taskModel.c(expValueMeanAdjustedTask == -1) = 1.*orphanHandleParameter;
% else
%     disp("OrphanHandling should be 'one', 'median', or 'zero': Stopping function")
%     return
% end
% 
% 
% 
% taskModel.c(expValueMeanAdjustedTask ~= -1) = 1 ./ expValueMeanAdjustedTask(expValueMeanAdjustedTask ~= -1);

solutionLP = optimizeCbModel(taskModel,osenseStr);
score = 0;

newVector = zeros(size(solutionLP.v)); 
newVector(abs(solutionLP.v) > 0.1 )= 1;
reactionsToKeepLP = newVector == 1;

 [ActiveReactionFluxes,scoreLP] = getActiveReactionsWithFluxesFromModelAndSaveModelMiniReactions(taskModel, solutionLP,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j);
end