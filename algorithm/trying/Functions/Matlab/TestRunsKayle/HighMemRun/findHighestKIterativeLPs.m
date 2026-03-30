function [highestK, optionsSaved] = findHighestKIterativeLPs(taskModel, expValueMeanAdjustedTask, i,j,irrev,osenseStr)

options = ["one","median", "zero"];
optionsTROR = [false,true];
score = 0;
%optionsSaved = "Unfeasible LP";

for ii=1:length(options)
    currentOption = options(ii);
    for jj = 1:length(optionsTROR)
        currentTror = optionsTROR(jj);
        [ActiveReactionFluxes, reactionsToKeepLP,scoreLP] =runLPProblemOneOverGeneExpressionForTROR_Optimization(taskModel, expValueMeanAdjustedTask, i,currentOption,j,irrev,osenseStr,currentTror);
        if scoreLP > score
            score = scoreLP;
            optionsSaved = {options(ii),optionsTROR(jj)};
        end
    end
end

highestK = score

end