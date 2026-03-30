function [highestK,optionsSaved] = findHighestKIterativeLPs_Vi_Times_Gi(taskModel,expValueMeanAdjustedTask,osenseStr)

options = ["one","median","zero"];
optionsTROR = [false,true];
score = 0;

for ii = 1:length(options)
    currentOption = options(ii);
    for jj = 1:length(optionsTROR)
        currentTror = optionsTROR(jj);
        [~,~,scoreLP] = runLPProblemOneOverGeneExpressionViTimesGi(taskModel,expValueMeanAdjustedTask,currentOption,osenseStr,currentTror);
        if scoreLP > score
            score = scoreLP;
            optionsSaved = {options(ii),optionsTROR(jj)};
        end
    end
end

highestK = score;

end