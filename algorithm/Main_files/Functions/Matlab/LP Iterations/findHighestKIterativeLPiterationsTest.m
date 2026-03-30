function [highestK,optionsSaved] = findHighestKIterativeLPiterationsTest(taskModel,expValueMeanAdjustedTask,osenseStr)

options = ["one","median","zero"];
optionsTROR = [false,true];
score = 0;
iterationTest = 0;

for ii = 1:length(options)
    currentOption = options(ii);
    for jj = 1:length(optionsTROR)
        currentTror = optionsTROR(jj);

        if currentOption == "one"
            while iterationTest < 10
                iterationTest = iterationTest + 0.01;
                [~,~,scoreLP] = runLPProblemOneOverGeneExpression_LPiterationsTest(taskModel,expValueMeanAdjustedTask,currentOption,osenseStr,currentTror,iterationTest);
                fprintf("option: %s\n",currentOption)
                fprintf("TROR option: %i\n",currentTror)
                fprintf("iterationTest: 1/%f\n",iterationTest)
                fprintf("LP score:%f\n",scoreLP)
                fprintf("\n")

                if scoreLP > score
                    score = scoreLP;
                    optionsSaved = {options(ii),optionsTROR(jj)};
                end
            end
        else
            [~,~,scoreLP] = runLPProblemOneOverGeneExpression_LPiterationsTest(taskModel,expValueMeanAdjustedTask,currentOption,osenseStr,currentTror,iterationTest);
            fprintf("option: %s\n",currentOption)
            fprintf("TROR option: %i\n",currentTror)
            fprintf("LP score:%f\n",scoreLP)
            fprintf("\n")

            if scoreLP > score
                score = scoreLP;
                optionsSaved = {options(ii),optionsTROR(jj)};
            end
        end
    end
end

highestK = score;

end