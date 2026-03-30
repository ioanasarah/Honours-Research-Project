function [AnswerSol,ReactionSolutions,activeReactions,AverageGeneExpressionAverageRxns,errorsInTask] = calculateSolutionAnswer(solution,model,epsilon, expressionRxns)

AnswerSol = solution.obj;




reactionLength = length(model.rxns);

ReactionSolutions = zeros(reactionLength, 8);
ReactionSolutions(1:reactionLength,1) = solution.int(1:reactionLength); %forward y
ReactionSolutions(1:reactionLength,2) = solution.int(reactionLength+1:reactionLength+reactionLength); %backward y
ReactionSolutions(1:reactionLength,4) = solution.cont(1*reactionLength+1:2*reactionLength); %forward Z
ReactionSolutions(1:reactionLength,5) = solution.cont(2*reactionLength+1:3*reactionLength); %backward z
ReactionSolutions(1:reactionLength,7) = solution.cont(1:reactionLength); %flux
ReactionSolutions(1:reactionLength,8) = expressionRxns; % reaction data

% Check if forward/backward reactions are both active (should not)
% set 1 for forward, -1 for backward
for i = 1:size(ReactionSolutions,1)
    if ReactionSolutions(i,1) == 1 && ReactionSolutions(i,2) == 1
        fprintf("rxn %s has both forward and backward reactions active", i)
    end

    if ReactionSolutions(i,1) == 1
        ReactionSolutions(i,3) = 1;
    elseif ReactionSolutions(i,2) == 1
        ReactionSolutions(i,3) = -1;
    end

    if ReactionSolutions(i,4) ~= 0
        ReactionSolutions(i,6) = ReactionSolutions(i,4);
    elseif ReactionSolutions(i,5) ~= 0
        ReactionSolutions(i,6) = -ReactionSolutions(i,5);
    end
end


%Test the data
ReactionSolutions = num2cell(ReactionSolutions);
ReactionSolutions = [ReactionSolutions, model.rxns];
ReactionSolutions(:, [1,2,4,5]) = [];

columnNames = {'Yi_F/R', 'Zi_F/R', 'Flux', 'Expression', 'rxnName'};
ReactionSolutions = cell2table(ReactionSolutions, 'VariableNames', columnNames);

%count active reactions
activeReactions = ReactionSolutions(round(abs([ReactionSolutions{:, 3}]),6) >= epsilon, :);
activeReactionsVcount = size(ReactionSolutions(abs([ReactionSolutions{:, 3}]) >= epsilon, :),1);
activeReactionsYcount = nnz(ReactionSolutions{:,1});
sumAverageRxns = sum(activeReactions(activeReactions{:,4} ~=-1  ,4));
AverageGeneExpressionAverageRxns = sumAverageRxns{1,1}/size(activeReactions(activeReactions{:,4} ~=-1,4),1);

if round(AverageGeneExpressionAverageRxns,3)  ~= round(AnswerSol,3) || activeReactionsVcount ~= activeReactionsYcount
    errorsInTask = true;
else
    errorsInTask = false;
end



%disp({activeReactionsVcount,sumAverageRxns{1,1}})



