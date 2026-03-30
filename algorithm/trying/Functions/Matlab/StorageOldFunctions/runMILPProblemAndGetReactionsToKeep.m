function [reactionsToKeep, score, MILPSolution] = runMILPProblemAndGetReactionsToKeep( MILPproblem, taskModel , epsilon)

% Possibly adjuts tolerances
varargin = parseSolverParameters('MILP');
adjustTol = true;
if adjustTol 
    varargin.intTol=1.0E-12;
    varargin.relMipGapTol=1.0E-12;
    varargin.absMipGapTol=1.0E-12;
    varargin.feasTol=1.0E-6; % This tolerance determines when a solution is
    % treated as feasible. In MILP, constraints and variable bounds may not
    % be satisfied exactly due to numerical errors. The feasTol parameter 
    % specifies the maximum violation allowed for a constraint or variable
    % bound to be considered as satisfied. If the violation of a constraint's
    % activity or a variable's bound is less in absolute magnitude than feasTol,
    % it is treated as satisfied. A smaller feasTol value indicates a stricter
    % tolerance for feasibility.
    varargin.optTol=1.0E-6;
end

% Perform test
solution = solveCobraMILP(MILPproblem,varargin);
% 
% % Get int solutions for checking
% ReactionSolutions = zeros(95, 8);
% ReactionSolutions(1:reactionLength,1) = solution.int(1:reactionLength); %forward y
% ReactionSolutions(1:reactionLength,2) = solution.int(reactionLength+1:reactionLength+reactionLength); %backward y
% ReactionSolutions(1:reactionLength,4) = solution.cont(1*reactionLength+1:2*reactionLength); %forward Z
% ReactionSolutions(1:reactionLength,5) = solution.cont(2*reactionLength+1:3*reactionLength); %backward z
% ReactionSolutions(1:reactionLength,7) = solution.cont(1:reactionLength); %flux
% ReactionSolutions(1:reactionLength,8) = expressionRxns; % reaction data
% 
% % Check if forward/backward reactions are both active (should not)
% % set 1 for forward, -1 for backward
% for i = 1:size(ReactionSolutions,1)
%     if ReactionSolutions(i,1) == 1 && ReactionSolutions(i,2) == 1 
%         fprintf("rxn %s has both forward and backward reactions active", i)
%     end
% 
%     if ReactionSolutions(i,1) == 1
%         ReactionSolutions(i,3) = 1;
%     elseif ReactionSolutions(i,2) == 1
%         ReactionSolutions(i,3) = -1;
%     end
% 
%     if ReactionSolutions(i,4) ~= 0
%         ReactionSolutions(i,6) = ReactionSolutions(i,4);
%     elseif ReactionSolutions(i,5) ~= 0
%         ReactionSolutions(i,6) = -ReactionSolutions(i,5);
%     end
% end
% 
% 
% %Test the data 
% ReactionSolutions = num2cell(ReactionSolutions);
% ReactionSolutions = [ReactionSolutions, model.rxns];
% ReactionSolutions(:, [1,2,4,5]) = [];
% 
% activeRxnsNum =  nnz([ReactionSolutions{:, 3}]);
% AverageGeneExpression = solution.obj/activeRxnsNum;
% 
% columnNames = {'Yi_F/R', 'Zi_F/R', 'Flux', 'Expression', 'rxnName'};
% ReactionSolutions = cell2table(ReactionSolutions, 'VariableNames', columnNames);
% 
% % tol = epsilon-(epsilon/10);
% activeReactions = ReactionSolutions(abs([ReactionSolutions{:, 3}]) >= epsilon, :);

reactionLength = length(taskModel.rxns);
reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= epsilon;

MILPSolution = solution;
score = solution.obj;

end