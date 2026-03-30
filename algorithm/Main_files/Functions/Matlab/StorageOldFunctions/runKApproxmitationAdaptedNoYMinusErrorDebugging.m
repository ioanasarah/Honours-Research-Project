function [solution,score,reactionsToKeep] = runKApproxmitationAdaptedNoYMinusErrorDebugging (taskModel, expValueMeanAdjustedTask, epsilon, rev2irrev)

timerval = tic;
kPercentageConstant = 0.5;
ktest = 0;    
kTolerance = 0.2; %reduces computation time when set higher, probably shouldn't be too high
kfeasibletest = 0;
kInfeasibletest = max(expValueMeanAdjustedTask);
counter = 0;
disp("create MILPproblem")

MILPproblem = createKMILPProblemAverageExpressionAdaptedNoYiMinusDebug(taskModel, expValueMeanAdjustedTask, epsilon, rev2irrev, ktest);
varargin = parseSolverParameters('MILP');
% varargin.feasTol = 10^-2;
% varargin.intTol = 10^-2;
% varargin.relMipGapTol = 10^-2;
% varargin.absMipGapTol = 10^-2;
% varargin.optTol = 10^-2;
% While loop

ktest = kInfeasibletest*kPercentageConstant;
while abs(kfeasibletest - kInfeasibletest) > kTolerance
    disp({kfeasibletest;kInfeasibletest;ktest})
    % Adjust MILPproblem to change to new ktest
    %MILPproblem = createKMILPProblemAverageExpression(taskModel, expValueMeanAdjusted, epsilon, rev2irrev, ktest);
    MILPproblem = adaptKMILPProblemAverageExpressionAdaptedNoYiMinus(MILPproblem, taskModel, expValueMeanAdjustedTask, ktest);
    % Calculate solution

    solution = solveCobraMILP(MILPproblem,varargin);
    counter = counter +1;

    disp(solution.stat)
    % if solution was feasible, set kfeasible to ktest
    
    if solution.stat ==1
        kfeasibletest = ktest;
        lastValidSolution = solution;
        lastValidK=ktest;
        kDiff = kInfeasibletest-kfeasibletest;
        ktest = kfeasibletest + (kPercentageConstant*kDiff);
    % if solution was infeasible, set kinfeasible to ktest
    elseif solution.stat ~= 1 
        kInfeasibletest = ktest;
        kDiff = kInfeasibletest-kfeasibletest;
        ktest = kInfeasibletest - ((1-kPercentageConstant)*kDiff);
    else
        fprintf("Solution stat was not 1 or -1, namely: %s", "solution.stat")
    end
    disp(length(find(solution.cont>=epsilon)))
end 

% incase narrowing ends on invalid solution due to rounding (?)
if exist('lastValidSolution')
    solution = lastValidSolution;
    score = lastValidK;
    reactionLength = length(taskModel.rxns);
    reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= epsilon;
else
    solution = 0;
    score = 0;
    reactionLength = 0;
    reactionsToKeep = 0;
end


% remove all reactions not in solution
% score = lastValidK;
%map to reversible model


% reactionLength = length(taskModel.rxns);
% reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= epsilon;
%disp(length(reactionsToKeep))
toc(timerval) 
end