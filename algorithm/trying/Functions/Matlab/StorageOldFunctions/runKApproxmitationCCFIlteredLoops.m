function [solution,score,reactionsToKeep] = runKApproxmitationCCFIlteredLoops (taskModelOldCore, expValueMeanAdjustedTask, epsilon)

timerval = tic;
kPercentageConstant = 0.6;
ktest = 0;    
kTolerance = 0.2; %reduces computation time when set higher, probably shouldn't be too high
kfeasibletest = 0;
kInfeasibletest = max(expValueMeanAdjustedTask);
counter = 0;
disp("create MILPproblem")
MILPproblemTaskRevCCFiltered = createKMILPProblemAverageExpressionLoops(taskModelOldCore, expValueMeanAdjustedTask, epsilon,1 );
MILPproblem = MILPproblemTaskRevCCFiltered;
% While loop

ktest = kInfeasibletest*kPercentageConstant;
while abs(kfeasibletest - kInfeasibletest) > kTolerance
    tic
    %ktest = 2;
    disp({kfeasibletest;kInfeasibletest;ktest})
    
    % Adjust MILPproblem to change to new ktest
    %MILPproblem = createKMILPProblemAverageExpression(taskModel, expValueMeanAdjusted, epsilon, rev2irrev, ktest);
    MILPproblem = adaptKMILPProblemAverageExpression(MILPproblem, taskModelOldCore, expValueMeanAdjustedTask, ktest);
    % Calculate solution
    solution = solveCobraMILP(MILPproblem);
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
    disp(length(find(abs(solution.cont)>=epsilon)))
    toc
end 

% incase narrowing ends on invalid solution due to rounding (?)
solution = lastValidSolution;

% remove all reactions not in solution
score = lastValidK;
%map to reversible model


reactionLength = size(taskSpecificMILP.S,2);
reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= epsilon;
toc(timerval) 
end