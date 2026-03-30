function [solution,score,reactionsToKeep] = runKApproxmitationIrrevTROR_Adjustment (taskModel, expValueMeanAdjustedTask, epsilon, rev2irrev, alpha, beta, i,j,osenseStr,irrev,orphanHandling)

% alpha = 0.5; % even addition of the TROR adjustment as regular obj funciton at 0.5
%              % alpha gives the weights to the TROR adjustment of the obj
%              % function and to the original objective function by setting
% 
% beta = 1;  % total addition to max score if all oprhan transport reactions would be 0
reactionLength = length(taskModel.rxns);
transportRxnIDs = findRealTransportReactions_AllTransport_EmptyGRRules(taskModel); 
disp(length(transportRxnIDs))
timerval = tic;
kPercentageConstant = 0.5;
ktest = 0;    
kTolerance = 0.1; %reduces computation time when set higher, probably shouldn't be too high

SMatrix = taskModel.S;
%bestK = findHighestKIterativeLPs(taskModel,expValueMeanAdjustedTask);
%[ActiveReactionFluxesLP, reactionsToKeepLP,scoreLP] = runLPProblemOneOverGeneExpression(taskModel, expValueMeanAdjustedTask, i,orphanHandling,j,irrev,osenseStr);


[kfeasibletest, optionsSaved] =  findHighestKIterativeLPs(taskModel, expValueMeanAdjustedTask, i,j,irrev,osenseStr);
disp(optionsSaved)

%kfeasibletest = 0;
kInfeasibletest =  ((2*alpha)* max(expValueMeanAdjustedTask)) + (2*(1-alpha)*beta);
counter = 0;
disp("create MILPproblem")

MILPproblem = createKMILPProblemAverageExpressionIrrevTROR_Ajustment(taskModel, expValueMeanAdjustedTask, epsilon, rev2irrev, ktest, transportRxnIDs, alpha, beta);
%MILPproblem = load('MILPProblemTask1LoopLaw.mat').MILPproblem;


varargin = parseSolverParameters('MILP');
% varargin.feasTol = 10^-2;
% varargin.intTol = 10^-2;
% varargin.relMipGapTol = 10^-2;
% varargin.absMipGapTol = 10^-2;
% varargin.optTol = 10^-2;
% While loop
save('MILPProblemTask1LoopLaw.mat', 'MILPproblem');
x0 = [];
x0ON = true;
patternStr = "C";
lastBposition = min(strfind(MILPproblem.vartype(3*reactionLength:end)', patternStr)+3*reactionLength)-2;



ktest = kfeasibletest+(kInfeasibletest-kfeasibletest)*kPercentageConstant;
while abs(kfeasibletest - kInfeasibletest) > kTolerance
    disp({kfeasibletest;kInfeasibletest;ktest})
    %ktest = -100;
    % Adjust MILPproblem to change to new ktest
    %MILPproblem = createKMILPProblemAverageExpression(taskModel, expValueMeanAdjusted, epsilon, rev2irrev, ktest);
    tic
    
    MILPproblem = adaptKMILPProblemAverageExpressionTROR_Adjusted(MILPproblem, taskModel, expValueMeanAdjustedTask, ktest,transportRxnIDs, alpha, beta,x0,x0ON);
    % Calculate solution
    varargin.printLevel  = 3;
    solution = solveCobraMILP(MILPproblem,varargin);
    counter = counter +1;
    toc
    disp(solution.stat)
    % if solution was feasible, set kfeasible to ktest
    
    if solution.stat ==1
        

        %x0 = [solution.cont(1:reactionLength);solution.int(1:reactionLength*2);solution.int(1+reactionLength*2:lastBposition-reactionLength);solution.cont(reactionLength+1:end)]
        x0 = solution.full;
        expressionK = findTrueKValue(solution, MILPproblem, expValueMeanAdjustedTask, alpha, beta,reactionLength,ktest,SMatrix);
        %solution1 = load('lastValidSolutionMILP_TROR_AdjustTask1Sample4Scored6.383alpha0.200beta1.000.mat').solution;
        %[newSolution, prunedModel, originalModel] = runRxnsMinimalizationUsingIMAT(taskModel, solution, expValueMeanAdjustedTask,reactionLength, epsilon);



        if expressionK> ktest
            kfeasibletest = expressionK;
        else
            disp('expressionK is not higher than ktest, this should not be possible')
            kfeasibletest = ktest;
        end
        lastValidSolution = solution;

        %lastValidK=ktest;
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
    %disp(length(find((abs(solution.cont(1:reactionLength)),6)>= (epsilon-0.01))))
end 

% incase narrowing ends on invalid solution due to rounding (?)
if exist('lastValidSolution')
    solution = lastValidSolution;
    
    selectedRows = round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-0.01);
    validValues = expValueMeanAdjustedTask(selectedRows) ~= -1;
    selectedExpression = expValueMeanAdjustedTask(selectedRows);  

    score = sum(selectedExpression(validValues)) / length(selectedExpression(validValues)) ;
    
    reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-0.01);
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