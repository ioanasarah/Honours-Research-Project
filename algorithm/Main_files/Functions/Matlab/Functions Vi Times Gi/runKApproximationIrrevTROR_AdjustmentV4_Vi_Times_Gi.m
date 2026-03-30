function [solution,score,reactionsToKeep,LogData] = runKApproximationIrrevTROR_AdjustmentV4_Vi_Times_Gi(TaskModel,expValueMeanAdjustedTask,epsilon,rev2irrev,alpha,beta,iTask,jSample,osenseStr,bounds)

varargin = parseSolverParameters('MILP');
varargin2 = parseSolverParameters('MILP');

kTolerance = 0.01;
SetMedianForORExpressions = true;
TransportReactionPenalty = 0;
varargin.feasTol = 10^-5;
x0ON = true;

timeLimit_ConstantInitialInfeas = 240*5;
exceededInitialInfeasConstraint = false;
difference_ConstantInfeasibleHeuristic = 0.1;

timeLimit_ConstantFirstFeas = 600*3;
timeLimit_ConstantHeursticMultipleFeasilibity = 3;
approximatedOptimalSolutionThroughTimeHeuristic = false;
difference_ConstantfeasibleHeuristic = 0.1; % Should probably for low kTolerance, be qual to kTolerancekTolerance;
foundAtleastOneFeasibleSolution = false;
amountOfTryMoreTries = 1;
couldNotFindAnySolution = false;
wentOverTimeLimitAmount = 0;
wentOverTimeLimit = false;
foundSolutionAfterGoingOverLimit = false;
wentOverTimeLimitTwice = false;

varargin.MIPFocus = 1; %Should focus on finding any feasible value
varargin.NodeMethod = 1; % 1 should set to dual simplex can try 2 for barrier
varargin.Method = 5; % 1 should set to dual simplex, can try 2 for barrier
varargin.Presolve = 2; % 2 is agressive presolving
varargin.Cuts = 3; % we got good results without this, 3 is agressive cuts
varargin.timeLimit = timeLimit_ConstantFirstFeas;
varargin.SolutionLimit = 1;
varargin.printLevel  = 3;
varargin2.printLevel  = 3;

useSpecialYflags = false;
adjustBounds = false;
adjustBoundsValue = bounds;

% alpha = 0.5; % even addition of the TROR adjustment as regular obj funciton at 0.5
%              % alpha gives the weights to the TROR adjustment of the obj
%              % function and to the original objective function by setting
%
% beta = 1;  % total addition to max score if all oprhan transport reactions would be 0

reactionLength = length(TaskModel.rxns);
timerval = tic;
ktest = 0;
[~, ~, transportRxnIDs] = findTransRxns(TaskModel);
SMatrix = TaskModel.S;

expValueMeanAdjustedTaskForComparison = expValueMeanAdjustedTask;
orphanHandleParameter = 0;

if SetMedianForORExpressions
    orphanHandleParameter = median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0));
    expValueMeanAdjustedTask(expValueMeanAdjustedTask == -1) = orphanHandleParameter;
end

sampleMedian = median(expValueMeanAdjustedTask(expValueMeanAdjustedTask>0));

[kfeasibletest, optionsSaved] =  findHighestKIterativeLPs_Vi_Times_Gi(TaskModel,expValueMeanAdjustedTask,osenseStr);
fprintf('\nInitial LP Guess: Highest K = %f\n',kfeasibletest)
disp(optionsSaved)
%kfeasibletest = kfeasibletest+(0.5*beta);

kInfeasibletest =  ((2*alpha)* (1000 * max(expValueMeanAdjustedTask) + TransportReactionPenalty)) + (2*(1-alpha)*beta) + TransportReactionPenalty;

disp("Create MILP Problem")

if adjustBounds
    tempUB = TaskModel.ub(1:length(TaskModel.rxns));
    tempLB = TaskModel.lb(1:length(TaskModel.rxns));
    tempUB(tempUB== 1000) = adjustBoundsValue;
    tempLB(tempLB== -1000) = -adjustBoundsValue;
    %CharnesCooperMILPproblemEdited = CharnesCooperMILPproblem;
    TaskModel.ub(1:length(TaskModel.rxns)) = tempUB;
    TaskModel.lb(1:length(TaskModel.rxns)) = tempLB;
end

disp("_______________________________________________________________________")
disp("___________________________SHOWING OPTMIZER SETTINGS___________________")
disp("_______________________________________________________________________")

fprintf('UseSpecialYflags = %s\n', mat2str(useSpecialYflags));
fprintf('timeLimit_ConstantFirstFeas = %i, timeLimit_ConstantHeursticMultipleFeasilibity = %i\n', timeLimit_ConstantFirstFeas, timeLimit_ConstantHeursticMultipleFeasilibity);
fprintf('difference_ConstantfeasibleHeuristic = %.1f, amountOfTryMoreTries = %i\n', difference_ConstantfeasibleHeuristic, amountOfTryMoreTries);
fprintf('AdjustBounds = %s, AdjustBoundsValue = %i\n', mat2str(adjustBounds), adjustBoundsValue);

fprintf('kTolerance = %.2f\n', kTolerance);
fprintf('SetMedianForORExpressions = %d\n', SetMedianForORExpressions);
fprintf('TransportReactionPenalty = %d\n', TransportReactionPenalty);
fprintf('varargin.feasTol = %.4f\n', varargin.feasTol);
fprintf('x0ON = %d\n', x0ON);
fprintf('timeLimit_ConstantInitialInfeas = %d\n', timeLimit_ConstantInitialInfeas);
fprintf('difference_ConstantInfeasibleHeuristic = %.1f\n', difference_ConstantInfeasibleHeuristic);
fprintf('timeLimit_ConstantFirstFeas = %d\n', timeLimit_ConstantFirstFeas);
fprintf('timeLimit_ConstantHeursticMultipleFeasilibity = %d\n', timeLimit_ConstantHeursticMultipleFeasilibity);
fprintf('difference_ConstantfeasibleHeuristic = %.1f\n', difference_ConstantfeasibleHeuristic);
fprintf('amountOfTryMoreTries = %d\n', amountOfTryMoreTries);
fprintf('varargin.MIPFocus = %d\n', varargin.MIPFocus);
fprintf('varargin.NodeMethod = %d\n', varargin.NodeMethod);
fprintf('varargin.Method = %d\n', varargin.Method);
fprintf('varargin.Presolve = %d\n', varargin.Presolve);
fprintf('varargin.Cuts = %d\n', varargin.Cuts);
fprintf('varargin.timeLimit = %d\n', varargin.timeLimit);
fprintf('varargin.SolutionLimit = %d\n', varargin.SolutionLimit);
fprintf('varargin.printLevel = %d\n', varargin.printLevel);
fprintf('varargin2.printLevel = %d\n', varargin2.printLevel);

parameterNames = {
    'UseSpecialYflags',
    'amountOfTryMoreTries',
    'AdjustBounds',
    'AdjustBoundsValue',
    'kTolerance',
    'SetMedianForORExpressions',
    'Median',
    'TransportReactionPenalty',
    'varargin.feasTol',
    'x0ON',
    'timeLimit_ConstantInitialInfeas',
    'difference_ConstantInfeasibleHeuristic',
    'timeLimit_ConstantFirstFeas',
    'timeLimit_ConstantHeursticMultipleFeasilibity',
    'difference_ConstantfeasibleHeuristic',
    'amountOfTryMoreTries',
    'varargin.MIPFocus',
    'varargin.NodeMethod',
    'varargin.Method',
    'varargin.Presolve',
    'varargin.Cuts',
    'varargin.timeLimit',
    'varargin.SolutionLimit',
    'varargin.printLevel',
    'varargin2.printLevel',
    'Alpha',
    'Beta'
    };

parameterValues = {
    mat2str(useSpecialYflags),
    amountOfTryMoreTries,
    mat2str(adjustBounds),
    adjustBoundsValue,
    kTolerance,
    SetMedianForORExpressions,
    orphanHandleParameter,
    TransportReactionPenalty,
    varargin.feasTol,
    x0ON,
    timeLimit_ConstantInitialInfeas,
    difference_ConstantInfeasibleHeuristic,
    timeLimit_ConstantFirstFeas,
    timeLimit_ConstantHeursticMultipleFeasilibity,
    difference_ConstantfeasibleHeuristic,
    amountOfTryMoreTries,
    varargin.MIPFocus,
    varargin.NodeMethod,
    varargin.Method,
    varargin.Presolve,
    varargin.Cuts,
    varargin.timeLimit,
    varargin.SolutionLimit,
    varargin.printLevel,
    varargin2.printLevel,
    alpha,
    beta
    };

% Open a text file for writing (create or overwrite)
fileID = fopen('parameter_values.txt', 'w');

% Check if the file was successfully opened
if fileID == -1
    error('Error opening the file for writing.');
end

% Loop through parameter names and values and write them to the file
for iiterator = 1:numel(parameterNames)
    fprintf(fileID, '%s = %s\n', parameterNames{iiterator}, mat2str(parameterValues{iiterator}));
end

% Close the file
fclose(fileID);

MILPproblem = createKMILPProblemAverageExpressionIrrevTROR_Ajustment(TaskModel, expValueMeanAdjustedTask, epsilon, rev2irrev, ktest, transportRxnIDs, alpha, beta);

currentDateTime = string(datetime('now'), 'yyyy-MM-dd_HH-mm-SS');
filenameMILP = sprintf('LastCreatedMILPproblem_T%iS%i_date_%s.mat',iTask,jSample,currentDateTime);
save(filenameMILP, 'MILPproblem')

LogData = "";

if adjustBounds
    tempUB = MILPproblem.ub(1:length(TaskModel.rxns));
    tempLB = MILPproblem.lb(1:length(TaskModel.rxns));
    tempUB(tempUB== 1000) = adjustBoundsValue;
    tempLB(tempLB== -1000) = -adjustBoundsValue;
    MILPproblem.ub(1:length(TaskModel.rxns)) = tempUB;
    MILPproblem.lb(1:length(TaskModel.rxns)) = tempLB;
end


%% OPTIMIZATION
MILPproblem.c(1:reactionLength) = 1; % v  minimziation of v seems best
MILPproblem.c(reactionLength+1:2*reactionLength) = 0; % y+ %even though minimization of yi makes more sense
MILPproblem.osense = +1; %(-1 max, +1 min)
fprintf("\nOsense = %i, c_vector v = %i, c_vector y = %i\n",MILPproblem.osense,MILPproblem.c(1),MILPproblem.c(reactionLength+1) )

%%
x0hint = [];

%% New search strategy:
% Use LP to find initial K guess
% then use small time limit to find some initial infeasible solutions,
% if timelimit reached but no infeasible solution, increase size
% if timelimit not reached and infeasible is found, perform again
% then take k_feas and start finding the first next value, in small steps
% use time heuristic for determingin how long we take for the next one
% (something like 4x the last one (needs to be tested)).
%fprintf("timeLimit_ConstantInitialInfeas = %i, kTolerance = %.2f\n",timeLimit_ConstantInitialInfeas, kTolerance )


%% Time driven method of shrinking K Infeasible

ktest = kInfeasibletest;
kdiff = abs(kfeasibletest - kInfeasibletest);
infeasibleCounter = 0;

while ~exceededInitialInfeasConstraint && kdiff >= kTolerance
    infeasibleCounter = infeasibleCounter+1;
    varargin2.timeLimit = timeLimit_ConstantInitialInfeas;
    disp({kfeasibletest;kInfeasibletest;ktest})
    disp(" ")
    rng(1) % unfortunately only gurobi method 4 and 5 are deterministic
    MILPproblem = adaptKMILPProblem_Vi_Times_Gi(MILPproblem, TaskModel, expValueMeanAdjustedTask, ktest,transportRxnIDs, alpha, beta,x0hint,x0ON,TransportReactionPenalty);
    logFile = convertStringsToChars(sprintf("infeasible_space_run_%03i.txt",infeasibleCounter));
    varargin2.logFile = logFile;
    solution = solveCobraMILP_Adjusted(MILPproblem,varargin2);
    if solution.stat == 1|| solution.stat == 10 % Optimal (1) or feasible (10) found
        disp("found a feasible solution in the infeasible timed heuristic approach, this should not be possible")
        kfeasibletest = ktest;
        lastValidSolution = solution;
        getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(TaskModel, lastValidSolution,expValueMeanAdjustedTaskForComparison, kfeasibletest , iTask,jSample, alpha, beta, epsilon);
        if solution.stat == 1
            LogData = LogData + sprintf("Infeasible Space: Optimal solution found at k = %.3f \n", kfeasibletest);
        elseif solution.stat == 10
            LogData = LogData + sprintf("Infeasible Space: Feasible solution found at k = %.3f \n", kfeasibletest);
        end
    elseif solution.stat == 0
        kInfeasibletest = ktest;
    elseif solution.stat == 3 % Time limit reached
        exceededInitialInfeasConstraint = true;
    end
    kdiff = abs(kfeasibletest - kInfeasibletest);
    ktest = kInfeasibletest-(kdiff*difference_ConstantInfeasibleHeuristic);
end

%% Finding at least 1 feasible solution using Heuristics of Gurobi
% Using maximazation objective function correctly:
% in previous attempts I would set the c for Y+ flags to be positive but in
% reality this doesn't do anything with K-test as that is part of a
% constraint. Instead one should add the 2*alpha(expr-ktest) values
% directly in the c vector and try to maximize this value itself. This
% might actually increase performance, similarly including beta should help
% it move to better solutions quicker

disp(varargin)
if useSpecialYflags
    disp("_____________________________________________________________________________")
    disp("Adjust Yflags for optimziation")
    disp("_____________________________________________________________________________")
    MILPproblem.osense = -1; %(-1 max, +1 min)
    MILPproblem.c(reactionLength+1:2*reactionLength) = 0; %Y+
    MILPproblem.c(2*reactionLength+1:3*reactionLength) = 0; %Yprime
    MILPproblem.c(1:reactionLength) = 0; % v  minimziation of v seems best
    betaAdjusted = beta/length(find(expValueMeanAdjustedTask ==-1));

    for iterator = 1:reactionLength
        if(expValueMeanAdjustedTask(iterator) ~= -1)
            MILPproblem.c(reactionLength+iterator) = alpha*2*(expValueMeanAdjustedTask(iterator)-ktest); %*2 is there since 0.5 is alpha default value
        else
            MILPproblem.c(reactionLength+iterator) = 0; %*2 is there since 0.5 is alpha default value
        end
        if(expValueMeanAdjustedTask(iterator) == -1)
            MILPproblem.c(2*reactionLength+iterator) = ((1-alpha)*2*betaAdjusted); %*2 is there since 0.5 is alpha default value
        else
            MILPproblem.c(2*reactionLength+iterator)= 0;
        end
    end

end
%% ONLY FOR TESTING
%kInfeasibletest = 6.6203;

%% Feasibility testing
kdiff = abs(kfeasibletest - kInfeasibletest);
ktest = kfeasibletest + (kdiff*difference_ConstantfeasibleHeuristic);
feasibleCounter = 0;

while (~approximatedOptimalSolutionThroughTimeHeuristic && ~couldNotFindAnySolution && kdiff >= kTolerance)
    feasibleCounter = feasibleCounter + 1;
    disp({kfeasibletest;kInfeasibletest;ktest})
    disp(" ")
    MILPproblem = adaptKMILPProblem_Vi_Times_Gi(MILPproblem, TaskModel, expValueMeanAdjustedTask, ktest,transportRxnIDs, alpha, beta,x0hint,x0ON,TransportReactionPenalty);
    timerStart = tic;
    rng(1) % unfortunately only gurobi method 4 and 5 are deterministic
    logFile = convertStringsToChars(sprintf("feasible_space_run_%03i.txt",feasibleCounter));
    varargin.logFile = logFile;
    solution = solveCobraMILP_Adjusted(MILPproblem,varargin);

    if solution.stat == 1|| solution.stat == 10 % Optimal (1) or feasible (10) found
        kfeasibletest = findTrueKValue(solution, MILPproblem, expValueMeanAdjustedTask, alpha, beta,reactionLength,ktest,SMatrix,TransportReactionPenalty,transportRxnIDs);
        foundAtleastOneFeasibleSolution = true;
        SolvingTime =  ceil(toc(timerStart));

        if SolvingTime*timeLimit_ConstantHeursticMultipleFeasilibity >timeLimit_ConstantFirstFeas && SolvingTime*timeLimit_ConstantHeursticMultipleFeasilibity> varargin.timeLimit
            varargin.timeLimit = SolvingTime*timeLimit_ConstantHeursticMultipleFeasilibity;
            fprintf("FeasSolveTime was adjusted to: %i \n", timeLimit_ConstantHeursticMultipleFeasilibity)
        end

        if wentOverTimeLimit && ~wentOverTimeLimitTwice
            wentOverTimeLimitAmount = 0;
            foundSolutionAfterGoingOverLimit = true;
        end

        lastValidSolution = solution;
        getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(TaskModel, lastValidSolution,expValueMeanAdjustedTaskForComparison, kfeasibletest , iTask,jSample, alpha, beta, epsilon);

        x0hint = solution.full;

        if solution.stat == 1
            LogData = LogData + sprintf("Feasible Space: Optimal solution found at k = %.3f \n", kfeasibletest);
        elseif solution.stat == 10
            LogData = LogData + sprintf("Feasible Space: Feasible solution found at k = %.3f \n", kfeasibletest);
        end

    elseif solution.stat == 0
        kInfeasibletest = ktest;

    elseif solution.stat == 3  % Time limit reached
        wentOverTimeLimitAmount = wentOverTimeLimitAmount + 1;

        if ~foundAtleastOneFeasibleSolution % If we haven't found a solution then we need to reduce kfeas (bad initial LP guess)
            kfeasibletest = kfeasibletest-1;
            wentOverTimeLimitAmount = wentOverTimeLimitAmount - 1;
        elseif wentOverTimeLimitAmount <=1
            %diff = abs(kfeasibletest - kInfeasibletest);
            %kInfeasibletest = ktest+kTolerance*2;

            if foundSolutionAfterGoingOverLimit == true % if we already went over the limit, and found another optimal solution
                wentOverTimeLimitTwice = true; % then we change our approach
                % so that we don't end up running forever by not resetting wentOverTimeLimit
            end
            wentOverTimeLimit = true; % If this is the first time, and we find a subsequent feasible solution,
            % that will reset our wentOverTimeLimit counter
            difference_ConstantfeasibleHeuristic = difference_ConstantfeasibleHeuristic/2; %change our step size
        elseif wentOverTimeLimitAmount <= amountOfTryMoreTries
            difference_ConstantfeasibleHeuristic = difference_ConstantfeasibleHeuristic/2; % change our step size again
        elseif foundAtleastOneFeasibleSolution
            approximatedOptimalSolutionThroughTimeHeuristic = true; %if we try more than amountOfTryMoreTries after reset,
            % and we found a solution, we stop
        elseif wentOverTimeLimitAmount <= amountOfTryMoreTries+1
            couldNotFindAnySolution = true;
        end
    end

    kdiff = abs(kfeasibletest - kInfeasibletest);
    ktest = kfeasibletest + (kdiff*difference_ConstantfeasibleHeuristic);
end
kdiff = abs(kfeasibletest - kInfeasibletest);
ktest = kfeasibletest + (kdiff*difference_ConstantfeasibleHeuristic);

varargin2.SolutionLimit = 1;
varargin2.timeLimit = varargin.timeLimit;

%% Feasibility testing using normal settings 
% (found a solution using heuristics for task 20 where the above approach
% didn't so it seems worth trying

if approximatedOptimalSolutionThroughTimeHeuristic
    exceededInitialInfeasConstraint = false;
    while ~exceededInitialInfeasConstraint && kdiff >= kTolerance
        disp({kfeasibletest;kInfeasibletest;ktest})
        disp(" ")
        rng(1) % unfortunately only gurobi method 4 and 5 are deterministic
        MILPproblem = adaptKMILPProblem_Vi_Times_Gi(MILPproblem, TaskModel, expValueMeanAdjustedTask, ktest,transportRxnIDs, alpha, beta,x0hint,x0ON,TransportReactionPenalty);
        logFile = convertStringsToChars(sprintf("feasible_space_after_feasible_run_%03i.txt",feasibleCounter));
        varargin2.logFile = logFile;
        solution = solveCobraMILP_Adjusted(MILPproblem,varargin2);
        fprintf("Solstat = %i\n", solution.stat)
        if solution.stat == 1|| solution.stat == 10
            difference_ConstantfeasibleHeuristic = difference_ConstantfeasibleHeuristic*2;
            kfeasibletest = ktest;
            lastValidSolution = solution;
            getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(TaskModel, lastValidSolution,expValueMeanAdjustedTaskForComparison, kfeasibletest , iTask,jSample, alpha, beta, epsilon);

            if solution.stat == 1
                LogData = LogData + sprintf("Feasible Space after feasible space: Optimal solution found at k = %.3f \n", kfeasibletest);
            elseif solution.stat == 10
                LogData = LogData + sprintf("Feasible Space after feasible space: Feasible solution found at k = %.3f \n", kfeasibletest);
            end

        elseif solution.stat == 0
            kInfeasibletest = ktest;
        elseif solution.stat == 3 % Time limit reached
            exceededInitialInfeasConstraint = true;
        end
        kdiff = abs(kfeasibletest - kInfeasibletest);
        ktest = kfeasibletest + (kdiff*difference_ConstantfeasibleHeuristic);
    end
end

%% Another round of infeasibility testing to decrease the upper bound

kdiff = abs(kfeasibletest - kInfeasibletest);
ktest = kInfeasibletest - (kdiff * difference_ConstantInfeasibleHeuristic);

if exceededInitialInfeasConstraint
    varargin2.timeLimit = varargin.timeLimit;
    exceededInitialInfeasConstraint = false;
    while ~exceededInitialInfeasConstraint && kdiff >= kTolerance
        disp({kfeasibletest;kInfeasibletest;ktest})
        disp(" ")
        rng(1) % unfortunately only gurobi method 4 and 5 are deterministic
        MILPproblem = adaptKMILPProblem_Vi_Times_Gi(MILPproblem, TaskModel, expValueMeanAdjustedTask, ktest,transportRxnIDs, alpha, beta,x0hint,x0ON,TransportReactionPenalty);
        logFile = convertStringsToChars(sprintf("infeasible_space_after_feasible_run_%03i.txt",feasibleCounter));
        varargin2.logFile = logFile;
        solution = solveCobraMILP_Adjusted(MILPproblem,varargin2);
        fprintf("Solstat = %i", solution.stat)
        if solution.stat == 1|| solution.stat == 10
            kfeasibletest = findTrueKValue(solution, MILPproblem, expValueMeanAdjustedTask, alpha, beta,reactionLength,ktest,SMatrix,TransportReactionPenalty,transportRxnIDs);
            lastValidSolution = solution;
            getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(TaskModel, ...
                lastValidSolution,expValueMeanAdjustedTaskForComparison, ...
                kfeasibletest,iTask,jSample,alpha,beta,epsilon);

            if solution.stat == 1
                LogData = LogData + sprintf("Infeasible Space after feasible space: Optimal solution found at k = %.3f \n", kfeasibletest);
            elseif solution.stat == 10
                LogData = LogData + sprintf("Infeasible Space after feasible space: Feasible solution found at k = %.3f \n", kfeasibletest);
            end

        elseif solution.stat == 0
            kInfeasibletest = ktest;
        elseif solution.stat == 3 % Time limit reached
            exceededInitialInfeasConstraint = true;
        end
        kdiff = abs(kfeasibletest - kInfeasibletest);
        ktest = kInfeasibletest-(kdiff*difference_ConstantInfeasibleHeuristic);
    end
end


% incase narrowing ends on invalid solution due to rounding (?)
if exist('lastValidSolution')
    solution = lastValidSolution;

    %selectedRows = round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-0.01);
    %validValues = expValueMeanAdjustedTask(selectedRows) ~= -1;
    %selectedExpression = expValueMeanAdjustedTask(selectedRows);

    %score = sum(selectedExpression(validValues)) / length(selectedExpression(validValues)) ;
    score = findTrueKValue(solution, MILPproblem, expValueMeanAdjustedTask, alpha, beta,reactionLength,ktest,SMatrix,TransportReactionPenalty,transportRxnIDs);
    reactionsToKeep = round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-(epsilon/100));
else
    solution = 0;
    score = 0;
    reactionsToKeep = 0;
end

elapsedTime = toc(timerval);
fprintf("\nK Approximation was done in %f seconds.\n", elapsedTime);

%% Recording variables interesting for run comparisons
LogData = LogData + sprintf("\nNumber of reactions in FastCC Irreversible Model: %d \nMedian Expression= %f \n\nTask Score is %f \nK infeasible= %f, K feasible = %f and last Ktest= %f \nGap = %f \nNumber of infeasible runs= %d \nNumber of feasible runs= %d \n", length(TaskModel.rxns), orphanHandleParameter, score, kInfeasibletest, kfeasibletest, ktest, kdiff, infeasibleCounter, feasibleCounter);
if kdiff > kTolerance
    f = fopen('__NOTconverged', 'wt');
    fclose(f);
else
    f = fopen('__converged', 'wt');
    fclose(f);
end

end
