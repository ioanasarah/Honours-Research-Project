i = 1;
j = 1;
epsilon = 1;

alpha = 0.5;
beta = 0;
osenseStr = "min";
irrev = true;
orphanHandling = "median";

cd("/home/jcornelis/Documents/MaCSBio-31.10/General/Functions/Metabolic_tasks/Jelle/Main Functions/Testing" + ...
    "T1S1a=0.5b=0")
MainFolder = pwd;
addpath(MainFolder)

% Load all modelse etc for testing of T1S1
load("expValueMeanAdjustedTaskIrrevKapprox.mat");
load("expValueMeanAdjustedTaskReversibleKApproximation.mat");
load("newIrrevModelKapprox.mat");
load("newReversibleModelForKApproxmiation.mat");
load("rev2irrevIrrevKapprox.mat");
load("CharnesCooperMILPproblem.mat");
load("expValueMeanAdjustedTask.mat")

fprintf("h")

%%
% Parameter to be tested
tested_values = [0, 1, 2, 3, 4, 5];
for test = tested_values

    currentDateTimeFolder = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    FolderNameRun = sprintf("RunNewCharnesCooper_Heur1_Method=%i_Date%s", test, currentDateTimeFolder);
    mkdir(FolderNameRun)
    copyfile('convertToEscherModelWithInputs.py', FolderNameRun)
    cd(FolderNameRun)

    diary logFile.txt

    varargin = parseSolverParameters('MILP');
    % Tested parameter
    varargin.Method = test;

    % Other parameters
    varargin.Heuristics = 0.1; %default 0.05
    varargin.printLevel  = 3;
    varargin.Presolve = 2;
    % varargin.Method = 4;
    varargin.Cuts = 3;
    varargin.NodeFileStart = 1;
    varargin.MIPFocus = 2;
    varargin.logFile = 'CharnesCooperLog.txt';
    varargin.NoRelHeurTime = 600;
    varargin.timeLimit = 3600;

    BoundsValue = 70;
    OrphanReactionValue = median(expValueMeanAdjustedTask(expValueMeanAdjustedTaskReversibleKApproximation>0)); % set to 0 if orphan reactions should not be included
    BoundsAdjust = true;
    
    %
    orphanReactionsIds = find(expValueMeanAdjustedTaskReversibleKApproximation ==-1);
    c_z = CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns));
    c_z(orphanReactionsIds) = OrphanReactionValue;
    CharnesCooperMILPproblem.c(4*length(newReversibleModelForKApproxmiation.rxns)+1:5*length(newReversibleModelForKApproxmiation.rxns)) = c_z;
    CharnesCooperMILPproblem.c(5*length(newReversibleModelForKApproxmiation.rxns)+1:6*length(newReversibleModelForKApproxmiation.rxns)) = c_z;  
    
    %
    CharnesCooperMILPproblemEdited = CharnesCooperMILPproblem;
    if BoundsAdjust
    tempUB = CharnesCooperMILPproblem.ub(1:length(newReversibleModelForKApproxmiation.rxns));
    tempLB = CharnesCooperMILPproblem.lb(1:length(newReversibleModelForKApproxmiation.rxns));
    tempUB(tempUB==1000) = BoundsValue;
    tempLB(tempLB== -1000) = -BoundsValue;
    
    CharnesCooperMILPproblemEdited.ub(1:length(newReversibleModelForKApproxmiation.rxns)) = tempUB;
    CharnesCooperMILPproblemEdited.lb(1:length(newReversibleModelForKApproxmiation.rxns)) = tempLB;            
    end
    %
    
    % Define parameter names and values
    parameterNames = {
    'Heuristics', 'printLevel', 'Presolve', 'Method', 'Cuts', 'NodeFileStart', 'MIPFocus', 'logFile', 'NoRelHeurTime', 'BoundsValue', 'BoundsAdjust','OrphanReactionValue', 'TimeLimit'
    };
    
    parameterValues = {
    varargin.Heuristics, varargin.printLevel, varargin.Presolve, varargin.Method, varargin.Cuts, varargin.NodeFileStart, varargin.MIPFocus, varargin.logFile, varargin.NoRelHeurTime,BoundsValue,BoundsAdjust,OrphanReactionValue, varargin.timeLimit
    };
    
    
    
    % Open a text file for writing (create or overwrite)
    fileID = fopen('parameter_values_CharnesCooper.txt', 'w');
    
    % Check if the file was successfully opened
    if fileID == -1
    error('Error opening the file for writing.');
    end
    
    % Loop through parameter names and values and write them to the file
    for i = 1:numel(parameterNames)
    fprintf(fileID, '%s = %s\n', parameterNames{i}, mat2str(parameterValues{i}));
    end
    
    % Close the file
    fclose(fileID);
    
    
    tic
    rng(1)
    CharnesCooperSolution = solveCobraMILP_Adjusted(CharnesCooperMILPproblemEdited, varargin);
    save('CharnesCooperSolution', 'CharnesCooperSolution')
    sumAverageExpressionCharnesCooper = sum(expValueMeanAdjustedTaskReversibleKApproximation(abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-0.1 &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
    scoreAverageCharnesCooper = sumAverageExpressionCharnesCooper/length(find( abs(CharnesCooperSolution.cont(1:length(newReversibleModelForKApproxmiation.rxns)))> epsilon-0.1 &expValueMeanAdjustedTaskReversibleKApproximation ~= -1));
    toc
    ActiveReactionFluxesCharnesCoopers = getActiveReactionsWithFluxesFromModelAndSaveModelCooper(newReversibleModelForKApproxmiation, CharnesCooperSolution,expValueMeanAdjustedTaskReversibleKApproximation,scoreAverageCharnesCooper , i, osenseStr , irrev,j);
    writecell(ActiveReactionFluxesCharnesCoopers,"ActiveReactionsTask_CharnesCooper.xlsx",'Sheet','ActReactionsTaskCharnes','Range','A1');
   
    
    command = sprintf('python convertToEscherModelWithInputs.py "%s"', "");  % Replace 'python' with 'python3' if needed
    [status, cmdout] = system(command);
    
    diary off
    nowCurrentFolder = pwd;
    if string(nowCurrentFolder) ~= string(MainFolder)
        cd(MainFolder)
    end
end

