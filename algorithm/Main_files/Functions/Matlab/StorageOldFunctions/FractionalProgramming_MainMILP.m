%%%%
% Clear vars for next run
clearvars

% Necessary evewry time
%initCobraToolbox
%changeCobraSolver('gurobi', 'all');

% Load some model
Oldmodel = readCbModel('e_coli_core.mat');
OldModelSaved = Oldmodel;



turnOffReactionsAndTurnOnOxphos = false;
setOxphosToyModel = false;

%%% Exchange Reactions
% Check the format of the model
if turnOffReactionsAndTurnOnOxphos
    if size(Oldmodel.rxns,2)>size(Oldmodel.rxns,1)
        Oldmodel.rxns=Oldmodel.rxns';
    end
    
    %Find all exchange/demand/sink reactions
    Exchange = {};
    
    for k=1:length(Oldmodel.rxns)
    
        if strlength(string(Oldmodel.rxns(k)))>2
            if  extractBetween(string(Oldmodel.rxns(k)),1,3) == 'EX_'
                Exchange(end+1) = Oldmodel.rxns(k);
            end
        end
        if sum(abs(Oldmodel.S(:,k))) == 1
            Exchange(end+1) = Oldmodel.rxns(k);
        end
    end
    Exchange=unique(Exchange);
    ExchangeInd = findRxnIDs(Oldmodel,Exchange);
    nonExchangesInd = setdiff(1:length(Oldmodel.lb),ExchangeInd) ;
    Oldmodel.lb(16) = 0;
    
    %Checking for lower bounds above 0 in exchanges.
    if sum(Oldmodel.lb(ExchangeInd) > 0) > 0
        msg = 'At least one exchange reaction has a lower bound higher than 0. This will make tasks with "ALLMETS" fail';
        warning(msg);
    end
    
    %Checking for lower bounds above 0 in non-exchanges.
    if sum(Oldmodel.lb(nonExchangesInd) > 0) > 0
        msg = 'At least one non-exchange reaction has a lower bound higher than 0. This could make ANY task fail';
        warning(msg);
    end
    
    
    %Close all exchange reactions
    previousLB = Oldmodel.lb;
    previousUB = Oldmodel.ub;
    Oldmodel.lb(findRxnIDs(Oldmodel,Exchange))=0;
    Oldmodel.ub(findRxnIDs(Oldmodel,Exchange))=0;
    
    
    if setOxphosToyModel 
        %%%% Toy example oxidative phosphorylation
        Oldmodel = changeRxnBounds(Oldmodel,'EX_glc__D_e',-10);
        Oldmodel = changeRxnBounds(Oldmodel,'EX_o2_e',-60); % set to -38
        Oldmodel = changeRxnBounds(Oldmodel,'EX_h2o_e',60);
        Oldmodel = changeRxnBounds(Oldmodel,'EX_co2_e',60);

    end
    % Add in ATP/ADP constraints and phosphate, change water
end 

model = Oldmodel;

% construct S
S = model.S;
S=full(S);
spy(S)

%Set Seed
rng('default')

% Add in random expression data
expressionRxns = 1*abs(randn(length(model.rxns),1)); %added +2 due to how loopLawConstraintsWorks (adds values, unsure what they actually relate to)

% Set oxidative phosphorylation reactions to very high for checking
oxphosRxnIds = findRxnIDs(Oldmodel,{'EX_glc__D_e','EX_o2_e','EX_h2o_e','EX_co2_e'});
expressionRxns(oxphosRxnIds) = 1;

% Test reactions for alternative flux
testRxnIds = findRxnIDs(Oldmodel,{'FRUpts2'});
expressionRxns(testRxnIds) = 1000000;
% expressionRxns = ones(length(model.rxns),1);

% Indice the reactions of the model
edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end

% Vars
epsilon = 0.01;

% Get S
S = model.S;
lb = model.lb; 
ub = model.ub;

% Create A size vars
YAmount = 2;
ZAmount = YAmount;
TAmount = 1;
reactionLength = length(edgeIndices);
M = 1;

% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t
A_ColumnLength = size(S,2) + YAmount*reactionLength + ZAmount*reactionLength + TAmount;

% Row lenght equals S matrix + UB/LB*rxn +1*rxn (for y+ + y- <= 1) + z*rxn + 2*t + 1 (for z
% constraint) 
A_RowLength = size(S,1) + 2*reactionLength + 1*reactionLength + ZAmount*reactionLength + TAmount*2 + 1;

%Create A
A = sparse(A_RowLength, A_ColumnLength);

%Fill A with S
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end


%% Order
%Y+ LB -epsilon
%Y- UB; + epsilon
%Y+ + Y- <=0
%Z+ link Y+
% Z- link Y-
% Z = 1 single row (twice as long for z+ and z-)
% 2 single entries for 0>=t>=M
%%


%Fill A matrix with values
for i = 1:length(edgeIndices)

    %Y+ LB -epsilon
    A(i+size(S,1),edgeIndices(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(edgeIndices(i))-epsilon;

    %Y- UB; + epsilon
    A(i+size(S,1)+length(edgeIndices),edgeIndices(i)) = 1;
    A(i+size(S,1)+length(edgeIndices),reactionLength+i+size(S,2)) = ub(edgeIndices(i))+epsilon;

    %Y+ + Y- <=0
    A(i+size(S,1)+2*length(edgeIndices),i + 1*size(S,2)) = 1;
    A(i+size(S,1)+2*length(edgeIndices),i+2*size(S,2)) = 1;

    %Z+ link Y+
    A(i+size(S,1)+3*length(edgeIndices),i + size(S,2)) = M ;
    A(i+size(S,1)+3*length(edgeIndices),i+3*size(S,2)) = 1;
    A(i+size(S,1)+3*length(edgeIndices),end) = -1;

    % Z- link Y-
    A(i+size(S,1)+4*length(edgeIndices),i + 2*size(S,2)) = M ;
    A(i+size(S,1)+4*length(edgeIndices),i+4*size(S,2)) = 1; 
    A(i+size(S,1)+4*length(edgeIndices),end) = -1;
   
    % Z = 1 single row (twice as long for z+ and z-)
    A(size(S,1)+5*length(edgeIndices) + 1, i + 3*size(S,2)) = 1;
    A(size(S,1)+5*length(edgeIndices) + 1, i + 4*size(S,2)) = 1;
    
    %% 2 single entries for 0>=t>=M
    % Question: Should this not be included in lb/ub variable of T?
    A(size(S,1)+5*length(edgeIndices) + 2, end) = 1;
    A(size(S,1)+5*length(edgeIndices) + 3, end) = 1;
    
end

spy(A)

% Creating lb and ub (equal to column length:
% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t )
lb_y = zeros(length(edgeIndices),1);
ub_y = ones(length(edgeIndices),1);
lb_z = zeros(length(edgeIndices),1);
ub_z = ones(length(edgeIndices),1);
lb_t = 0;
ub_t = M;

lb = [lb;lb_y;lb_y;lb_z;lb_z;lb_t];
ub = [ub;ub_y;ub_y;ub_z;ub_z;ub_t];

% Creating vartype (equal to column length:
% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t )
vartype_v(1:size(S,2),1) = 'C'; % Vi
vartype_y(1:YAmount*reactionLength,1) = 'B'; %Yi + - 
vartype_z(1:YAmount*reactionLength,1) = 'B'; %ZI + -
vartype_t(1,1) = 'C'; % t

vartype = [vartype_v;vartype_y;vartype_z;vartype_t];

% Creating c (equal to column length:
% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t )
c_v = zeros(reactionLength,1); % we don't maximize v
%%c_y = ones(reactionLength,1); % maximize, Should equal to expression data
c_y = expressionRxns; % maximize, equals to expression data
c_z = zeros(reactionLength,1); % we don't maximize z
c_t = zeros(1,1); % we don't maximize t

c = [c_v;c_y;c_y;c_z;c_z;c_t];

% Creating B
b_s = zeros(size(S,1),1); % SV = 0 equality
lb_indices = lb(edgeIndices); % >= Y_LB
ub_indices = ub(edgeIndices); % <= Y_UB
b_YPlusY = ones(size(S,2),1);  % <= Y+ + Y- <=1
b_ZlinkY1 = M*ones(size(S,2),1); % z_link_y <= M
b_ZlinkY2 = M*ones(size(S,2),1); % z_link_y <= M (Same row as above)
b_z = 1;
b_t_M_LB = 0;
b_t_M_UB = M;
b = [b_s;lb_indices;ub_indices;b_YPlusY;b_ZlinkY1;b_ZlinkY2;b_z;b_t_M_LB;b_t_M_UB];
%b_ZlinkY2;b_z;b_t_M_LB;b_t_M_UB

%% Order
%SV = 0 
%Y+ LB -epsilon
%Y- UB; + epsilon
%Y+ + Y- <=0
%Z+ link Y+
%Z- link Y-
% Z = 1 single row (twice as long for z+ and z-)
% t >= 0 single row
% t <= M single row
%%

% Creating csense
csense1(1:size(S,1)) = 'E'; % SV = 0 equality
csense2(1:length(edgeIndices)) = 'G'; % >= Y_LB
csense3(1:length(edgeIndices)) = 'L'; % <= Y_UB
csense4(1:length(edgeIndices)) = 'L'; % <= Y+ + Y- <=1
csense5(1:length(edgeIndices)) = 'L'; % z_link_y <= M  
csense6(1:length(edgeIndices)) = 'L'; % z_link_y <= M (Same row as above)
csense7(1) = 'E'; % z = 1
csense8(1) = 'G'; % t >= 0
csense9(1) = 'L'; % t <= M
csense = [csense1 csense2 csense3 csense4 csense5 csense6 csense7 csense8  csense9];

% create MILPproblem
MILPproblem.A = A; % A Matrix
MILPproblem.b = b; % Constraint vector
MILPproblem.c = c; % Objective coeficients (v and Y)
MILPproblem.lb = lb; % Lower bounds 
MILPproblem.ub = ub; % Upper bounds
MILPproblem.csense = csense; % Equality/lower/higher 
MILPproblem.vartype = vartype; % Vartype ('C', 'B', 'I')
MILPproblem.osense = -1; %(-1 max, +1 min)
MILPproblem.x0 = []; % Initial guess value





TurnOnLoopConstraints = false; 


% add LoopConstraints
if TurnOnLoopConstraints
    MILPproblem = addLoopLawConstraints(MILPproblem,model);
    % Adjust csense since loop constraints does something weird
    tempcsense = MILPproblem.csense;
    tempcsense2 = tempcsense(1,:);
    for i = 2:size(tempcsense,1)
        tempcsense2(end+1) = tempcsense(i);
    end
    MILPproblem.csense = tempcsense2;
end


% peforming test
solution = solveCobraMILP(MILPproblem);
disp(solution.stat)

if TurnOnLoopConstraints
    t_solution = solution.cont(end-2:end); %with looplaw constraints, which one is t? > third last
else
    t_solution = solution.cont(end);
end
disp(t_solution)

% Get int solutions for checking
ReactionSolutions = zeros(95, 4);
ReactionSolutions(1:reactionLength,1) = solution.int(1:reactionLength); %forward
ReactionSolutions(1:reactionLength,2) = solution.int(reactionLength+1:reactionLength+reactionLength); %backward
ReactionSolutions(1:reactionLength,4) = solution.cont(1:reactionLength);

% Check if forward/backward reactions are both active (should not)
% set 1 for forward, -1 for backward
for i = 1:size(ReactionSolutions,1)
    if ReactionSolutions(i,1) == 1 && ReactionSolutions(i,2) ==1 
        fprintf("rxn %s has both forward and backward reactions active", i)
    end
     
    if ReactionSolutions(i,1) == 1
        ReactionSolutions(i,3) = 1;
    elseif ReactionSolutions(i,2) == 1
        ReactionSolutions(i,3) = -1;
    end
end

ReactionSolutions = num2cell(ReactionSolutions);
ReactionSolutions = [ReactionSolutions, model.rxns];

activeRxns =  nnz([ReactionSolutions{:, 3}]);
AverageGeneExpression = solution.obj/activeRxns;


% 
% 
% 
% 
% 
% 
% 
% rxnIndex = 1:size(model.S,2);
% vartype = 'MI*P';
% MILPproblem = addLoopLawConstraints(LPProblem,model, rxnIndex,1); % Method is set to 1, for forward & Reversible
% 
% A = full(MILPproblem.A);
% 
% 
% %%%% Binary Search
% % Start by checking current LargestK, if infeasible break
% ktest = 0;
% MILPproblem = FractionalProgramming_UpdateMILPProblem(MILPproblem,expressionRxns,ktest);
% 
% %adjust cSense
% tempcsense = MILPproblem.csense;
% tempcsense2 = tempcsense(1,:);
% for i = 2:size(tempcsense,1)
%     tempcsense2(end+1) = tempcsense(i);
% end
% MILPproblem.csense = tempcsense2;
% 
% 
% kTolerance = 0.0001;
% kfeasibletest = 0;
% kInfeasibletest = 9000; %max(expressionRxns);
% counter = 0;
% 
% % While loop
% while abs(kfeasibletest - kInfeasibletest) > kTolerance
%     % Calculate new ktest
%     ktest = (kfeasibletest+kInfeasibletest)/2;
%     disp({kfeasibletest;kInfeasibletest;ktest})
% 
%     % Adjust MILPproblem to change to new ktest
%     MILPproblem = FractionalProgramming_UpdateMILPProblem(MILPproblem,expressionRxns,ktest);
%     % spy(MILPproblem.A)
%     % Calculate solution
%     solution = solveCobraMILP(MILPproblem);
%     counter = counter +1;
% 
%     disp(solution.stat)
%     % if solution was feasible, set kfeasible to ktest
%     if solution.stat ==1
%         kfeasibletest = ktest;
%         lastValidSolution = solution;
%         lastValidK=ktest;
%     % if solution was infeasible, set kinfeasible to ktest
%     elseif solution.stat == -1 
%         kInfeasibletest = ktest;
%     else
%         fprintf("Solution stat was not 1 or -1, namely: %s", "solution.stat")
%     end
%     disp(length(find(solution.cont)))
% end 
% 
% % incase narrowing ends on invalid solution due to rounding (?)
% solution = lastValidSolution;
% 
% % remove all reactions not in solution
% x = solution.cont; 
% %map to reversible model
% 
% %Convert IrrevModelSolution to Reversible Model
% ContSolutionIrrev = solution.cont(1:end-1);
% ContSolutionReversible = solution.cont(1:length(Oldmodel.rxns));
% for i =length(Oldmodel.rxns)+1:length(modelIrrev.rxns)-1 %-1 to remove pseudoreaction 'net flux'
%     if abs(ContSolutionIrrev(i))>0.00001 %add tol var
%         fluxValue = -1*ContSolutionIrrev(i); %invert value for reversed reactions
%         name = extractBetween(string(modelIrrev.rxns(i)),1,strlength(string(modelIrrev.rxns(i)))-1);
%         for j = 1:length(Oldmodel.rxns)
%             if name ==  extractBetween(string(modelIrrev.rxns(j)),1,strlength(string(modelIrrev.rxns(j)))-1);
%                 if ContSolutionReversible(j) ==0
%                 ContSolutionReversible(j) = fluxValue ;
%                 else
%                 fprintf("Reaction %s has both forward and backward non-zero value",name)
%                 end
%             end
%         end
%     end
% end
% 
% 
% %backwards is minus
% %%convert to cvs
% rxnRemList = modelIrrev.rxns(abs(x) < 0.0001); %add in variable 
% tissueModel = removeRxns(modelIrrev,rxnRemList);
% 
% % run optmizeCBModel to check if still valid model
% CBSolution = optimizeCbModel(tissueModel, 'min');
