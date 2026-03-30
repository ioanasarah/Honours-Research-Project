%%%%
% Clear vars for next run
clearvars

% Necessary evewry time
%initCobraToolbox
%changeCobraSolver('gurobi', 'all');

% Load some model
% export to Escher FBA
% # Name, flux 
% PGI, 5.445
% PFK, 5.303
% Scan epsilon values 1-4 to 1e-6
% And feasTol values 1e-6 to 1e-9
% calculate manual vs obj/average


Oldmodel = readCbModel('e_coli_core.mat');
OldModelSaved = Oldmodel;


%% create super simple toymodel
%1-2-3-4-5
%1(2)-6(1)-7(2)-8(nan)-5(3)
%set in orphan reaction to test/toy model

turnOffReactionsAndTurnOnOxphos = true;
setOxphosToyModel = true;

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


%% Test expression data
% Set oxidative phosphorylation reactions to very high for checking
oxphosRxnIds = findRxnIDs(Oldmodel,{'EX_glc__D_e','EX_o2_e','EX_h2o_e','EX_co2_e'});
expressionRxns(oxphosRxnIds) = 0;

% Test reactions for alternative flux
testRxnLowIds = findRxnIDs(Oldmodel,{'EX_acald_e', 'ACALD','ACALDt'});
% expressionRxns(testRxnLowIds) = 0.000000001;

testRxnHighIds = findRxnIDs(Oldmodel,{'ACONTa'});
% expressionRxns(testRxnHighIds) = 1000;

%%





% Indice the reactions of the model
edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end

% Vars
epsilon = 0.001;

% Possibly adjuts tolerances
varargin = parseSolverParameters('MILP');
adjustTol = true;
if adjustTol 
    varargin.intTol=1.0E-12; % default
    varargin.relMipGapTol=1.0E-12; % default
    varargin.absMipGapTol=1.0E-12; % default
    varargin.feasTol=1.0E-9; % This tolerance determines when a solution is
    % treated as feasible. In MILP, constraints and variable bounds may not
    % be satisfied exactly due to numerical errors. The feasTol parameter 
    % specifies the maximum violation allowed for a constraint or variable
    % bound to be considered as satisfied. If the violation of a constraint's
    % activity or a variable's bound is less in absolute magnitude than feasTol,
    % it is treated as satisfied. A smaller feasTol value indicates a stricter
    % tolerance for feasibility.
    varargin.optTol=1.0E-6; % default
end


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

% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1yprime*rxn + zprime*rxn+ 1*t
A_ColumnLength = size(S,2) + YAmount*reactionLength + ZAmount*reactionLength + 1*reactionLength  + TAmount;

% Row lenght equals S matrix + UB/LB*rxn +1rxn for yprime +1*rxn (for y+ + y- <= 1) + 3z*rxn + +   2*t + 1 (for z
% constraint) 
A_RowLength = size(S,1) + 2*reactionLength + 1*reactionLength + 1*reactionLength +4*2*reactionLength + 2*reactionLength+ TAmount*2 + 1;

%Create A0
A = sparse(A_RowLength, A_ColumnLength);

%Fill A with S
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

%% Order
%Y+ LB -epsilon
%Y- UB; + epsilon
% Y+ + Y- <=0
% Y+ + Y- + Yprime = 1
%Z>=0  'G' 0 
%Z<= M*Y == Z-(M*Y) <= 0 'L' 0
% M*y - z + t <= M 'L' M 
% M*y + z - t <= M 'L' M
%Z>=0  'G' 0 
%Z<= M*Y == Z-(M*Y) <= 0 'L' 0
% M*y - z + t <= M 'L' M 
% M*y + z - t <= M 'L' M
% Yprime LB
% Yprime UB
% Z = 1 single row (twice as long for z+ and z-)
% 2 single entries for 0>=t>=M
%%
orphanReactions = {1,2,3,4};
%Fill A matrix with values
for i = 1:length(edgeIndices)

    %Y+ LB -epsilon
    A(i+size(S,1),edgeIndices(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(edgeIndices(i)) - epsilon;

    %Y- UB; + epsilon
    A(i+size(S,1)+length(edgeIndices),edgeIndices(i)) = 1;
    A(i+size(S,1)+length(edgeIndices),reactionLength+i+size(S,2)) = ub(edgeIndices(i)) + epsilon;

    %Y+ + Y- <= 0
    A(i+size(S,1)+2*length(edgeIndices),i + 1*size(S,2)) = 1;
    A(i+size(S,1)+2*length(edgeIndices),i + 2*size(S,2)) = 1;

    %Y+ + Y- + Yprime = 1
    A(i+size(S,1)+3*length(edgeIndices),i + 1*size(S,2)) = 1;
    A(i+size(S,1)+3*length(edgeIndices),i + 2*size(S,2)) = 1;
    A(i+size(S,1)+3*length(edgeIndices),i + 3*size(S,2)) = 1;


    %% Z Plus blocks
    %Z>=0  'G' 0 
    A(i+size(S,1)+4*length(edgeIndices),i + 4*size(S,2)) = 1;

    %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
    A(i+size(S,1)+5*length(edgeIndices),i + 1*size(S,2)) = -M;
    A(i+size(S,1)+5*length(edgeIndices),i + 4*size(S,2)) = 1;

    % M*y + z-t <= M 'L' M
    A(i+size(S,1)+6*length(edgeIndices),i + 1*size(S,2)) = M;
    A(i+size(S,1)+6*length(edgeIndices),i + 4*size(S,2)) = 1;
    A(i+size(S,1)+6*length(edgeIndices),end) = -1;

    % M*y - z + t <= M 'L' M 
    A(i+size(S,1)+7*length(edgeIndices),i + 1*size(S,2)) = M;
    A(i+size(S,1)+7*length(edgeIndices),i + 4*size(S,2)) = -1;
    A(i+size(S,1)+7*length(edgeIndices),end) = 1;

    %% Z minus blocks
    %Z>=0  'G' 0 
    A(i+size(S,1)+8*length(edgeIndices),i + 4*size(S,2)) = 1;

    %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
    A(i+size(S,1)+9*length(edgeIndices),i + 2*size(S,2)) = -M;
    A(i+size(S,1)+9*length(edgeIndices),i + 5*size(S,2)) = 1;

    % M*y + z-t <= M 'L' M
    A(i+size(S,1)+10*length(edgeIndices),i + 2*size(S,2)) = M;
    A(i+size(S,1)+10*length(edgeIndices),i + 5*size(S,2)) = 1;
    A(i+size(S,1)+10*length(edgeIndices),end) = -1;

    % M*y - z + t <= M 'G' -M 
    A(i+size(S,1)+11*length(edgeIndices),i + 2*size(S,2)) = M;
    A(i+size(S,1)+11*length(edgeIndices),i + 5*size(S,2)) = -1;
    A(i+size(S,1)+11*length(edgeIndices),end) = 1;


    %% Other blocks
    % Yprime LB
    A(i+size(S,1)+12*length(edgeIndices),edgeIndices(i)) = 1;
    A(i+size(S,1)+12*length(edgeIndices),i+ 3*size(S,2)) = lb(edgeIndices(i));

    % Yprime UB
    A(i+size(S,1)+13*length(edgeIndices),edgeIndices(i)) = 1;
    A(i+size(S,1)+13*length(edgeIndices),i+ 3*size(S,2)) = ub(edgeIndices(i));
    
    % Z = 1 single row (twice as long for z+ and z-)
    %%% Set orphan to 0 in addition to c_z (obj.coef) of orphan to 0
    A(size(S,1)+14*length(edgeIndices) + 1, i + 4*size(S,2)) = 1;
    A(size(S,1)+14*length(edgeIndices) + 1, i + 5*size(S,2)) = 1;

    %% 2 single entries for 0>=t>=M
    % Question: Should this not be included in lb/ub variable of T?
    A(size(S,1)+14*length(edgeIndices) + 2, end) = 1;
    A(size(S,1)+14*length(edgeIndices) + 3, end) = 1;    
    
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

lb = [lb;lb_y;lb_y;lb_y;lb_z;lb_z;lb_t];
ub = [ub;ub_y;ub_y;ub_y;ub_z;ub_z;ub_t];

% Creating vartype (equal to column length:
% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t )
vartype_v(1:size(S,2),1) = 'C'; % Vi
vartype_y(1:3*reactionLength,1) = 'B'; %Yi + - 
vartype_z(1:2*reactionLength,1) = 'C'; %Zi + -
vartype_t(1,1) = 'C'; % t

vartype = [vartype_v;vartype_y;vartype_z;vartype_t];

% Creating c (equal to column length:
% Column length equals S matrix + 2Y*rxn + 2z*rxn + 1*t )
c_v = zeros(reactionLength,1); % we don't maximize v
%%c_y = ones(reactionLength,1); % don't maximize y
c_y = zeros(reactionLength,1); % don't maximize y
c_yprime = zeros(reactionLength,1); % don't maximize yprime
c_z = expressionRxns; % maximize z+ and z- equal to expression
c_zprime = zeros(reactionLength,1); % don't maximize zprime
c_t = zeros(1,1); % we don't maximize t

c = [c_v;c_y;c_y;c_yprime;c_z;c_z;c_t];

% Creating B
b_s = zeros(size(S,1),1); % SV = 0 equality
lb_indices = lb(edgeIndices); % >= Y_LB
ub_indices = ub(edgeIndices); % <= Y_UB
b_YPlusY = ones(size(S,2),1);  % <= Y+ + Y- <=1
b_Yprime = ones(size(S,2),1); % y+ + y- + yprime = 1 
b_z1 = zeros(size(S,2),1); %Z>=0  'G' 0 
b_z2 = zeros(size(S,2),1); %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
b_z3 = M*ones(size(S,2),1); % M*y - z + t <= M 'L' M 
b_z4 = M*ones(size(S,2),1); % M*y + z - t <= M 'L' M

b_yprime_LB = lb(edgeIndices);  % Yprime >= LB
b_yprime_UB = ub(edgeIndices); % Yprime  >= UB
b_z = 1; 
b_t_M_LB = 0; 
b_t_M_UB = M;
b = [b_s;lb_indices;ub_indices;b_YPlusY;b_Yprime;b_z1; b_z2 ;b_z3; b_z4; b_z1; b_z2; b_z3; b_z4;  b_yprime_LB;b_yprime_UB;   b_z;b_t_M_LB;b_t_M_UB];
%

%% Order
%SV= 0
%Y+ LB -epsilon
%Y- UB; + epsilon
% Y+ + Y- <=0
% Y+ + Y- + Yprime = 1
%Z>=0  'G' 0 
%Z<= M*Y == Z-(M*Y) <= 0 'L' 0
% M*y - z + t <= M 'L' M 
% M*y + z - t <= M 'L' M
%Z>=0  'G' 0 
%Z<= M*Y == Z-(M*Y) <= 0 'L' 0
% M*y - z + t <= M 'L' M 
% M*y + z - t <= M 'L' M
% Yprime LB
% Yprime UB
% Z = 1 single row (twice as long for z+ and z-)
% 2 single entries for 0>=t>=M
%%

% Creating csense
csense1(1:size(S,1)) = 'E'; % SV = 0 equality
csense2(1:length(edgeIndices)) = 'G'; % Y+ >= Y_LB
csense3(1:length(edgeIndices)) = 'L'; % Y- <= Y_UB
csense4(1:length(edgeIndices)) = 'L'; % Y+ + Y- <=1
csense5(1:length(edgeIndices)) = 'E'; % Y+ + Y- + Yprime = 1
csense6(1:length(edgeIndices)) = 'G'; %Z>=0  'G' 0 
csense7(1:length(edgeIndices)) = 'L'; %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
csense8(1:length(edgeIndices)) = 'L'; % M*y - z + t <= M 'L' M 
csense9(1:length(edgeIndices)) = 'L'; % M*y + z - t <= M 'L' M
csense10(1:length(edgeIndices)) = 'G'; %Yprime >lb
csense11(1:length(edgeIndices)) = 'L'; %Yprime >lb
csense12(1) = 'E'; % z = 1
csense13(1) = 'G'; % t >= 0
csense14(1) = 'L'; % t <= M
csense = [csense1 csense2 csense3 csense4 csense5 csense6 csense7 csense8 csense9 csense6 csense7 csense8 csense9 csense10   csense11 csense12 csense13 csense14];

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


TurnOnLoopConstraints = true; 

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



% Perform test
solution = solveCobraMILP(MILPproblem,varargin);
a = getCobraSolverParamsOptionsForType('MILP');
disp(solution.stat)

if TurnOnLoopConstraints
    t_solution = solution.cont(end-2:end); %with looplaw constraints, which one is t? > third last
else
    t_solution = solution.cont(end);
end
disp(t_solution)

% Get int solutions for checking
ReactionSolutions = zeros(95, 8);
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

activeRxnsNum =  nnz([ReactionSolutions{:, 3}]);
AverageGeneExpression = solution.obj/activeRxnsNum;

columnNames = {'Yi_F/R', 'Zi_F/R', 'Flux', 'Expression', 'rxnName'};
ReactionSolutions = cell2table(ReactionSolutions, 'VariableNames', columnNames);

% tol = epsilon-(epsilon/10);
activeReactions = ReactionSolutions(abs([ReactionSolutions{:, 3}]) >= epsilon, :);

solution.obj