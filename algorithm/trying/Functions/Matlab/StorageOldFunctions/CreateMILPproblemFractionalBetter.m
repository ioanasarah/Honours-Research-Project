
model = readCbModel('e_coli_core.mat');
if size(model.rxns,2)>size(model.rxns,1)
    model.rxns=model.rxns';
end

%Find all exchange/demand/sink reactions
Exchange = {};

for k=1:length(model.rxns)

    if strlength(string(model.rxns(k)))>2
        if  extractBetween(string(model.rxns(k)),1,3) == 'EX_'
            Exchange(end+1) = model.rxns(k);
        end
    end
    if sum(abs(model.S(:,k))) == 1
        Exchange(end+1) = model.rxns(k);
    end
end

Exchange=unique(Exchange);
ExchangeInd = findRxnIDs(model,Exchange);
nonExchangesInd = setdiff(1:length(model.lb),ExchangeInd) ;
model.lb(16) = 0;

%Checking for lower bounds above 0 in exchanges.
if sum(model.lb(ExchangeInd) > 0) > 0
    msg = 'At least one exchange reaction has a lower bound higher than 0. This will make tasks with "ALLMETS" fail';
    warning(msg);
end

%Checking for lower bounds above 0 in non-exchanges.
if sum(model.lb(nonExchangesInd) > 0) > 0
    msg = 'At least one non-exchange reaction has a lower bound higher than 0. This could make ANY task fail';
    warning(msg);
end

%Close all exchange reactions
previousLB = model.lb;
previousUB = model.ub;
model.lb(findRxnIDs(model,Exchange))=0;
model.ub(findRxnIDs(model,Exchange))=0;

%%%% Toy example oxidative phosphorylation
model = changeRxnBounds(model,'EX_glc__D_e',-10);
model = changeRxnBounds(model,'EX_o2_e',-60);
model = changeRxnBounds(model,'EX_h2o_e',60);
model = changeRxnBounds(model,'EX_co2_e',60);

%Set 1: GLCpts, PGI, PFK, FBA, TPI, GAPD, PGK, PGM, ENO,PYK, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD, H20t, co2t 28 Seedparameter = 1

%activeRxns= {'GLCpts', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK', 'PDH', 'CS', 'MDH','FUM','SUCDi','SUCOAS','AKGDH','ICDHyr','ACONTa','ACONTb','NADH16','CYTBD', 'O2t','ATPS4r', 'ATPM', 'NADTRHD','H2Ot', 'CO2t',};
%seedParameter =1;  %{[32]}    {[2.3234e+03]} with PDH orphan
%  {[32]}    {[2.0250e+03]} with 'PDH','SUCDi','ACONTa'orphan

%Set 2: GLCpts, PGI, PFK, FBA, TPI, GAPD, PGK, PGM, ENO,PYK, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD , FRD7,H20t, co2t 29 Seedparameter = 2

activeRxns= {'GLCpts', 'PGI', 'FRD7' 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK', 'PDH', 'CS', 'MDH','FUM','SUCDi','SUCOAS','AKGDH','ICDHyr','ACONTa','ACONTb','NADH16','CYTBD', 'O2t','ATPS4r', 'ATPM', 'NADTRHD','H2Ot', 'CO2t',};
seedParameter = 2;  % Doesnt turn on FRD7 for some reason


%Set 3: GLCpts, PGI, PFK, FBA, TPI, GAPD, PGK, PGM, ENO,PYK, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD , FRD7,H20t, co2t, FBP 30 Seedparameter = 3
activeRxns= {'GLCpts', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','FBP', 'PGK', 'PGM', 'ENO', 'PYK', 'PDH', 'CS', 'MDH','FUM','SUCDi','SUCOAS','AKGDH','ICDHyr','ACONTa','ACONTb','NADH16','CYTBD', 'O2t','ATPS4r', 'ATPM', 'NADTRHD','H2Ot', 'CO2t',};
seedParameter = 5;  % Doesnt turn on FRD7 for some reason, excludes ATPM except when:
% expressionRxns(16) = 607;
% expressionRxns(65) = 604;

%Set 4:  GLCpts, PGI, G6PDH2r, PGL, GND, RPI, RPE, TKT1, TALA, TKT2, GAPD, PGK, PGM, ENO, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD,FRD7 34 Seedparameter = 4

% activeRxns= {'GLCpts', 'PGI', 'G6PDH2r', 'TKT2', 'TALA' , 'PGL', 'RPE', 'RPI', 'TKT1',  'GND', 'GAPD','PGK', 'PGM', 'ENO', 'PDH', 'CS', 'MDH','FUM','SUCDi','SUCOAS','AKGDH','ICDHyr','ACONTa','ACONTb','NADH16','CYTBD', 'O2t','ATPS4r', 'ATPM', 'NADTRHD','H2Ot', 'CO2t',};
% seedParameter = 6;

%Set 5:  GLCpts, PGI, G6PDH2r, PGL, GND, RPI, RPE, TKT1, TALA, TKT2, GAPD, PGK, PGM, ENO, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD 33 Seedparameter = 5

activeRxns= {'GLCpts', 'PGI', 'G6PDH2r', 'TKT2', 'TALA' , 'PGL', 'RPE','RPI', 'TKT1',  'GND', 'GAPD','PGK', 'PGM', 'ENO', 'PDH', 'CS', 'MDH','FUM','SUCDi','SUCOAS','AKGDH','ICDHyr','ACONTa','ACONTb','NADH16','CYTBD', 'O2t','ATPS4r', 'ATPM', 'NADTRHD','H2Ot', 'CO2t',};
seedParameter = 6;

%Set 6: GLCpts, PGI, PFK, FBA, TPI, GAPD, PGK, PGM, ENO, PDH, CS, MDH,
%FUM, SUCDi, SucAOs, AKGHD, CDHyr, ACONTA, ACONTB, NADH16, CYTBD, O2T,
%ATPS4, ATPM, NADTRHD , FRD7,H20t, co2t, FBP, gludy, glusy, GLNS, PPC, ME2
%34 Seedparameter = 6


% activeRxns= {'PFK';'PFL';'PGK';'PGL';'ACALD';'AKGt2r';'PGM';'PPC';
%     'ACONTa';'ACONTb';'ATPM';'PPS';'ADK1';'AKGDH';'ATPS4r';
%      'CO2t';'CS';'SUCCt2_2';'CYTBD';'D_LACt2';'ENO';'ETOHt2r';'SUCDi';'SUCOAS';
%     'THD2';'FBA';'FUM';'G6PDH2r';'GAPD';'H2Ot';'ICDHyr'; 'MDH';'NADH16';'O2t';'PDH'};


orphanRxns = {'BIOMASS_Ecoli_core_w_GAM', 'TALA'}; %PDH','SUCDi','ACONTa'
expressionRxns = createTestData_ProvideAnswer(model,activeRxns, orphanRxns, seedParameter);
activeExpressionInd = findRxnIDs(model,activeRxns);
%expressionRxns(5) = 250;
%expressionRxns(16) = 607;
%expressionRxns(65) = 604;

varargin = parseSolverParameters('MILP');
epsilon = 0.01;

% Count the number of common variables
commonVars = intersect(orphanRxns, activeRxns);
lengthOrphansInActive =  numel(commonVars);

answer = (sum(expressionRxns(activeExpressionInd))+lengthOrphansInActive)/(length(activeExpressionInd)-lengthOrphansInActive);

espilon = 1;
% Find orphan ID's
orphanReactionsIds = find(expressionRxns ==-1);


% Indice the reactions of the model
edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
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



nnz_A = length(edgeIndices) * 28;
ind = zeros(nnz_A, 1);
j = zeros(nnz_A, 1);
v = zeros(nnz_A, 1);

count = 0;

for i = 1:length(edgeIndices)

    %Y+ LB -epsilon
    %A(i+size(S,1),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1);
    j(count) = edgeIndices(i);
    v(count) = 1;

    %A(i+size(S,1),i+size(S,2)) = lb(edgeIndices(i)) - epsilon;
    count = count + 1;
    ind(count) = i + size(S, 1);
    j(count) = i + size(S, 2);
    v(count) = lb(edgeIndices(i)) - epsilon;

    %Y- UB; + epsilon
    %A(i+size(S,1)+length(edgeIndices),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + length(edgeIndices);
    j(count) = edgeIndices(i);
    v(count) = 1;

    %A(i+size(S,1)+length(edgeIndices),reactionLength+i+size(S,2)) = ub(edgeIndices(i)) + epsilon;
    count = count + 1;
    ind(count) = i + size(S, 1) + length(edgeIndices);
    j(count) = reactionLength + i + size(S, 2);
    v(count) = ub(edgeIndices(i)) + epsilon;




    %Y+ + Y- <= 0
    %A(i+size(S,1)+2*length(edgeIndices),i + 1*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 2 * length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+2*length(edgeIndices),i + 2*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 2 * length(edgeIndices);
    j(count) = i + 2 * size(S, 2);
    v(count) = 1;


    %Y+ + Y- + Yprime = 1
    %A(i+size(S,1)+3*length(edgeIndices),i + 1*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 3 * length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+3*length(edgeIndices),i + 2*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 3 * length(edgeIndices);
    j(count) = i + 2 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+3*length(edgeIndices),i + 3*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 3 * length(edgeIndices);
    j(count) = i + 3 * size(S, 2);
    v(count) = 1;


    %% Z Plus blocks
    %Z>=0  'G' 0
    %A(i+size(S,1)+4*length(edgeIndices),i + 4*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 4 * length(edgeIndices);
    j(count) = i + 4 * size(S, 2);
    v(count) = 1;

    %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
    %A(i+size(S,1)+5*length(edgeIndices),i + 1*size(S,2)) = -M;
    count = count + 1;
    ind(count) = i + size(S, 1) + 5* length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = -M;

    %A(i+size(S,1)+5*length(edgeIndices),i + 4*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 5 * length(edgeIndices);
    j(count) = i + 4 * size(S, 2);
    v(count) = 1;


    % M*y + z-t <= M 'L' M
    %A(i+size(S,1)+6*length(edgeIndices),i + 1*size(S,2)) = M;
    count = count + 1;
    ind(count) = i + size(S, 1) + 6 * length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = M;

    %A(i+size(S,1)+6*length(edgeIndices),i + 4*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 6 * length(edgeIndices);
    j(count) = i + 4 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+6*length(edgeIndices),end) = -1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 6 * length(edgeIndices);
    j(count) = size(S, 2) + 5 * length(edgeIndices) + 1; % double check this %A_ColumnLength = size(S,2) + 2*reactionLength + 2*reactionLength + 1*reactionLength  + 1;
    v(count) = -1;

    % M*y - z + t <= M 'L' M
    %A(i+size(S,1)+7*length(edgeIndices),i + 1*size(S,2)) = M;
    count = count + 1;
    ind(count) = i + size(S, 1) + 7 * length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = M;

    %A(i+size(S,1)+7*length(edgeIndices),i + 4*size(S,2)) = -1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 7 * length(edgeIndices);
    j(count) = i + 4 * size(S, 2);
    v(count) = -1;

    %A(i+size(S,1)+7*length(edgeIndices),end) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 7* length(edgeIndices);
    j(count) = size(S, 2) + 5 * length(edgeIndices) + 1;
    v(count) = 1;

    %% Z minus blocks
    %Z>=0  'G' 0
    %A(i+size(S,1)+8*length(edgeIndices),i + 4*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+8*length(edgeIndices);
    j(count) = i + 4*size(S,2);
    v(count) = 1;

    %Z<= M*Y == Z-(M*Y) <= 0 'L' 0
    %A(i+size(S,1)+9*length(edgeIndices),i + 2*size(S,2)) = -M;
    count = count + 1;
    ind(count) = i+size(S,1)+9*length(edgeIndices);
    j(count) = i + 2*size(S,2);
    v(count) = -M;

    %A(i+size(S,1)+9*length(edgeIndices),i + 5*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+9*length(edgeIndices);
    j(count) = i + 5*size(S,2);
    v(count) = 1;


    % M*y + z-t <= M 'L' M
    %A(i+size(S,1)+10*length(edgeIndices),i + 2*size(S,2)) = M;
    count = count + 1;
    ind(count) = i+size(S,1)+10*length(edgeIndices);
    j(count) = i + 2*size(S,2);
    v(count) = M;

    %A(i+size(S,1)+10*length(edgeIndices),i + 5*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+10*length(edgeIndices);
    j(count) = i + 5*size(S,2);
    v(count) = 1;

    %A(i+size(S,1)+10*length(edgeIndices),end) = -1;
    count = count + 1;
    ind(count) = i+size(S,1)+10*length(edgeIndices);
    j(count) = size(S, 2) + 5 * length(edgeIndices) + 1;
    v(count) = -1;

    % M*y - z + t <= M 'G' -M
    %A(i+size(S,1)+11*length(edgeIndices),i + 2*size(S,2)) = M;
    count = count + 1;
    ind(count) = i+size(S,1)+11*length(edgeIndices);
    j(count) = i + 2*size(S,2);
    v(count) = M;

    %A(i+size(S,1)+11*length(edgeIndices),i + 5*size(S,2)) = -1;
    count = count + 1;
    ind(count) = i+size(S,1)+11*length(edgeIndices);
    j(count) = i + 5*size(S,2);
    v(count) = -1;

    %A(i+size(S,1)+11*length(edgeIndices),end) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+11*length(edgeIndices);
    j(count) = size(S, 2) + 5 * length(edgeIndices) + 1;
    v(count) = 1;


    %% Other blocks
    % Yprime LB
    %A(i+size(S,1)+12*length(edgeIndices),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+12*length(edgeIndices);
    j(count) = edgeIndices(i);
    v(count) = 1;

    %A(i+size(S,1)+12*length(edgeIndices),i+ 3*size(S,2)) = lb(edgeIndices(i));
    count = count + 1;
    ind(count) = i+size(S,1)+12*length(edgeIndices);
    j(count) = i+ 3*size(S,2);
    v(count) = lb(edgeIndices(i));


    % Yprime UB
    %A(i+size(S,1)+13*length(edgeIndices),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+13*length(edgeIndices);
    j(count) = edgeIndices(i);
    v(count) =  1;

    %A(i+size(S,1)+13*length(edgeIndices),i+ 3*size(S,2)) = ub(edgeIndices(i));
    count = count + 1;
    ind(count) = i+size(S,1)+13*length(edgeIndices);
    j(count) = i+ 3*size(S,2);
    v(count) =  ub(edgeIndices(i));

    % Z = 1 single row (twice as long for z+ and z-)
    %%% Set orphan to 0 in addition to c_z (obj.coef) of orphan to 0
    if ismember(i, orphanReactionsIds)
        % A(size(S,1)+14*length(edgeIndices) + 1, i + 4*size(S,2)) = 0;
        count = count + 1;
        ind(count) = size(S,1)+14*length(edgeIndices) + 1;
        j(count) = i + 4*size(S,2);
        v(count) =  0;

        %A(size(S,1)+14*length(edgeIndices) + 1, i + 5*size(S,2)) = 0;
        count = count + 1;
        ind(count) = size(S,1)+14*length(edgeIndices) + 1;
        j(count) = i + 5*size(S,2);
        v(count) =  0;

    else
        %A(size(S,1)+14*length(edgeIndices) + 1, i + 4*size(S,2)) = 1;
        count = count + 1;
        ind(count) = size(S,1)+14*length(edgeIndices) + 1;
        j(count) = i + 4*size(S,2);
        v(count) =  1;

        %A(size(S,1)+14*length(edgeIndices) + 1, i + 5*size(S,2)) = 1;
        count = count + 1;
        ind(count) = size(S,1)+14*length(edgeIndices) + 1;
        j(count) = i + 5*size(S,2);
        v(count) =  1;

    end
end
%% 2 single entries for 0>=t>=M
% Question: Should this not be included in lb/ub variable of T?
%A(size(S,1)+14*length(edgeIndices) + 2, end) = 1;
count = count + 1;
ind(count) = size(S,1)+14*length(edgeIndices) + 2;
j(count) = size(S, 2) + 5 * length(edgeIndices) + 1;
v(count) =  1;

%A(size(S,1)+14*length(edgeIndices) + 3, end) = 1;
count = count + 1;
ind(count) = size(S,1)+14*length(edgeIndices) + 3;
j(count) = size(S, 2) + 5 * length(edgeIndices) + 1;
v(count) =  1;

ind = ind(1:count);
j = j(1:count);
v = v(1:count);

A = sparse(ind, j, v, size(S, 1) + 14 * length(edgeIndices) + 3, reactionLength + 5 * size(S, 2) + 1);
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

spy(A);
