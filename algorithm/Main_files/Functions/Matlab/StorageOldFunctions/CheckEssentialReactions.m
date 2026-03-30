%%%%
% Clear vars for next run
clearvars

% Necessary evewry time
%initCobraToolbox

% Load some model
model = readCbModel('e_coli_core.mat');

%Set epsilon
epsilon = 0.001; 

% % Create irr model
% modelIrrev = convertToIrreversible(Oldmodel);

%[essentialRxns,usedRxns] = essentialRxnsTasks2023_1(model,1);

% Add in random expression data
expressionRxns = abs(randn(length(model.rxns),1));

% Create lineIndices
edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end



saveUsed = true;


relaxation = 1e-06;

[solMin, modelIrrevFMOri]= minimizeModelFlux(model,'min'); %Compute the minimal set of reactions
if relaxation == Inf
    modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',Inf,'u');
    modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',1e-06,'l');%Will prevent negative or null netFlux
else
    modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',solMin.f + relaxation,'u');
    modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',max(solMin.f - relaxation, 1e-06),'l');%Will prevent negative or null netFlux
end
%modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',1e-06,'l');
%Define the list of reactions to test
rxnsToCheck=modelIrrevFMOri.rxns(abs(solMin.x)>10^-6);
%[transRxns] = findTransRxns(model,1);
%rxnsToCheck=rxnsToCheck_raw(~ismember(rxnsToCheck_raw,transRxns));

% Loop that set to 0 each reaction to test and check if the problem
% still has a solution
%essentialRxns={};
isEssential = false(size(rxnsToCheck));
if saveUsed
    usedRxns = false(size(modelIrrevFMOri.rxns));
end
for i=1:numel(rxnsToCheck)
    %disp(num2str(i))
    modelIrrevFM = modelIrrevFMOri;
    modelIrrevFM.lb(findRxnIDs(modelIrrevFM,rxnsToCheck(i)))=0;
    modelIrrevFM.ub(findRxnIDs(modelIrrevFM,rxnsToCheck(i)))=0;
    modelIrrevFM.csense(1:length(modelIrrevFM.mets),1) = 'E';
    modelIrrevFM.osense = -1;
    modelIrrevFM.A=modelIrrevFM.S;
    sol=solveCobraLP(modelIrrevFM);
    disp(sol.stat)
    if sol.stat==0 || isempty(sol.full)
        isEssential(i) = true;

    elseif saveUsed
        tosavevarr = find(sol.full);
        usedRxns(sol.full > 1e-06) = true;
    end

end
essentialRxns = rxnsToCheck(isEssential);
rxns_kept=unique(essentialRxns);
rxns_final=cell(size(rxns_kept));
if saveUsed
    usedRxns = modelIrrevFMOri.rxns(usedRxns);
    usedRxns = [usedRxns;essentialRxns;rxnsToCheck];
    usedRxns = unique(usedRxns);
end
%% Analysis part
for i=1: length(rxns_kept)
    string=rxns_kept{i};
    if strcmp('_f', string(end-1:end))==1
        rxns_final{i}= string(1:end-2);
    elseif strcmp('_b', string(end-1:end))==1
        rxns_final{i}= string(1:end-2);
    elseif strcmp('_r', string(end-1:end))==1
        rxns_final{i}= string(1:end-2);
    else
        rxns_final{i}=string;
    end
end
essentialRxns=unique(rxns_final);
essentialRxns(strcmp(essentialRxns, 'netFlux')) = [];
isTemp = contains(essentialRxns,'temporary');
if sum(isTemp) ~=length(isTemp)
    essentialRxns = essentialRxns(~isTemp);
end

%% Analysis part2
if saveUsed
    usedFinal = cell(size(usedRxns));
    for i=1: length(usedRxns)
        string=usedRxns{i};
        if strcmp('_f', string(end-1:end))==1
            usedFinal{i}= string(1:end-2);
        elseif strcmp('_b', string(end-1:end))==1
            usedFinal{i}= string(1:end-2);
        elseif strcmp('_r', string(end-1:end))==1
            usedFinal{i}= string(1:end-2);
        else
            usedFinal{i}=string;
        end
    end
    usedRxns=unique(usedFinal);
    usedRxns(strcmp(usedRxns, 'netFlux')) = [];
    isTemp = contains(usedRxns,'temporary');
    if sum(isTemp) ~=length(isTemp)
        usedRxns = usedRxns(~isTemp);
    end
end




















% % Get S (irr is sparse, instead of double, shouldn't matter but might for
% % null calc)
% S = model.S;
% SForNull = full(S);
% NullOfIrrMod = null(SForNull);
% lb = model.lb;  % +0.01;
% ub = model.ub;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating Toy model  
% Toy Data
edgeData = [1 1 1; 1 2 1; 2 3 2.3;  3 4 3.2;   4 5 5.2;    5 6 22;    2 7 1;    7 8 1;    8 6 1; 6 9 1; 9 9 1];
nNodes = 9;
nEdges = size(edgeData, 1);
% 
% %Create S for Toy model
% % S = zeros(nNodes+1, nEdges+1); %with pseudorxn
% S = zeros(nNodes, nEdges);
% for i = 1:nEdges
%     startNode = edgeData(i, 1);
%     endNode = edgeData(i, 2);
%     weight = edgeData(i, 3);
% 
%     S(startNode, i) = -1;
%     S(endNode, i) = 1;
% end
% 
% %Add pseudoMet/Reaction
% for i = 1:size(S,2)
%     S(1,i)= 1;
%     % S(end,i) = 1;
% end
% % S(end,end) = -1;
% % S(1,1) = -1;
% % S(end,1) = 0;
% % S(1,end) = 0;
% %S(1,end+1) = -1;
% 
% % Create Rxns
% expressionRxns = edgeData(1:end, 3);
% %expressionRxns(end+1) = 1; %add reaction data for pseudoreaction
% % Create edgeIndices
% edgeIndices = zeros(size(S,2),1);
% for i = 1:size(S,2)
%     edgeIndices(i) = i;
% end
% 
% lb = zeros(length(edgeIndices),1);
% ub = 1000* ones(length(edgeIndices),1);  % +0.01;
% lb(1) = -1;
% % lb(end) = 1;
%  ub(1) = 0;
% % ub(end) = 1;
% % End creation Toy Model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %set K for testing
% k = 1; %for testing
% 
% % Create ExpressionMinusKVector 
% ExpressionMinusKVector = zeros(length(edgeIndices),1);
% for i= 1:length(edgeIndices)
%     ExpressionMinusKVector(i) = expressionRxns(i) - k;
% end
% 
% % Creating A matrix
% A = sparse(size(S,1)+2, size(S,2)+length(edgeIndices));
% [m,n,s] = find(S);
% for i = 1:length(m)
%     A(m(i),n(i)) = s(i);
% end
% 
% %Fill A matrix with values
% for i = 1:length(edgeIndices)
%     A(1+size(S,1),edgeIndices(i)) = 0; %tried both
%     A(1+size(S,1),i+size(S,2)) = ExpressionMinusKVector(i);% expressionRxns(i) - k ;
% 
%     A(2+size(S,1),i) = 1;
%     A(2+size(S,1),i+size(S,2)) = 1;
% end
% 
% spy(A)
% 
% % Creating lb and ub
% lb_y = zeros(length(edgeIndices),1);
% ub_y = ones(length(edgeIndices),1);
% lb = [lb;lb_y;];
% ub = [ub;ub_y;];
% 
% % Creating vartype
% vartype1(1:size(S,2),1) = 'C';
% vartype2(1:length(edgeIndices),1) = 'B';
% vartype = [vartype1;vartype2];
% 
% % Creating c
% c_v = 0 * ones(size(S,2),1); %this should be adapted to activity; see following lines
% c_v(end) = 1;%necessary to be able to actually minimize flux? See minimizeModelFlux
% %c_v(1) = 1; %set begin reaction to 1
% c_y = ones(length(edgeIndices),1);% set coefficients of decision variables to 1
% c = [c_v;c_y;];
% 
% %Creating B; b probably needs to be remade in binary search loop
% b_s = zeros(size(S,1),1); %need to add in k here, unsure if all values or single value
% % lb_indices = lb(edgeIndices); %  unsure if correct
% % ub_indices = ub(edgeIndices);
% kComparison = 0;
% EnforceMinimumFlux = 5;
% b = [b_s;kComparison; EnforceMinimumFlux];
% 
% % Creating csense
% csense1(1:size(S,1)) = 'E';
% csense4(1) = 'G'; % sum(ExpressionMinusKVector *yi) => 0 (kComparison)
% csense5(1) = 'G'; % 
% csense = [csense1  csense4  csense5];
% 
% % create MILPproblem
% MILPproblem.A = A; % A Matrix
% MILPproblem.b = b; % Constraint vector
% MILPproblem.c = c; % Objective coeficients (v and Y)
% MILPproblem.lb = lb; % Lower bounds 
% MILPproblem.ub = ub; % Upper bounds
% MILPproblem.csense = csense; % Equality/lower/higher 
% MILPproblem.vartype = vartype; % Vartype ('C', 'B', 'I')
% MILPproblem.osense = 1; %(-1 max, +1 min)
% MILPproblem.x0 = []; % Initial guess value
% 
% 
% % % Set extra params for optimizeCbModel and run optimzeCbModel
% % osenseStr = 'min';
% % minNorm = 0;
% % disp("Running fluxMinimization")
% % solutionCB = optimizeCbModel(modelIrrev, osenseStr, minNorm);
% % disp(solutionCB.f);
% % disp(solutionCB.stat);
% 
% %run SolveCobraMILP
% disp("Running MILP")
% solution = solveCobraMILP(MILPproblem);
% disp(solution.obj);
% disp(solution.stat);
