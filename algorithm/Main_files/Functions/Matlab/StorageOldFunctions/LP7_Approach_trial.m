tic
[ATask, fastCCModel, V] = fastcc(taskModelOldCore,epsilon);
toc
Anew = zeros(size(taskModelOldCore.rxns));
indices = ismember(taskModelOldCore.rxns, taskModelOldCore.rxns(ATask));
Anew(indices) =1;
reactionsNotPartOfTask = Anew==0;
newReversibleModelForKApproxmiation = removeRxns(taskModelOldCore, taskModelOldCore.rxns(reactionsNotPartOfTask));
expValueMeanAdjustedTaskReversibleKApproximation = expValueMeanAdjustedTaskOld(Anew ==1);


% save('newReversibleModelForKApproxmiation', 'newReversibleModelForKApproxmiation')
% save('expValueMeanAdjustedTaskReversibleKApproximation', 'expValueMeanAdjustedTaskReversibleKApproximation')
%load('expValueMeanAdjustedTaskReversibleKApproximation')
[newIrrevModelKapprox, matchRev, rev2irrevIrrevKapprox, irrev2rev] = convertToIrreversible(newReversibleModelForKApproxmiation);
expValueMeanAdjustedTaskIrrevKapprox = expValueMeanAdjustedTaskReversibleKApproximation;

for g=1:size(rev2irrevIrrevKapprox,1)
    if size(rev2irrevIrrevKapprox{g},2) == 2
        expValueMeanAdjustedTaskIrrevKapprox(rev2irrevIrrevKapprox{g}(2)) = expValueMeanAdjustedTaskReversibleKApproximation(rev2irrevIrrevKapprox{g}(1));
    end
end


load("newIrrevModelKapprox.mat")
load("expValueMeanAdjustedTaskIrrevKapprox")
load("rev2irrevIrrevKapprox.mat")


model = newIrrevModelKapprox;
J = 1:length(model.rxns);
LPproblem = buildLPproblemFromModel(model);

% Set median expression
expression = expValueMeanAdjustedTaskIrrevKapprox ;
expression(expression<=0) = median(expression(expression>0));

nj = numel(J);
[m,n] = size(model.S);
[m2,n2] = size(LPproblem.A);

% objective
f = -[zeros(n2,1); ones(nj,1);zeros(n2,1)];

kTolerance = 0.001;

kFeasible = 0;
kinfeasible = max(expression);
kDifference = kInfeasible - kFeasible;
Kconstant = 0.6;
ktest = kFeasible + (Kconstant * kDifference);


while kDifference > kTolerance
    % equalities
    Aeq = [LPproblem.A, sparse(m2,nj),sparse(m2,nj) ];
    beq = LPproblem.b;

    %ktest = 4.4396;%4.4396

    % inequalities
    Ij = sparse(nj,n2);
    Ij(sub2ind(size(Ij),(1:nj)',J(:))) = -1;
    Aineq = [sparse([Ij, speye(nj)]), sparse(nj,nj)] ;


    AdditionalConstraint = sparse(nj,n2);
    indices = sub2ind(size(AdditionalConstraint), (1:nj)', J(:));
    AdditionalConstraint(indices) = model.lb;
    additionalInEquality= sparse([speye(nj),sparse(nj,nj),AdditionalConstraint]);

    AdditionalConstraint2 = sparse(nj,n2);
    indices2 = sub2ind(size(AdditionalConstraint2), (1:nj)', J(:));
    AdditionalConstraint2(indices2) = model.ub;
    additionalInEquality2= sparse([speye(nj),sparse(nj,nj),AdditionalConstraint2]);

    additionalInEquality3= sparse([sparse(nj,nj),speye(nj),speye(nj)]);

    additionalInEquality4 = sparse([sparse(nj,nj),sparse(nj,nj),sparse(nj,nj)]);

    % Iterate through rev2irrevIrrevKapprox
    for i = 1:numel(rev2irrevIrrevKapprox)
        % Check if the element contains two values
        if isequal(numel(rev2irrevIrrevKapprox{i}), 2)
            % Get the values
            val1 = rev2irrevIrrevKapprox{i}(1);
            val2 = rev2irrevIrrevKapprox{i}(2);

            % Set the corresponding elements in additionalInEquality4 to 1
            additionalInEquality4(val1, nj + val1) = 1;
            additionalInEquality4(val1, nj + val2) = 1;
        else
            additionalInEquality4(i, nj+i)=0;
        end

    end
    spy(additionalInEquality4)


    bineq = zeros(nj,1);
    additionalEquality = model.lb;
    additionalEquality2 = model.ub;
    additionalEquality3 = ones(nj,1);
    additionalEquality4 = ones(nj,1);
    AineqK = sparse([zeros(n2,1);expression-ktest;zeros(n2,1) ]');
    bineqK = 0;

    % bounds
    %lb = [LPproblem.lb; -inf*ones(nj,1)]; %% Original
    lb = [LPproblem.lb; zeros(nj,1) ; zeros(nj,1)];
    ub = [LPproblem.ub; ones(nj,1)*epsilon ; ones(nj,1)*epsilon];

    % Set up LP problem
    LP7problem.A=[Aeq;Aineq;additionalInEquality;additionalInEquality2;additionalInEquality3;additionalInEquality4;AineqK];
    LP7problem.b=[beq;bineq;additionalEquality;additionalEquality2;additionalEquality3;additionalEquality4; bineqK];
    LP7problem.lb=lb;
    LP7problem.ub=ub;
    LP7problem.c=f;
    LP7problem.osense=-1;%minimise
    LP7problem.csense = [LPproblem.csense; repmat('L',nj,1); repmat('G',nj,1);repmat('L',nj,1); repmat('E',nj,1); repmat('L',nj,1); repmat('G',1,1) ];
    LP7problem.vartype = [char(repmat("C",1*length(model.rxns),1));char(repmat("B",2*length(model.rxns),1))];

    spy(LP7problem.A)

    % tic
    % solutionLP = solveCobraLP(LP7problem);
    % toc

    tic
    LP7problem = addLoopLawConstraints(LP7problem, model, 1:length(model.rxns));
    toc

    fprintf("%f",ktest)
    tic
    rng(1)
    solutionMILP = solveCobraMILP(LP7problem)
    toc

    if solutionMILP.stat ==1 || solutionMILP.stat ==10

        MILP_v_variables = solutionMILP.cont(1:length(model.rxns));
        MILP_fluxes = MILP_v_variables(find(MILP_v_variables>epsilon-(epsilon/10)));
        MILP_zVar = MILP_z_variables(find(MILP_v_variables>epsilon-(epsilon/10)));
        MILP_amount_Fluxes = length(find(MILP_v_variables>epsilon-(epsilon/10)));
        MILP_rxns= model.rxns(find(MILP_v_variables>epsilon-(epsilon/10)));
        MILP_rxnNames= model.rxnNames(find(MILP_v_variables>epsilon-(epsilon/10)));
        MILP_rxnFormulas = printRxnFormula(model,model.rxns(find(MILP_v_variables>epsilon-(epsilon/10))),false );
        MILP_Expression = expValueMeanAdjustedTaskIrrevKapprox(find(MILP_v_variables>epsilon-(epsilon/10)));

        kFeasible = ktest;
    elseif solutionMILP.stat ==0
        kinfeasible = ktest;
    end
    kDifference = kinfeasible - kFeasible;
    disp(kDifference)
    ktest = kFeasible + (Kconstant * kDifference);
end
% spy(milpproblemfromtask1.A)
% tic
% solutionMILP = solveCobraMILP(LP7problem)
% toc
%
% MILP_v_variables = solutionMILP.cont(1:length(model.rxns));
% MILP_z_variables = solutionMILP.cont(length(model.rxns)+1:2*length(model.rxns));
% % to inspect:
% %Avector = full(LP7problem.A(25864,1:2*length(model.rxns))');
%
%
% milpproblemfromtask1 = load("LastCreatedMILPproblem_Task1Sample1_date_2023-11-28___17-55-20.mat").MILPproblem;
% spy(milpproblemfromtask1.A)


