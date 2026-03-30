function expressionK = findTrueKValue_LP7_Approach(solutionMILP, expValueMeanAdjustedTask, model, epsilon)


reactionLength = length(model.rxns);
SMatrix = model.S;
ORindices = expValueMeanAdjustedTask == -1;
expressionIndices = expValueMeanAdjustedTask ~= -1;
kRowNumber = 1+size(SMatrix,1)+5*(reactionLength); 

MILP_v_variables = solutionMILP.cont(1:length(model.rxns));
MILP_fluxes = MILP_v_variables(find(MILP_v_variables>epsilon-(epsilon/100))&& expValueMeanAdjustedTaskIrrevKapprox>0 );
MILP_Expression = expValueMeanAdjustedTaskIrrevKapprox(find(MILP_v_variables>epsilon-(epsilon/100)) && expValueMeanAdjustedTaskIrrevKapprox>0 );

sumExpressionValues = sum(MILP_Expression);
expressionK = (sumExpressionValues)/length(MILP_Expression);


end