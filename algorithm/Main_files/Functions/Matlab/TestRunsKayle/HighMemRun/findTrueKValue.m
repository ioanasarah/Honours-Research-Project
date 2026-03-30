function expressionK = findTrueKValue(solution, MILPproblem, expValueMeanAdjustedTask, alpha, beta,reactionLength,ktest,SMatrix,TransportReactionPenalty,transportRxnIDs);




ORindices = expValueMeanAdjustedTask == -1;
expressionIndices = expValueMeanAdjustedTask ~= -1;
kRowNumber = 1+size(SMatrix,1)+5*(reactionLength); 




logicalIndicesExpression1 = find((solution.int(1:reactionLength) == 1) & (expressionIndices == 1));
logicalIndicesExpression = logicalIndicesExpression1+ (reactionLength);

ExpressionNonTransportReactions = find((solution.int(1:reactionLength) == 1) & (transportRxnIDs == 0) &(expressionIndices == 1)); %these get a bonus (penalty to others)


logicalIndicesOR1 = find((solution.int(reactionLength+1:2*reactionLength) == 1) & (ORindices == 1));
logicalIndicesOR = logicalIndicesOR1+(reactionLength*2);
valuesYPrimeOR = sum(full(MILPproblem.A(kRowNumber,logicalIndicesOR)));

% inactiveYprimeORTransportReactions1 = find((solution.int(1:reactionLength) == 1) & (transportRxnIDs == 1)&(ORindices == 1));
% inactiveYprimeORTransportReactions = inactiveYprimeORTransportReactions1+(reactionLength*2);
% valuesinactiveYprimeORTransportReactions = sum(full(MILPproblem.A(kRowNumber,logicalIndicesOR)));

sumExpressionValues = 2*alpha*sum(expValueMeanAdjustedTask(logicalIndicesExpression1));
expressionK = (( (valuesYPrimeOR + ( sumExpressionValues+(TransportReactionPenalty*length(ExpressionNonTransportReactions)) ) )  / (2*alpha) ) )   /length(logicalIndicesExpression);


end