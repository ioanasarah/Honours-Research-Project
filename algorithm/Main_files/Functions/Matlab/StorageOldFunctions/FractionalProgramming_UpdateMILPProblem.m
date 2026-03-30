function MILPproblem = FractionalProgramming_UpdateMILPProblem(MILPproblem,expressionRxns,ktest)

ExpressionMinusKVector = zeros((size(MILPproblem.A,2)/3),1);
disp(length(ExpressionMinusKVector))
for i= 1:(size(MILPproblem.A,2)/3)
    ExpressionMinusKVector(i) = expressionRxns(i) - ktest;
end

ExpressionMinusKVector1 = [ExpressionMinusKVector; ExpressionMinusKVector];

for i = 1:(size(ExpressionMinusKVector))
    MILPproblem.A(end,i+(size(MILPproblem.A,2)/3)) = ExpressionMinusKVector(i);    
end


 