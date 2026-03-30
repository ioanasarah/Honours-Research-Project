function MILPproblem = adjustMILPproblem(MILPproblem,expressionRxns,currentPathLength,ktest)

ExpressionMinusKVector = zeros((size(MILPproblem.A,2)/2),1);
for i= 1:(size(MILPproblem.A,2)/2)
    ExpressionMinusKVector(i) = expressionRxns(i) - ktest;
end

for i = 1:(size(MILPproblem.A,2)/2)
    MILPproblem.A(end-1,i+(size(MILPproblem.A,2)/2)) = ExpressionMinusKVector(i);    
end
MILPproblem.b(end) = currentPathLength;

