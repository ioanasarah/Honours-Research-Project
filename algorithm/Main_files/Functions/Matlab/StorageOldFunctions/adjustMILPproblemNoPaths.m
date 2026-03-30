function MILPproblem = adjustMILPproblemNoPaths(MILPproblem,expressionRxns,ktest)

ExpressionMinusKVector = zeros((size(MILPproblem.A,2)/2),1);
for i= 1:(size(MILPproblem.A,2)/2)
    ExpressionMinusKVector(i) = expressionRxns(i) - ktest;
end

for i = 1:(size(MILPproblem.A,2)/2)
    MILPproblem.A(end,i+(size(MILPproblem.A,2)/2)) = ExpressionMinusKVector(i);    
end


 