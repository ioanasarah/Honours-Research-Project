function MILPproblem = adaptKMILPProblemAverageExpressionLoops(MILPproblem, model, expressionRxns, k);

S = model.S;

edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end
% 
% ind = [];
% j = [];
v1 = [];
v2 = [];
count = 0;

for i = 1:length(edgeIndices)

    count = count + 1;
    % ind(count) = 1+size(S,1)+7*length(edgeIndices);
    % j(count) = i+ size(S,2);
    if(expressionRxns(i) ~= -1)
        v1(i) = expressionRxns(i)-k;
    else
        v1(i) = 0;
    end
    
    count = count + 1;
    % ind(count) = 1+size(S,1)+7*length(edgeIndices);
    % j(count) = i+ 2*size(S,2);
    if(expressionRxns(i) ~= -1)
        v2(i) = expressionRxns(i)-k;
    else
        v2(i) = 0;
    end

end

v = [v1,v2];

% Replace the desired part of the original matrix with the adjusted part
MILPproblem.A(1+size(S,1)+5*length(edgeIndices), 1+1*size(S,2):length(edgeIndices)+2*size(S,2)) = v;
spy(MILPproblem.A)
