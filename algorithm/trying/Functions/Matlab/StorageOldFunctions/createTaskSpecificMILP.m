function taskSpecificMILP = createTaskSpecificMILP(coreMILPProblem, taskModelOldCore, expressionRxns,epsilon)

MILPproblem = coreMILPProblem;
% S matrix
% UB 
% LB
% Expression
k=0;

S = taskModelOldCore.S;

lb = taskModelOldCore.lb;
ub = taskModelOldCore.ub;

% len1=length(taskModelOldCore.lb);
% b = setdiff(taskModelOldCore.lb, coreMILPProblem.lb(1:len1));
% idx = find(~ismember(coreMILPProblem.lb(1:len1),taskModelOldCore.lb));
% differentElement = taskModelOldCore.lb(idx);

reactionLength = length(taskModelOldCore.rxns);

edgeIndices = zeros(length(taskModelOldCore.rxns),1);
for i = 1:length(taskModelOldCore.rxns)
    edgeIndices(i) = i;
end
% 
% ind = [];
% j = [];
v1 = [];
v2 = [];
count = 0;

% adjust ub lb 
MILPproblem.lb(1:size(S,2)) = lb;
MILPproblem.ub(1:size(S,2)) = ub;

% replace Smatrix 
MILPproblem.A(1:size(S,1),1:size(S,2)) = S;

% adjust ub lb A matrix values

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

    %Y+ + Y- + Yprime = 1
    %A(i+size(S,1)+3*length(edgeIndices),i + 1*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 2 * length(edgeIndices);
    j(count) = i + 1 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+3*length(edgeIndices),i + 2*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 2 * length(edgeIndices);
    j(count) = i + 2 * size(S, 2);
    v(count) = 1;

    %A(i+size(S,1)+3*length(edgeIndices),i + 3*size(S,2)) = 1;
    count = count + 1;
    ind(count) = i + size(S, 1) + 2 * length(edgeIndices);
    j(count) = i + 3 * size(S, 2);
    v(count) = 1;

    %% Other blocks
    % Yprime LB
    %A(i+size(S,1)+12*length(edgeIndices),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+3*length(edgeIndices);
    j(count) = edgeIndices(i);
    v(count) = 1;

    %A(i+size(S,1)+12*length(edgeIndices),i+ 3*size(S,2)) = lb(edgeIndices(i));
    count = count + 1;
    ind(count) = i+size(S,1)+3*length(edgeIndices);
    j(count) = i+ 3*size(S,2);
    v(count) = lb(edgeIndices(i));


    % Yprime UB
    %A(i+size(S,1)+13*length(edgeIndices),edgeIndices(i)) = 1;
    count = count + 1;
    ind(count) = i+size(S,1)+4*length(edgeIndices);
    j(count) = edgeIndices(i);
    v(count) =  1;

    %A(i+size(S,1)+13*length(edgeIndices),i+ 3*size(S,2)) = ub(edgeIndices(i));
    count = count + 1;
    ind(count) = i+size(S,1)+4*length(edgeIndices);
    j(count) = i+ 3*size(S,2);
    v(count) =  ub(edgeIndices(i));

    %Single Row expression-k>0    
    count = count + 1;
    ind(count) = 1+size(S,1)+5*length(edgeIndices);
    j(count) = i+ size(S,2);
    if(expressionRxns(i) ~= -1)
        v(count) = expressionRxns(i)-k;
    else
        v(count) = 0;
    end

    count = count + 1;
    ind(count) = 1+size(S,1)+5*length(edgeIndices);
    j(count) = i+ 2*size(S,2);
    if(expressionRxns(i) ~= -1)
        v(count) = expressionRxns(i)-k;
    else
        v(count) = 0;
    end

end
ind = ind(1:count);
j = j(1:count);
v = v(1:count);

A = sparse(ind, j, v, size(S, 1) + 5* length(edgeIndices) + 1, reactionLength + 3 * size(S, 2));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end


 MILPproblem.A(1:size(S, 1) + 5* length(edgeIndices) + 1, 1:reactionLength + 3 * size(S, 2)) = A;

%Adjust B
MILPproblem.b(size(S,1)+1:size(S,1)+1*length(edgeIndices)) = lb;
MILPproblem.b(size(S,1)+1+1*length(edgeIndices):size(S,1)+2*length(edgeIndices)) = ub;
MILPproblem.b(size(S,1)+1+3*length(edgeIndices):size(S,1)+4*length(edgeIndices)) = lb;
MILPproblem.b(size(S,1)+1+4*length(edgeIndices):size(S,1)+5*length(edgeIndices)) = ub;

taskSpecificMILP = MILPproblem;
spy(MILPproblem.A)
end