function MILPproblem = constructMILPproblem(model)

edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end
epsilon = 0.001;
% Get S
S = model.S;
lb = model.lb; 
ub = model.ub;

A = sparse(size(S,1)+2*length(edgeIndices)+1, size(S,2)+length(edgeIndices));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

% Create EMPTY ExpressionMinusKVector 
ExpressionMinusKVector = zeros(length(edgeIndices),1);

%Fill A matrix with values
for i = 1:length(edgeIndices)
    A(i+size(S,1),edgeIndices(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(edgeIndices(i))-epsilon;

    A(i+size(S,1)+length(edgeIndices),edgeIndices(i)) = 1;
    A(i+size(S,1)+length(edgeIndices),i+size(S,2)) = 0;

    A(1+size(S,1)+2*length(edgeIndices),edgeIndices(i)) = 0; 
    A(1+size(S,1)+2*length(edgeIndices),i+size(S,2)) = ExpressionMinusKVector(i);
end

% Creating lb and ub
lb_y = zeros(length(edgeIndices),1);
ub_y = ones(length(edgeIndices),1);
lb = [lb;lb_y;];
ub = [ub;ub_y;];

% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:length(edgeIndices),1) = 'B';
vartype = [vartype1;vartype2];

% Creating c
c_v = 0*ones(size(S,2),1); 
%c_v(25) = 0; %is this necessary? Who knows 142 and 25 (net flux or
%biomass) to 1 in regular loop leads to different total values 
c_y = ones(length(edgeIndices),1);
c = [c_v;c_y;];

%Creating B; b probably needs to be remade in binary search loop
b_s = zeros(size(S,1),1);
lb_indices = lb(edgeIndices); 
ub_indices = ub(edgeIndices);
kComparison = 0;
b = [b_s;lb_indices;ub_indices;kComparison;];

% Creating csense
csense1(1:size(S,1)) = 'E';
csense2(1:length(edgeIndices)) = 'G';
csense3(1:length(edgeIndices)) = 'L';
csense4(1) = 'G'; % sum(ExpressionMinusKVector *yi) => 0 (kComparison)
csense = [csense1 csense2 csense3 csense4  ];

% create MILPproblem
MILPproblem.A = A; % A Matrix
MILPproblem.b = b; % Constraint vector
MILPproblem.c = c; % Objective coeficients (v and Y)
MILPproblem.lb = lb; % Lower bounds 
MILPproblem.ub = ub; % Upper bounds
MILPproblem.csense = csense; % Equality/lower/higher 
MILPproblem.vartype = vartype; % Vartype ('C', 'B', 'I')
MILPproblem.osense = 1; %(-1 max, +1 min)
MILPproblem.x0 = []; % Initial guess value

end