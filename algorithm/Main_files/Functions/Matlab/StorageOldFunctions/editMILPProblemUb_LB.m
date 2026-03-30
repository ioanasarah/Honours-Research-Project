function MILPproblem = editMILPProblemUb_LB(MILPproblem, model, epsilon)

% Indice the reactions of the model
edgeIndices = zeros(length(model.rxns),1);
for i = 1:length(model.rxns)
    edgeIndices(i) = i;
end

% Get S
S = model.S;
lb = model.lb;
ub = model.ub;

% Create A size vars
reactionLength = length(edgeIndices);
M = 1;

for i = 1:length(edgeIndices)
    %V UB LB
    MILPproblem.A(i+size(S,1),i+size(S,2)) = lb(edgeIndices(i)) - epsilon;
    MILPproblem.A(i+size(S,1)+length(edgeIndices),reactionLength+i+size(S,2)) = ub(edgeIndices(i)) + epsilon;
    
    %Yprime UB LB
    MILPproblem.A(i+size(S,1)+12*length(edgeIndices),i+ 3*size(S,2)) = lb(edgeIndices(i));
    MILPproblem.A(i+size(S,1)+13*length(edgeIndices),i+ 3*size(S,2)) = ub(edgeIndices(i));
end


%%%% CHECK IF THIS WORKS
lengthUBIndices = 1:length(model.ub);
MILPproblem.ub(lengthUBIndices) = ub;
MILPproblem.lb(lengthUBIndices) = lb;

% b_s = zeros(size(S,1),1); % SV = 0 equality
lb_indices = lb(edgeIndices); % >= Y_LB
ub_indices = ub(edgeIndices); % <= Y_UB
lb_ub_together = [lb_indices;ub_indices];
MILPproblem.b(size(S,1)+1:size(S,1)+2*length(edgeIndices)) = lb_ub_together;


end