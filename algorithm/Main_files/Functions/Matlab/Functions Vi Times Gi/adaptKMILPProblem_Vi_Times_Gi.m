function MILPproblem = adaptKMILPProblem_Vi_Times_Gi(MILPproblem,taskModel,expressionRxns,ktest,transportRxnIDs,alpha,beta,x0hint,x0ON,TransportReactionPenalty)

[~, ~, transportRxnIDs] = findTransRxns(taskModel);

S = taskModel.S;

betaAdjusted = beta/length(find(expressionRxns ==-1)); % adds up to beta if all transport rxns are INactive
gammaAdjusted = TransportReactionPenalty/length(find(transportRxnIDs ==1 & expressionRxns == -1));
edgeIndices = zeros(length(taskModel.rxns),1);
for i = 1:length(taskModel.rxns)
    edgeIndices(i) = i;
end
% 
% ind = [];
% j = [];
v_flux = [];
v1 = [];
v2 = [];
count_vflux = 0;
count = 0;
count2 = 0;
 

for i = 1:length(edgeIndices)

    count_vflux = count_vflux + 1;
    if(expressionRxns(i) ~= -1)
        v_flux(count_vflux) = expressionRxns(i);
    else
        if sampleMedian ~= 0
            v_flux(count_vflux) = sampleMedian; 
        else
            v_flux(count_vflux) = 0;
        end
    end

    count = count + 1;
    % ind(count) = 1+size(S,1)+7*length(edgeIndices);
    % j(count) = i+ size(S,2);
    if(expressionRxns(i) ~= -1)
        v1(count) = -ktest;
    else
        if sampleMedian ~= 0
            v1(count) = -ktest;
        else
            v1(count) = 0;
        end
    end

    count2 = count2 + 1;
    % ind(count) = 1+size(S,1)+7*length(edgeIndices);
    % j(count) = i+ 2*size(S,2);
    if(expressionRxns(i) == -1)
        if transportRxnIDs(i) == 1 %if transport reaction and OR, give extra incentive to turn off
            v2(count2) = ((1-alpha)*2*betaAdjusted)+gammaAdjusted;
        else
            v2(count2) = ((1-alpha)*2*betaAdjusted); %*2 is there since 0.5 is alpha default value
        end
    else
        v2(count2) = 0;
    end

end

v_flux = [v_flux];
v1 = [v1];
v2 = [v2];


% Replace the desired part of the original matrix with the adjusted part
MILPproblem.A(1+size(S,1)+5*length(edgeIndices), 1+0*size(S,2):length(edgeIndices)+0*size(S,2)) = v_flux;
MILPproblem.A(1+size(S,1)+5*length(edgeIndices), 1+1*size(S,2):length(edgeIndices)+1*size(S,2)) = v1;
MILPproblem.A(1+size(S,1)+5*length(edgeIndices), 1+2*size(S,2):length(edgeIndices)+2*size(S,2)) = v2;

if x0ON
    MILPproblem.x0hint = x0hint;
    MILPproblem.VarHintVal = x0hint;
end
spy(MILPproblem.A)
