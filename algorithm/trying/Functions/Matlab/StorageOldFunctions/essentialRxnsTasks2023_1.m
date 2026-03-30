function [essentialRxns,usedRxns] = essentialRxnsTasks2023_1(model,saveUsed,relaxation)
%Fixed a reinitialization problem (Eva)
%Fixed the minimize/maximize error and added the usedRxns system (Bastien)
%
%Calculates the minimal number of reactions needed to be active for each
% task using Parsimoneous enzyme usage Flux Balance Analysis
%
% USAGE:
%
%    [essentialRxns,usedRxns] = essentialRxnsTasks(mode,saveUsed,relaxation)
%
% INPUT:
%    model:            a model structure
%OPTIONAL INPUT:
%    saveUsed:       Boolean, should all the reactions used for the task be
%    save. Defaults to false.
%
% OUTPUT:
%    essentialRxns:    cell array with the names of the essential
%                      reactions for the task
%   usedRxns:           cell array with the names of the
%                      reactions used for the task
%  
%
% NOTE:
%
%    Essential reactions are those which, when constrained to 0, result in an
%    infeasible problem.
%
% .. Authors:
%   		- Originally written for RAVEN toolbox by Rasmus Agren, 2013-11-17
%   		- Adapted for cobratoolbox and modified to rely on pFBA by Richelle Anne, 2017-05-18
if ~exist('saveUsed','var') || nargin == 1
    saveUsed = false;
end
if nargin < 3 
    relaxation = 1e-06;
end
	[solMin, modelIrrevFMOri]= minimizeModelFlux(model,'min'); %Compute the minimal set of reactions
    if relaxation == Inf
        modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',Inf,'u');
        modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',1e-06,'l');%Will prevent negative or null netFlux
    else
        modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',solMin.f + relaxation,'u');
        modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',max(solMin.f - relaxation, 1e-06),'l');%Will prevent negative or null netFlux
    end
    %modelIrrevFMOri = changeRxnBounds(modelIrrevFMOri,'netFlux',1e-06,'l');
    %Define the list of reactions to test
	rxnsToCheck=modelIrrevFMOri.rxns(abs(solMin.x)>10^-6);
    %[transRxns] = findTransRxns(model,1);
    %rxnsToCheck=rxnsToCheck_raw(~ismember(rxnsToCheck_raw,transRxns));

    % Loop that set to 0 each reaction to test and check if the problem
    % still has a solution
	%essentialRxns={};
    isEssential = false(size(rxnsToCheck));
    if saveUsed
        usedRxns = false(size(modelIrrevFMOri.rxns));
    end
    for i=1:numel(rxnsToCheck)
        %disp(num2str(i))
        modelIrrevFM = modelIrrevFMOri;
        modelIrrevFM.lb(findRxnIDs(modelIrrevFM,rxnsToCheck(i)))=0;
        modelIrrevFM.ub(findRxnIDs(modelIrrevFM,rxnsToCheck(i)))=0;
        modelIrrevFM.csense(1:length(modelIrrevFM.mets),1) = 'E';
        modelIrrevFM.osense = -1;
        modelIrrevFM.A=modelIrrevFM.S;
        sol=solveCobraLP(modelIrrevFM);
        if sol.stat==0 || isempty(sol.full)
            isEssential(i) = true;
            
        elseif saveUsed
            usedRxns(sol.full > 1e-06) = true;
        end
        
    end
    essentialRxns = rxnsToCheck(isEssential);
    rxns_kept=unique(essentialRxns);
    rxns_final=cell(size(rxns_kept));
    if saveUsed
        usedRxns = modelIrrevFMOri.rxns(usedRxns);
        usedRxns = [usedRxns;essentialRxns;rxnsToCheck];
        usedRxns = unique(usedRxns);
    end
    %% Analysis part
    for i=1: length(rxns_kept)
        string=rxns_kept{i};
        if strcmp('_f', string(end-1:end))==1
            rxns_final{i}= string(1:end-2);
        elseif strcmp('_b', string(end-1:end))==1
            rxns_final{i}= string(1:end-2);
        elseif strcmp('_r', string(end-1:end))==1
            rxns_final{i}= string(1:end-2);
        else
            rxns_final{i}=string;
        end
    end
    essentialRxns=unique(rxns_final);
    essentialRxns(strcmp(essentialRxns, 'netFlux')) = [];
    isTemp = contains(essentialRxns,'temporary');
    if sum(isTemp) ~=length(isTemp)
        essentialRxns = essentialRxns(~isTemp);
    end
    
    %% Analysis part2
    if saveUsed  
        usedFinal = cell(size(usedRxns));
        for i=1: length(usedRxns)
            string=usedRxns{i};
            if strcmp('_f', string(end-1:end))==1
                usedFinal{i}= string(1:end-2);
            elseif strcmp('_b', string(end-1:end))==1
                usedFinal{i}= string(1:end-2);
            elseif strcmp('_r', string(end-1:end))==1
                usedFinal{i}= string(1:end-2);
            else
                usedFinal{i}=string;
            end
        end
        usedRxns=unique(usedFinal);
        usedRxns(strcmp(usedRxns, 'netFlux')) = [];
        isTemp = contains(usedRxns,'temporary');
        if sum(isTemp) ~=length(isTemp)
            usedRxns = usedRxns(~isTemp);
        end
    end
    
    
end