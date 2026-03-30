function createExcelOfReactions(i,j,withMilp,LPTable,MinimumReactionsTable,UsedReactionsTable,MILPCooperTable,MILPCooperIrrevTable,MILPKApproxTable)

nanString = ' ';
if nargin < 9
    MILPKApproxTable = repmat({nanString}, 1, 6);
end


if withMilp 
    withMilpstr = "withMILP";
else
    withMilpstr = "onlyLP";
end

numActiveReactions = max([size(MILPKApproxTable, 1), size(LPTable, 1), size(MinimumReactionsTable, 1),size(UsedReactionsTable, 1),size(MILPCooperTable, 1),size(MILPCooperIrrevTable, 1)]);

if size(LPTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(LPTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(LPTable, 2));
    LPTable = vertcat(LPTable, tempArray);
end

if size(MILPKApproxTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(MILPKApproxTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(MILPKApproxTable, 2));
    MILPKApproxTable = vertcat(MILPKApproxTable, tempArray);
end

if size(MinimumReactionsTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(MinimumReactionsTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(MinimumReactionsTable, 2));
    MinimumReactionsTable = vertcat(MinimumReactionsTable, tempArray);
end

if size(UsedReactionsTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(UsedReactionsTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(UsedReactionsTable, 2));
    UsedReactionsTable = vertcat(UsedReactionsTable, tempArray);
end

if size(MILPCooperTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(MILPCooperTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(MILPCooperTable, 2));
    MILPCooperTable = vertcat(MILPCooperTable, tempArray);
end

if size(MILPCooperIrrevTable, 1) < numActiveReactions
    paddingRows = numActiveReactions - size(MILPCooperIrrevTable, 1);
    tempArray = repmat({nanString}, paddingRows, size(MILPCooperIrrevTable, 2));
    MILPCooperIrrevTable = vertcat(MILPCooperIrrevTable, tempArray);
end

disp('hey')






% Extract the reaction names from ActiveReactionFluxesUsedReactions
usedReactions = UsedReactionsTable(:,2);

% Find the reactions only in ActiveReactionFluxesUsedReactions
OnlyInUsedReactions = usedReactions(~ismember(usedReactions, LPTable(:, 4)) & ~ismember(usedReactions, MILPKApproxTable(:, 4)));

% Find the reactions only in ActiveReactionFluxesLP compared to ActiveReactionFluxesUsedReactions
OnlyInLPVersusUsed = LPTable(~ismember(LPTable(:, 4), usedReactions), 4);

% Find the reactions only in ActiveReactionFluxesMILP compared to ActiveReactionFluxesUsedReactions
OnlyInMILPVersusUsed = MILPKApproxTable(~ismember(MILPKApproxTable(:, 4), usedReactions), 4);

% Find the reactions missing in ActiveReactionFluxesMILP compared to ActiveReactionFluxesUsedReactions
MissingInMILPVersusUsed = usedReactions(~ismember(usedReactions, MILPKApproxTable(:, 4)));

% Find the reactions missing in ActiveReactionFluxesLP compared to ActiveReactionFluxesUsedReactions
MissingInLPVersusUsed = usedReactions(~ismember(usedReactions, LPTable(:, 4)));

% Find the reactions only in ActiveReactionFluxesMILP compared to ActiveReactionFluxesLP
OnlyInMILPVsLP = MILPKApproxTable(~ismember(MILPKApproxTable(:, 4), LPTable(:, 4)), 4);

% Find the reactions missing in ActiveReactionFluxesMILP compared to ActiveReactionFluxesLP
MissingInMILPVsLP = LPTable(~ismember(LPTable(:, 4), MILPKApproxTable(:, 4)), 4);

 
padSize = numActiveReactions - numel(OnlyInUsedReactions);
OnlyInUsedReactions = [OnlyInUsedReactions; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(OnlyInLPVersusUsed);
OnlyInLPVersusUsed = [OnlyInLPVersusUsed; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(OnlyInMILPVersusUsed);
OnlyInMILPVersusUsed = [OnlyInMILPVersusUsed; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(MissingInMILPVersusUsed);
MissingInMILPVersusUsed = [MissingInMILPVersusUsed; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(MissingInLPVersusUsed);
MissingInLPVersusUsed = [MissingInLPVersusUsed; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(OnlyInMILPVsLP);
OnlyInMILPVsLP = [OnlyInMILPVsLP; repmat({nanString}, padSize, 1)];

padSize = numActiveReactions - numel(MissingInMILPVsLP);
MissingInMILPVsLP = [MissingInMILPVsLP; repmat({nanString}, padSize, 1)];

MinimumReactionsTable = MinimumReactionsTable(:, [1,2, 4, 6]);
MILPKApproxTable = MILPKApproxTable(:, [2,3, 4, 6]);
LPTable = LPTable(:, [1,2, 4, 6]);

MILPCooperIrrevTable = MILPCooperIrrevTable(:, [2,3, 4, 6]);
MILPCooperTable = MILPCooperTable(:, [2,3, 4, 6]);


columnWidths = [20, 20, 25, 20, 20, 25, 20, 20, 25, 20];
combinedTable = table( MinimumReactionsTable,UsedReactionsTable,LPTable, OnlyInUsedReactions, OnlyInLPVersusUsed,MissingInLPVersusUsed, MILPKApproxTable,MILPCooperTable,MILPCooperIrrevTable, OnlyInMILPVersusUsed, MissingInMILPVersusUsed, OnlyInMILPVsLP, MissingInMILPVsLP);
   
dispName = sprintf('ActiveReactionsTask%iSample%i%s.xlsx', i, j, withMilpstr);
% Write the table to an Excel file
writetable(combinedTable, dispName, 'WriteVariableNames', true, 'AutoFitWidth', false);


end