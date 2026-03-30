function TaskModel = adjustTaskModel(TaskModelOld, diffInReactions)

amountOfPlaceholderReactions = diffInReactions;

% Close all exchange reactions
%model = closeExchangeReactions(model);

tModel = TaskModelOld;
for i = 1:amountOfPlaceholderReactions
    INPUT = sprintf('FakeRxn%i', i);
    [tModel]=addReaction(tModel,['PlaceHolderEmpty',INPUT],[' <=> ',INPUT],[],[],0,0);
end

listReactions = tModel.mets((end-amountOfPlaceholderReactions)+1:end);
[tModel, rxnRemoveList] = removeMetabolites(tModel, listReactions, false);

TaskModel = tModel;

end