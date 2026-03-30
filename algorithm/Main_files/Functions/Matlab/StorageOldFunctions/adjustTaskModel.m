function taskModel = adjustTaskModel(taskModelOld, diffInReactions)

amountOfPlaceholderReactions = diffInReactions;
%CLose all exchange reactions
%model = closeExchangeReactions(model);

tModel = taskModelOld;
for i =1:amountOfPlaceholderReactions
    INPUT = sprintf('FakeRxn%i', i);
    [tModel]=addReaction(tModel,['PlaceHolderEmpty',INPUT],[' <=> ',INPUT],[],[],0,0);
end

listReactions = tModel.mets((end-amountOfPlaceholderReactions)+1:end);
[tModel, rxnRemoveList] = removeMetabolites(tModel, listReactions, false);

taskModel = tModel;

end