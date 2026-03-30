function emptyTaskModel =createEmptyTaskModel(model)

amountOfPlaceholderReactions = 100;
%CLose all exchange reactions
%model = closeExchangeReactions(model);

tModel = model;
for i =1:amountOfPlaceholderReactions
    INPUT = sprintf('FakeRxn%i', i);
    [tModel]=addReaction(tModel,['PlaceHolderEmpty',INPUT],[' <=> ',INPUT],[],[],0,0);
end

listReactions = tModel.mets((end-amountOfPlaceholderReactions)+1:end);
[tModel, rxnRemoveList] = removeMetabolites(tModel, listReactions, false);

emptyTaskModel = tModel;


end