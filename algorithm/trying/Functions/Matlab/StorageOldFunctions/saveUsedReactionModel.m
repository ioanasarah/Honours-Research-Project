function saveUsedReactionModel (taskModelForCharnesCooper,i,j)

dispName =  sprintf("UsedReactionsTask%iSample%i",i,j);
newModelName = string("newModel"+dispName+".mat");
save(newModelName, "taskModelForCharnesCooper");
end