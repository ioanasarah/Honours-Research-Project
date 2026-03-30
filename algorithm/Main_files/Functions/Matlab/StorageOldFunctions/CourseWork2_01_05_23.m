%%%%
% Clear vars for next run
%clearvars

%reconModelName = 'recon2_model.xml';
%model = readCbModel(reconModelName);
%type(reconModelName)
%mlStruct = parseXML(reconModelName);


disp("________________________________")


blockedReactions = findBlockedReaction(model);




listRxns = {};
model = creategrRulesField(model);
% model.grRules = model.rules;
% for i =1:length(model.rxnisVersionOfobo__46__ecoID)
%     myString = string((model.rxnisVersionOfobo__46__ecoID(i)));
%     extract = extractBetween(myString,strlength(myString)-2,strlength(myString));
%     if extract == '383'
%         disp(i)
%     end
% end
[reactions, list1] = findRxnsFromGenes(model, '383.1');
genes = printGPRForRxns(model, reactions);

