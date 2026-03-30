function ActiveReactionFluxes = getActiveReactionsWithFluxesFromModelAndIrrevTROR_inter_save(taskModel, solution,expValueMeanAdjustedTask,score , i, osenseStr , irrev,j, alpha, beta, epsilon)


if irrev 
    irrevStr = "Irrev";
else
    irrevStr = "NoIrrev";
end


if osenseStr == "min"
    osenseStr = "Min";
elseif osenseStr == "max"
    osenseStr = "Max";
end


reactionLength = length(taskModel.rxns);
selectedRows = round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-(epsilon/100));
validValues = expValueMeanAdjustedTask(selectedRows) ~= -1;
selectedExpression = expValueMeanAdjustedTask(selectedRows);

averageScore = sum(selectedExpression(validValues)) / length(selectedExpression(validValues)) ;



currentDateTime = datestr(now, 'yyyy-mm-dd___HH-MM-SS');
dispName =  sprintf("Task%i_Sample%i_Scored%.3f_average%.3f",i,j,score,averageScore);

reactionLength = length(taskModel.rxns);
%disp(reactionLength)
%disp(length(solution.cont))

if ~isstruct(solution)
     ActiveReactionFluxes=0;
else
yiFlagsPlus =  solution.int(1:reactionLength);
% yiFlagsMinus = solution.int(reactionLength+1:2*reactionLength);
% yiFlagsCombined = zeros(reactionLength,1);
% 
% for j = 1:reactionLength
%     if yiFlagsPlus(j) == 1 && yiFlagsMinus(j) == 1
%         fprintf("rxn %s has both forward and backward reactions active\n", j)
%     end
% 
%     if yiFlagsPlus(j) == 1
%         yiFlagsCombined(j) = 1;
%     elseif yiFlagsMinus(j) == 1
%         yiFlagsCombined(j) = -1;
%     end
% end

newVector = zeros(size(solution.cont(1:reactionLength))); 
newVector(round(abs(solution.cont(1:reactionLength)),6)>= (epsilon-(epsilon/100))) = 1;
selectedRows = newVector == 1;

printFlag = false;


tmodel = taskModel;

%tmodel.mets = tmodel.metNames;
rxnFormulas = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

b = printRxnFormula(tmodel,tmodel.rxns(selectedRows),printFlag );

pattern = '(\d+\s)?MAM\d+\[[a-z]\]';
pattern3 = 'MAM\d+\[[a-z]\]';
pattern2 = '\[(.*?)\]';

rxnFormulas = cell(length(b),1);
for bIterator = 1:length(b)
    c = string(b(bIterator));
    fullFormula = "";
    splitString = split(c, '->');
    metabolitesBeforeForIdx = regexp(splitString{1}, pattern3, 'match');
    metabolitesBefore = regexp(splitString{1}, pattern, 'match');
    beforeMetsArray = cell(length(metabolitesBefore),1);
    for cIterator = 1: length(metabolitesBefore)
        idx = findMetIDs(taskModel,metabolitesBeforeForIdx(cIterator));
        metName = taskModel.metNames(idx);
        numberStr = regexprep(metabolitesBefore{cIterator}, '(\d+\s)?MAM\d+\[[a-z]\]', '$1');
        numberOfMetabolites = str2double(numberStr);
        extension = regexp(metabolitesBefore{cIterator}, pattern2, 'tokens');
        if  isnan(numberOfMetabolites)  
             newName = sprintf("%s[%s]", string(metName), string(extension));
        elseif mod(numberOfMetabolites, 1) == 0
            newName = sprintf("%.0f %s[%s]",numberOfMetabolites , string(metName), string(extension));

        else
           newName = sprintf("%.2f %s[%s]",numberOfMetabolites , string(metName), string(extension));
        end
       
        beforeMetsArray(cIterator) = cellstr(newName);
    end

    metabolitesAfter = regexp(splitString{2}, pattern, 'match');
    metabolitesAfterForIDX = regexp(splitString{2}, pattern3, 'match');
    afterMetsArray = cell(length(metabolitesAfter),1);
    for cIterator = 1: length(metabolitesAfter)
        idx = findMetIDs(taskModel,metabolitesAfterForIDX(cIterator));
        metName = taskModel.metNames(idx);
        numberStr = regexprep(metabolitesAfter{cIterator}, '(\d+\s)?MAM\d+\[[a-z]\]', '$1');
        numberOfMetabolites = str2double(numberStr);
        extension = regexp(metabolitesAfter{cIterator}, pattern2, 'tokens');
        if  isnan(numberOfMetabolites)  
             newName = sprintf("%s[%s]", string(metName), string(extension));
        elseif mod(numberOfMetabolites, 1) == 0
            newName = sprintf("%.0f %s[%s]",numberOfMetabolites , string(metName), string(extension));
           
        else
            newName = sprintf("%.2f %s[%s]",numberOfMetabolites , string(metName), string(extension));
        end      
        afterMetsArray(cIterator) = cellstr(newName);
    end

    separator1 = ' + ';
    separator2 = ' -> ';
    beforePart = strjoin(beforeMetsArray, separator1);
    afterPart = strjoin(afterMetsArray, separator1);
    fullFormula = [beforePart, separator2, afterPart];
    fullFormula = cellstr(fullFormula);
    rxnFormulas(bIterator) =fullFormula; 
end



%disp(rxnFormulas(1:10))





newCellArray = [num2cell(yiFlagsPlus(selectedRows)), num2cell(solution.cont(selectedRows)), num2cell(expValueMeanAdjustedTask(selectedRows)), taskModel.rxns(selectedRows), taskModel.rxnNames(selectedRows),rxnFormulas];


removeRows = newVector == 0;

newModel = removeRxns(taskModel, taskModel.rxns(removeRows));

solutionName = string("lValSol_Intermediate_Sol"+dispName+".mat");
newModelName = string("newModel_Intermediate_Solution"+dispName+".mat");
fluxSolutionName = string("fluxTable_Intermediate_Solution"+dispName+".mat");

save(solutionName, "solution");
save(newModelName, "newModel");
save(fluxSolutionName, "newCellArray");

ActiveReactionFluxes = newCellArray;
writecell(ActiveReactionFluxes,string("ActiveReactionsTask"+dispName+".xlsx"),'Sheet','ActiveReactionsTask','Range','A1');

end