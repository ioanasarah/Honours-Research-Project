clearvars
model = readCbModel('e_coli_core.mat');

%spy(model.S);

% %Create Nullmodel
% NullS = null(model.S);
% 
% %Extract rows corresponding to ATP and ADP
% atp_idx = find(strcmp(model.metFormulas, 'C10H12N5O13P3'));
% adp_idx = find(strcmp(model.metFormulas, 'C10H12N5O10P2'));
% amp_idx = find(strcmp(model.metFormulas, 'C10H12N5O7P'));
% atp_adp_rows = model.S([atp_idx, amp_idx adp_idx], :);
% 
% 
% %are rows linearly dependent
% if rank(atp_adp_rows) < size(atp_adp_rows, 1)
%     disp('coupled')
% else
%     disp('notcoupled')
% end
% 
% 
% model_changedBound = changeRxnBounds(model, {'EX_glc__D_e'}, -20); %-20 is flux 
% model_changedBoundKnockout =  changeRxnBounds(model_changedBound, {'EX_glc__D_e'}, 0);
% model_opt = optimizeCbModel(model_changedBound);
% listWithrxnNames = {};
% model_opt.v = round(model_opt.v, 4); %'ICDHyr
% for i = 1:length(model_changedBound.rxns)
%     %disp(model_opt.v(i))
%     if model_opt.v(i) == model_opt.v(29) % 29  = citrate synthase
%     listWithrxnNames(end+1) = model_changedBound.rxns(i);
% 
%     end 
% end
% disp(listWithrxnNames)



%turn off reactions all carbon reactions
CarbonReactions = {'EX_ac_e','EX_etoh_e','EX_fru_e', 'EX_lac__D_e','EX_glc__D_e','EX_pyr_e','EX_succ_e'};
model_changedBoundsNoCarbon = changeRxnBounds(model, CarbonReactions, 0);

%create Results_Table
Results_Table = table(CarbonReactions', ...
    zeros(length(CarbonReactions),1), ...
    zeros(length(CarbonReactions),1), ...
    'VariableNames', {'Names', 'Aerobic', 'Anaerobic'});
Results_Table.Aerobic = string(Results_Table.Aerobic);
Results_Table.Anaerobic = string(Results_Table.Anaerobic);


%loop through
for i = 1:length(CarbonReactions)
    
    % Open O2 reaction
    model_changedBoundsNoCarbon = changeRxnBounds(model_changedBoundsNoCarbon, {'EX_o2_e'}, -1000,'l');
    model_changedBoundsNoCarbon = changeRxnBounds(model_changedBoundsNoCarbon, {'EX_o2_e'}, 1000,'u');
    
    % Open single carbon reaction and optimize model, display biomass
    model_changedBoundsNoCarbon = changeRxnBounds(model_changedBoundsNoCarbon, CarbonReactions(i), -20);
    model_opt = optimizeCbModel(model_changedBoundsNoCarbon);    

    if(model_opt.stat ==1)
        Results_Table.Aerobic(i) = model_opt.v(25)
    else 
        Results_Table.Aerobic(i) = "Infeasible"
    end

    %close O2 reaction
    model_changedBoundsNoCarbon = changeRxnBounds(model_changedBoundsNoCarbon, {'EX_o2_e'}, 0);

    % optimize model, display biomass
    model_opt = optimizeCbModel(model_changedBoundsNoCarbon);
    
    if(model_opt.stat ==1)
        Results_Table.Anaerobic(i) = model_opt.v(25)
    else 
        Results_Table.Anaerobic(i) = "Infeasible"
    end

    % close carbonreaction
    model_changedBoundsNoCarbon = changeRxnBounds(model_changedBoundsNoCarbon, CarbonReactions(i),0);
    
end


