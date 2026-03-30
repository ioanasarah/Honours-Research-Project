%%%%%%% Check the mass and charge balance of a reaction

initCobraToolbox(false) 
% only run this section once, use the 'run section' button for the next parts of the script
% rather than 'run', otherwise cobra reinitialises every time
%% Load Task Model
% The uploaded task model does not necessarily have to match the task the reaction is being tested for

model = readCbModel('newModel_Intermediate_SolutionTask250_Sample1_Scored5.000_average4.008_alpha0.500_beta0.000_2024-04-29___00-32-35.mat')

%% Add new reaction to be tested to the model
% Chancge the reaction to the desired formula. In case of error, check that
% all spaces are included between IDs and + or ->

NewModel = addReaction(model, 'NET0001', 'reactionFormula', 'MAM01975[m] + MAM01370[c] + 2 MAM01371[m] + MAM01371[c] + 2 MAM02040[c] + MAM01596[m] ->  MAM03121[c] + MAM01862[c] + 2 MAM01285[m] + MAM01334[c] + 2 MAM02751[m] + MAM02759[c] + 2 MAM02039[m] + 2 MAM02039[c]')

% NewModel.rxns should now have one more cell than Model.rxns containing NET0001

%% Mass Charge balance Analysis
% check if the reaction is actually mass- and charge-balanced
% I. e., it should not create or destroy any atoms or charges

% Define the reaction ID to check
reactionID = 'NET0001';
% Find the index of the reaction in the model
rxnIndex = find(strcmp(NewModel.rxns, reactionID));

[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(NewModel);

imBalancedMass(68) % change the number to the cell of NewModel.rxns that contains the new reaction
% positive imbalances --> production from nothing
% negative imbalances --> consumption without product output

imBalancedCharge(68) % change the number to the cell of NewModel.rxns that contains the new reaction