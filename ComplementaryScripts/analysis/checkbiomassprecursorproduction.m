function [mapping] = checkbiomassprecursorproduction
%This function is to check in coremodel, which biomass precursors cannot be synthsised
%Inputs: Coremodel and a model which can produce all biomass precursors


cd ../ModelFiles/SSmodels
load('coremodel.mat')
load('ARN.mat') % this model can be either model which can produce all biomass precursors


Coremodel = Coremodel;
reducedModel = reducedModel;
rxnlist = setdiff(reducedModel.rxns,Coremodel.rxns);
%caculate sol for reduced model and coremodel
sol_core = optimizeCbModel(Coremodel,'max','one');
sol_core.f
sol_red = optimizeCbModel(reducedModel,'max','one');
sol_red.f
bio_rxn = {'r_4048';'r_4049';'r_4050';'r_4598';'r_4599'};% all biomass pseudoreactions except protein, we will manually add amino acid production, since the precursor in protein_pseudoreaction are aa_chargerd tRNAs
Coremodel1 = ravenCobraWrapper(Coremodel);
produced = canProduce(Coremodel1,Coremodel1.mets);
[~,bio_rxn_index] = ismember(bio_rxn,Coremodel.rxns);
mets = [];
for i = 1:length(bio_rxn_index)
    mets_temp = find(Coremodel.S(:,bio_rxn_index(i))< 0 )
    mets = [mets;mets_temp];
end
aa = {'s_1267','s_0956','s_0966','s_0970','s_0974','s_0982','s_1000','s_0992','s_1004','s_1007','s_1017','s_1022','s_1026','s_1030','s_1033','s_1036','s_1041','s_1046','s_1049','s_1052','s_1057'};
[~,aa_index] = ismember(aa,Coremodel1.mets);
aa_index = transpose(aa_index);
mets = [mets;aa_index];
Coremodel.metNames(mets);
result = produced(mets);
metsindex = mets(~result);
mets_query = Coremodel1.mets(metsindex);
metname = Coremodel.metNames(metsindex);
rxnlist = setdiff(reducedModel.rxns,Coremodel.rxns);


%minimize the rxnlist that we should test
model = reducedModel;
removedrxn = [];
for i = 1:length(rxnlist)
    modelOut = removeRxns(model, rxnlist{i});
    sol = optimizeCbModel(modelOut,'max','one');
    %sol.f
    if sol.f > 0
        model = modelOut;
        removedrxn = [removedrxn;rxnlist{i}];
    end
end


% rxnx_exist is the list that make the difference between coremodel and the
% reduced moel
rxnx_exist = setdiff(rxnlist,removedrxn);
mapping = [];
model_out = model;
for i = 1:length(rxnx_exist)
    model = model_out;
    model = removeRxns(model, rxnx_exist{i});
    for j = 1:length(mets_query)
    [model1, rxns]=addExchangeRxns(model,'out',mets_query(j));
    model1 = changeObjective(model1,model1.rxns(end),+1);
    sol = optimizeCbModel(model1,'max','one');
    sol.f
    if sol.f <= 0 
        mapping = [mapping;rxnx_exist(i),mets_query(j),metname(j)];
    else 
        model_out = model; 
    end
    end
end
save('CoremodelBiomassProductionResults.mat','mapping')