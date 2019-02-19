function produced = checkProduce(model,mets)
%   Checks which metabolites that can be produced from a model using the
%   specified constraints. This is a less advanced but faster version of
%   checkProduction.
%
%   model       a model structure
%   mets        either a cell array of metabolite IDs, a logical vector 
%               with the same number of elements as metabolites in the model,
%               or a vector of indexes to check for (opt, default model.mets)
%
%   produced    vector with true if the corresponding metabolite could be
%               produced
%
%   Usage: produced=checkProduce(model,mets)
%
%   Adapted from canProduce in Raven toolbox
%    Feiran 2018.09.03
%

if nargin<2
    mets=model.mets;
end
produced= [];
for i = 1:length(mets)
    newModel = model;
    [newModel, rxnIDexists] = addReaction(newModel,['EX_',mets{i}],'metaboliteList',mets(i),...
        'stoichCoeffList',-1,...
        'reversible',0,...
        'checkDuplicate',1);
    if isempty(rxnIDexists)
        newModel.c(1:end) = 0
        newModel.c(end) = 1
    else
        newModel.c(1:end) = 0
        newModel.c(rxnIDexists) = 1
    end
    sol = optimizeCbModel(newModel);
    if sol.obj > 10^-5
        produced=[produced;mets(i),sol.x(rxnIDexists)];
    else
        produced=[produced;mets(i),0];
    end
end
    
     
    

