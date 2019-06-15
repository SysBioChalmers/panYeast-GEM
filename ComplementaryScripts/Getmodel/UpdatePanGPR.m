%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpdatePanGPR
% Add changes from the Pangenome new anootation for new genes + manual curation on those changes
% Input: model, PanGenes.tsv,SGDgeneNames.tsv.
% As for the reference of new GPR, please find detailed information in:
% ComplementaryData/databases/Pangenes.tsv
% NOTE: changeGeneAssociation.m is a function from cobra
%
% Feiran Li 2018.09.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')

%Load yeast model:
cd yeast-GEM
model    = load('ModelFiles/mat/yeastGEM.mat');
model    = model.model;
yeastVer = model.modelID(strfind(model.modelID,'_v')+1:end);
cd ..

% Load ortholog information
fid           = fopen('../../ComplementaryData/SpecificModelData/PanGenes.tsv');
orth     = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
Pan.genes     = orth{1};
Pan.orthlogs = orth{2};
fclose(fid);

%replace the orthologs with the genes that existed in the model to generate
%new GPRs, isoenzymes in the Pan genes for the existed genes in the model
NEWGPRList = AddOrthologGPRrules(model,Pan.genes,Pan.orthlogs);


for i = 1:length(NEWGPRList(:,1))
    model    = changeGeneAssociation(model, model.rxns{NEWGPRList{i,1}}, ['( ' NEWGPRList{i,3} ' )']);
end
% Delete unused genes (if any)
model = removeUnusedGenes(model);


% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end

% add gene standard name for new genes
fid = fopen('../../ComplementaryData/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

%Remove the cloned repos:
rmdir('yeast-GEM', 's')

% Save model:
model = rmfield(model,'grRules');
saveYeastModel(model)

