function model = changeGenesonlyinGEMtoPanID
%This functio is to change genes only exist in GEM to PanID. There are only
%one geneID kept in PanID for duplicate genes.
%input: a mapping list of genes; a cobra model
%output: new model




cd ..
model = loadYeastModel;
%Load mapping list:
fid = fopen('../ComplementaryData/SpecificModelData/genesonlyinGEM.tsv');
mapplist = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
mapping.rxnIDs  = mapplist{3};
mapping.new_GPR  = mapplist{5};
fclose(fid);

[~,rxnindex] = ismember(mapping.rxnIDs,model.rxns);
for i = 1:length(rxnindex)
    model = changeGeneAssociation(model, model.rxns{rxnindex(i)},mapping.new_GPR{i});
end

% Add gene standard name for new genes
fid = fopen('../ComplementaryData/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end


% Save model:
model = removeUnusedGenes(model);
saveYeastModel(model)

end
