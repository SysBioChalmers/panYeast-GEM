function [reducedModel,resultfile] = SpecificModel(strain)
% This function is to generate strain specific models from a panmodel.
% Input of this function are: a panmodel, which will be loaded using the
% function as loadYeastModel and a gene_strainMatrix, which will be loaded
% from ../ComplementaryData/genesMatrix_PresenceAbsence_new.xlsx.
%Usage: [reducedModel,resultfile] = SpecificModel(strain)
%       of [reducedModel,resultfile] = SpecificModel  %note this one
%       generate SSmodels for 1011 strians
%load model
cd ..
model = loadYeastModel;
%change model to raven format
modelr = ravenCobraWrapper(model);
model = modelr;

%load presenceAvsence data
fid2 = fopen('../ComplementaryData/SpecificModelData/genesMatrix_PresenceAbsence_new.csv');
format = repmat('%s ',1,1012);
format = strtrim(format);
data = textscan(fid2,format,'Delimiter',',','HeaderLines',0);
for i = 1:length(data)
genesMatrix(:,i) = data{i};
end
StrianData.genes = genesMatrix(2:end,1);
StrianData.strains = genesMatrix(1,2:end)';
StrianData.levels = cellfun(@str2double,genesMatrix(2:end,2:end));

if nargin<1
    strain = StrianData.strains;
end
%If the supplied object is a character array, then convert it to a cell
%array
if ischar(strain)
    strain={strain};
end

%Check that the strain exists
if ~ismember(upper(strain),upper(StrianData.strains))
    EM='The strain name does not match';
    dispEM(EM);
end


resultfile = [];

% %generate the genelist that don't exist
for j = 1:length(strain)
[~,ID] = ismember(strain(j),StrianData.strains);
lvl = StrianData.levels(:,ID);
lvl_tmp = lvl == 0;
idx = ismember(upper(StrianData.genes),upper(model.genes));
genelist = StrianData.genes(lvl_tmp & idx);



%generate the specific model according to type

reducedModel = removeGenes(model,genelist,true,true,true);


%add annotaion and generate result file
reducedModel.id=[strain{j},'specific model genereted from panYeast'];
cd ../ModelFiles/SSmodels/
save([strain{j},'.mat'],'reducedModel')
cd ../../ComplementaryScripts/
resultfile = [resultfile;strain(j),length(genelist),length(reducedModel.genes),length(reducedModel.rxns),length(reducedModel.mets)];
end

%cd ..
fid2 = fopen('../ComplementaryData/Results/specificModelResultFile.tsv','w');
formatSpec = '%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fid2,formatSpec,'strain','genelistnumber','genes','rxns','mets','sol.x');
for i = 1:length(resultfile(:,1))
    fprintf(fid2,formatSpec,char(resultfile(i,1)),num2str(resultfile{i,2}),num2str(resultfile{i,3}),num2str(resultfile{i,4}),num2str(resultfile{i,5}));
end
fclose(fid2);
end