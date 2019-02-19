%generate figures for model rxn numbers and met numbers

% load data
fid           = fopen('../ComplementaryData/Results/specificModelResultFile.tsv');
StrainModels     = textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
SSModel.genes     = StrainModels{3};
SSModel.rxns = StrainModels{4};
SSModel.mets    = StrainModels{5};
SSModel.genelist = StrainModels{2};
fclose(fid);

% Mets
group = categorical(SSModel.mets);
summary(group)


%figure classfication of strains
figure
histogram(group)

ylabel('Models','FontSize',24,'FontName','Arial')
xlabel('Mets in different ssModels','FontSize',24,'FontName','Arial')

set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.11 0.11 0.75 0.8]);


% Rxns
group = categorical(SSModel.rxns);
summary(group)


%figure classfication of strains
figure
histogram(group)

ylabel('Models','FontSize',24,'FontName','Arial')
xlabel('Rxns in different ssModels','FontSize',24,'FontName','Arial')

set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.11 0.11 0.75 0.8]);

% Genes
group = categorical(SSModel.genes);
summary(group)


%figure classfication of strains
figure
histogram(group)

ylabel('Models','FontSize',24,'FontName','Arial')
xlabel('Genes in different ssModels','FontSize',24,'FontName','Arial')

set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.11 0.11 0.75 0.8]);