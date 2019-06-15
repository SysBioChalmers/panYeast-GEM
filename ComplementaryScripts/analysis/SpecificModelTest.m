function [results] = SpecificModelTest(strain)
% This function wil run a simple test for SSmodels
% This funciton is to check whether several critical amino acids and
% byproducts can be produced or not
% Input should be SSmodels, be sure to run SpecificModel.m before run this
% one.
% Usage: [results] = SpecificModelTest(strain)
%      or [results] = SpecificModelTest  %note this one will give the
%      result for 1011 SSmodels

if nargin<1
%load presenceAvsence data
genesMatrix = readtable('../ComplementaryData/SpecificModelData/genesMatrix_PresenceAbsence_new.xlsx');
StrianData.genes = genesMatrix.geneID;
StrianData.strains = genesMatrix.Properties.VariableNames(2:end)';
StrianData.levels = table2array(genesMatrix(:,2:end));
strain = StrianData.strains;
end

Results = [];

for i = 1 : length(strain)
        filename = [strain{i},'.mat'];
        cd ../../ModelFiles/SSmodels/
        load(filename);
        cd ../../ComplementaryScripts/
        model = reducedModel;
        %change from raven format to cobra format
        %model = ravenCobraWrapper(reducedModel);
        %for cobra model
        %mets = {'s_1459[e]','s_0067[e]','s_0726[e]','s_0681[e]','s_0933[e]','s_1267[e]','s_0956[e]','s_0966[e]','s_0970[e]','s_0974[e]','s_0982[e]','s_1000[e]','s_0992[e]','s_1004[e]','s_1007[e]','s_1017[e]','s_1022[e]','s_1026[e]','s_1030[e]','s_1033[e]','s_1036[e]','s_1041[e]','s_1046[e]','s_1049[e]','s_1052[e]','s_1057[e]'};
        %for raven model
        mets = {'s_1459','s_0067','s_0726','s_0681','s_0933','s_1267','s_0956','s_0966','s_0970','s_0974','s_0982','s_1000','s_0992','s_1004','s_1007','s_1017','s_1022','s_1026','s_1030','s_1033','s_1036','s_1041','s_1046','s_1049','s_1052','s_1057'};
        % metNames
        %mets = {'succinate [extracellular]','(S)-malate [extracellular]','fumarate [extracellular]','ethanol [extracellular]','isobutanol [extracellular]','ornithine [extracellular]','L-alanine [extracellular]','L-arginine [extracellular]','L-asparagine [extracellular]','L-aspartate [extracellular]','L-cysteine [extracellular]','L-glutamine [extracellular]','L-glutamate [extracellular]','L-glycine [extracellular]','L-histidine [extracellular]','L-isoleucine [extracellular]','L-leucine [extracellular]','L-lysine [extracellular]','L-methionine [extracellular]','L-phenylalanine [extracellular]','L-proline [extracellular]','L-serine [extracellular]','L-threonine [extracellular]','L-tryptophan [extracellular]','L-tyrosine [extracellular]','L-valine [extracellular]'}
        produced = checkProduce(model,mets);
        producedFlux = transpose(produced(:,2));

        %%%change media
        model = reducedModel;
        %model = minimal_Y6(model)
        sol = optimizeCbModel(model);
        if sol.f > 0
        solmin = sol.f;
        else
        solmin = 0;
        end
        [model,pos] = changeMedia_Yeast8(model,'D-glucose exchange','YEP');
        sol = optimizeCbModel(model);
        if sol.f ~= 0
            Results = [Results;filename,producedFlux,solmin,sol.f];
        else
            Results = [Results;filename,producedFlux,solmin,0];
        end
        cd ../ComplementaryData/Results
        save('SpecificModelTestResults.mat','Results')
end
end
