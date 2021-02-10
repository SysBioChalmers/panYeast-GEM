%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeMedia(model,media,flux)
%
% function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources.
%
% INPUT:
%   - model:  An enzyme constrained model
%   - meadia: Media type ('YEP' for complex, 'MAA' minimal with Aminoacids,
%                          'Min' for minimal media)
%   - flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% OUTPUT:
%   - model: The ECmodel with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,pos] = changeMedia_Yeast8(model,c_source,media,flux)

% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange (reversible)'

% first block any uptake
pos           = getUptakeIndexes(model,c_source);
model.ub(pos) = 0;

% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;

%Find substrate production rxn and block it:
pos_rev = strcmpi(model.rxnNames,c_source(1:strfind(c_source,...
                                            ' (reversible)')-1));
model.ub(pos_rev) = 0;
%For growth on fructose and mannose the transport takes place in a passive
%way. [Boles & Hollenberg, 2006]
if strcmp(c_source,'D-fructose exchange')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
elseif strcmp(c_source,'D-mannose exchange')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
end
%The media will define which rxns to fix:
if strcmpi(media,'YEP')
    N = 25;     %Aminoacids + Nucleotides
elseif strcmpi(media,'MAA')
    N = 21;     %Aminoacids
elseif strcmpi(media,'Min')
    N = 1;      %Only the carbon source
end

%UB parameter (manually optimized for glucose on Min+AA):
b = 0.08;
%UB parameter (manually optimized for glucose complex media):
c = 2;


%Define fluxes in case of ec model:
if nargin < 5   %Limited protein    
    if N>1
       flux    = -b*ones(1,N);
       if N>21
           flux(22:25) = -c;
       end
    end
    flux(1) = -1;%set glucose uptake to be 1
end

%Fix values as LBs:
for i = 1:N
    model.lb(pos(i)) = flux(i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUptakeIndexes(model,c_source)
    pos(1)  = find(strcmpi(model.rxnNames,c_source));
    pos(2)  = find(strcmpi(model.rxnNames,'L-alanine exchange'));
    pos(3)  = find(strcmpi(model.rxnNames,'L-arginine exchange'));
    pos(4)  = find(strcmpi(model.rxnNames,'L-asparagine exchange'));
    pos(5)  = find(strcmpi(model.rxnNames,'L-aspartate exchange'));
    pos(6)  = find(strcmpi(model.rxnNames,'L-cysteine exchange'));
    pos(7)  = find(strcmpi(model.rxnNames,'L-glutamine exchange'));
    pos(8)  = find(strcmpi(model.rxnNames,'L-glutamate exchange'));
    pos(9)  = find(strcmpi(model.rxnNames,'L-glycine exchange'));
    pos(10) = find(strcmpi(model.rxnNames,'L-histidine exchange'));
    pos(11) = find(strcmpi(model.rxnNames,'L-isoleucine exchange'));
    pos(12) = find(strcmpi(model.rxnNames,'L-leucine exchange'));
    pos(13) = find(strcmpi(model.rxnNames,'L-lysine exchange'));
    pos(14) = find(strcmpi(model.rxnNames,'L-methionine exchange'));
    pos(15) = find(strcmpi(model.rxnNames,'L-phenylalanine exchange'));
    pos(16) = find(strcmpi(model.rxnNames,'L-proline exchange'));
    pos(17) = find(strcmpi(model.rxnNames,'L-serine exchange'));
    pos(18) = find(strcmpi(model.rxnNames,'L-threonine exchange'));
    pos(19) = find(strcmpi(model.rxnNames,'L-tryptophan exchange'));
    pos(20) = find(strcmpi(model.rxnNames,'L-tyrosine exchange'));
    pos(21) = find(strcmpi(model.rxnNames,'L-valine exchange'));
    pos(22) = find(strcmpi(model.rxnNames,'2''-deoxyadenosine exchange'));
    pos(23) = find(strcmpi(model.rxnNames,'2''-deoxyguanosine exchange'));
    pos(24) = find(strcmpi(model.rxnNames,'thymidine exchange'));
    pos(25) = find(strcmpi(model.rxnNames,'deoxycytidine exchange'));
    pos(26) = find(strcmpi(model.rxnNames,'D-glucose exchange'));
end
