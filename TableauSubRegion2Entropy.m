function [Entropy] = TableauSubRegion2Entropy(Tableau, SubRegionIndex)
% Compute the Renyi entropy of a subsystem (SubRegion)

% SubRegionIndex is a subset of {1,2,3,...,L}.
% Tableau is the full Tableau, with its non-empty rows commuting.
% Tableau should have size L x 2L

L = size(Tableau,2)/2;
CompRegionIndex = setdiff(1:L,SubRegionIndex); % The complementary Region
SubLc = numel(CompRegionIndex); % The length of the complementary region
perm_ind = [CompRegionIndex,SubRegionIndex]; % new arrangement of sites
Tableau_perm = Tableau(:,[perm_ind,L+perm_ind]);

% extract the Renyi entropy from the height function at x = SubLc:
h = heightfunc(Tableau_perm,'pure'); 
switch SubLc
    case 0
        Entropy = 0 ;
    otherwise
        Entropy = h(1,SubLc);
end
end


