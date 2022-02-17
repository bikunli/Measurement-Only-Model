function [h] = heightfunc(Tableau,option)
% input an arbitary tableau
% output a vector h, representing h(x)
% h(1) = h(x=1)
% This function uses the method introduced  in paper: New J. Phys. 7, 170 (2005)
% Version 4.0, Date: 06/09/2021

[n,n2] = size(Tableau);
if n2 ~= 2*n
    error('The size of Tableau_photon is improper!');
end


h = zeros(1,n);

if nargin == 1 % general case
    rv_ind = n:-1:1;
    Tableau_reversed = Tableau(:,[rv_ind,n+rv_ind]);
    [T_echelon,~] = ECHELON_TABLEAU(Tableau_reversed);
    B = Tableau2Bigram(T_echelon);
    
    for i = 1:n
        h(i) = i - numel( find(B(:,1)>= (n-i+1)) );
    end
elseif nargin == 2
    switch option
        case 'pure' % if pure state
            [T_echelon,~] = ECHELON_TABLEAU(Tableau);
            B = Tableau2Bigram(T_echelon);
            for i = 1:n
                h(i) = n - i - numel( find(B(:,1)> i) );
            end
        otherwise
            h = heightfunc(Tableau);
    end
end
end

function [Bigram] = Tableau2Bigram(Tableau)
% Any tableau has a Bigram, which is the output of this function
% Bigram: there are n rows, each row is [l(g),r(g)]
n = numel(find(~all(Tableau == 0,2))); % non trivial row of tableau
Bigram = zeros(n,2); %
for i_r = 1:n
    Bigram(i_r,1) = find(Tableau(i_r,1:n)+Tableau(i_r,n+(1:n))~=0, 1, 'first');
    Bigram(i_r,2) = find(Tableau(i_r,1:n)+Tableau(i_r,n+(1:n))~=0, 1, 'last');
end
end