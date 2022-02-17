function [newrow_list] = tableau_newrow_FactEns(ensemble,num)
% generate a new row for tableau
% The new row is taken from a "Factorizable ensemble",
% with parameters (r,q_x,q_y,q_z)

L = ensemble.L;
if strcmp(ensemble.bc,'open') == 1
    RandomPosition = randi([1,L-ensemble.r+1],num,1);
else
    RandomPosition = randi([1,L],num,1);
end
RandomPaulistring = randsrc(num,ensemble.r,[1, 2, 3; ensemble.qx, ensemble.qy, ensemble.qz]); % 1,2,3 for X,Y,Z
pauli_to_Z2Z2 = [1, 1, 0;
    0, 1, 1]; % each colomn is the corresponding representation in Z2xZ2 (~F_4)
% the first row for the left segment, the second row for the right segment (of tableau)

newrow_list = zeros(num,2*L);
for i = 1:num
    tableau.pre_elements = pauli_to_Z2Z2(:,RandomPaulistring(i,:));
    % initialization of outcome
    tableau.left = zeros( 1, L );
    tableau.right = zeros( 1, L );
    % evaluation
    tableau.left(1:ensemble.r) = tableau.pre_elements(1,:);
    tableau.right(1:ensemble.r) = tableau.pre_elements(2,:);
    % random shifting the measure, sending it to the designated RandomPosition
    tableau.left = circshift(tableau.left,[0,RandomPosition(i) - 1]);
    tableau.right = circshift(tableau.right,[0,RandomPosition(i) - 1]);
    newrow_list(i,:) = cat(2, tableau.left,tableau.right); % glue left and right segments
end
end