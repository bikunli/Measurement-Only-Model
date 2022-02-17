function [new_Tableau] = UpdateTableau(old_Tableau, newrow)
% newrow: new tableau row to be imposed as 'measurement'
% old_Tableau: pre-measurement state
% new_Tableau: post-measurement state

old_Tableau = logical(old_Tableau);
% We assume old_Tableau has size Lx2L
L = size(old_Tableau,1);

% find the zero rows number, it could be non-zero if the stablizer state is
% not pure：
old_NonZeroRowsNUM = numel(find(~all(old_Tableau == 0 ,2))); 
old_Rank = gfrank(old_Tableau,2);

if old_Rank ~= old_NonZeroRowsNUM
    warning('The old Tablau has old_Rank ~= old_NonZeroRowsNum');
end

% find out the rows that are anticommutative with newrow
AntiCommuteIndicator = isanticommutative(old_Tableau,newrow);
AC_index = find(AntiCommuteIndicator);
temp_Tableau = old_Tableau;

if old_NonZeroRowsNUM == L && numel(AC_index) == 0
    new_Tableau = old_Tableau;
    
elseif old_NonZeroRowsNUM < L && numel(AC_index) == 0 % Possibility 1
    temp_Tableau(old_NonZeroRowsNUM+1,:) = newrow;
    if old_Rank == gfrank(temp_Tableau,2)
        new_Tableau = old_Tableau;
    elseif old_Rank + 1 == gfrank(temp_Tableau,2)
        new_Tableau = temp_Tableau;
    end
    
elseif numel(AC_index) == 1 %&& old_NonZeroRowsNUM < L   % Possibility 2
    temp_Tableau(AC_index(1),:) = newrow;
    new_Tableau = temp_Tableau;
    
elseif numel(AC_index) > 1 %&& old_NonZeroRowsNUM < L  % Possibility 3
    % first do the 'gauge transformation'
%     temp_Tableau(AC_index(2:end),:) = mod( temp_Tableau(AC_index(2:end),:) + temp_Tableau(AC_index(1:(end-1)),:) ,2);
    temp_Tableau(AC_index(2:end),:) = mod( temp_Tableau(AC_index(2:end),:) + temp_Tableau(AC_index(1),:) ,2);
    temp_Tableau(AC_index(1),:) = newrow;
    new_Tableau = temp_Tableau;
end

% See if improve efficiency
% new_Tableau = sparse(new_Tableau);

end


function [yes1no0] = isanticommutative(row1,row2)
% if two rows of tableau are anti-commutative, 1 for yes, 0 for no
% This function follows the Appendix C of arXiv: 1905.05195v5, but the
% notation of X and Z definition in the tableau is reversed.
% row1 is considered as the "new row" that left multiplies on "previous" row2.
if size(row1,2)~= size(row2,2)
    warning('Two rows should have the same length!')
end
if mod(size(row1,2),2) == 1
    warning('Each input row should have even number of elements!');
end
L = size(row1,2)/2;
P = [zeros(L),eye(L);eye(L),zeros(L)];
yes1no0 = mod(row1*P*(row2.'),2); 
end


function rk = gfrank(a, p)
%GFRANK Compute the rank of a matrix over a Galois field.
% 
%   RK = GFRANK(A) calculates the rank of the matrix A in GF(2).
%
%   RK = GFRANK(A, P) calculates the rank of the matrix A in GF(P).
%
%   Note: This function performs computations in GF(P) where P is prime. To
%   work in GF(2^M), use the RANK function with Galois arrays.
%
%   See also GFLINEQ, GF/RANK.

%   The method used here is similar to the Gaussian elimination. The
%   algorithm has taken advantage of the binary computation, and double
%   sided elimination has been used.

%   Copyright 1996-2017 The MathWorks, Inc.

narginchk(1, 2);

if nargin == 1
    p = 2; 
else % nargin == 2
    validateattributes(p,{'double','single'}, ...
        {'scalar','real','nonnegative','integer','finite'}, ...
        'gfrank', 'the P input');
    if ~isprime(p)
        error(message('comm:gfrank:PNotPrime'));
    end
end

% make matrix A a long matrix.
[n, m] = size(a);
if n < m
    a = a';
end

x = a(:);
if ((max(x) >=p) || (min(x) < 0) || any(any(floor(x)~=x)))
    error(message('comm:gfrank:InvalidElementsInAB'))
end
[n, m] = size(a);
k = 1;
i = 1;
ii = [];
kk = [];

% forward major element selection
while (i <= n) && (k <= m)
    % make the diagonal element into 1, or jump one line.
    while (a(i,k) == 0) && (k <= m)
        ind = find(a(i:n, k) ~= 0);
        if isempty(ind) && (k == m)
            break;
        elseif isempty(ind)
            k = k + 1;
        else
            indx = find(a(i:n, k) == 1);
            if isempty(indx)
               ind_major = ind(1);
            else
               ind_major = indx(1);
            end
            j = i + ind_major - 1;
            tmp = a(i, :);
            a(i,:) = a(j, :);
            a(j, :) = tmp;
        end
    end

    % clear all nonzero elements in the column except the major element.
    if (a(i,k) ~= 0)
     % to make major element into 1
        if (a(i,k) ~= 1)
           a(i,:) = rem(a(i,k)^(p-2) * a(i,:), p);
        end
        ind = find(a(:,k) ~= 0)';
        ii = [ii i];
        kk = [kk k];
        vec = (k:m);
        for j = ind
            if j > i
                % to make sure the column will be zero except the major element.
                a(j, vec) = rem(a(j, vec) + a(i, vec) * (p - a(j, k)), p);
            end
        end
        k = k + 1;
    end
    i = i + 1;
end

rk = find(sum(a')>0,1,'last');
if isempty(rk)
    rk = 0;
end

end