function pars = listPartitions(L, N)
%LISTPARTITIONS returns an array of all possible non-zero partitions of
%a given number of elements.
%   [PARS] = LISTPARTITIONS(L, N) returns an array of all non-zero
%   N-partitions of L elements. In other words, list all possible {x_i},
%   x_1 + x_2 + ... + x_N = L and x_i >= 1.

if N == 1
    pars = [L];
else
    pars = zeros(0, N);
    for x1 = 1 : L - N + 1
        pars1 = listPartitions(L - x1, N - 1);
        pars = [pars; x1 * ones(size(pars1,1),1), pars1];
    end
end

% while
%     segment1(incState) = segment1(incState) + 1;
%     if sum(segment1) < L - 1, 
%         segment1(incState
%     else
%         incState = incState - 1; 
%     end
% end
