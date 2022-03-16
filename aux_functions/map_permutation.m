function [Y] = map_permutation(X, permut,q,K)
    % X is the input data
    % permut is a permutation array of indexes from X
    % q is the alphabet size
    % K is the input dimensionality

    % The line below creates an array like this:
    % [q^1 q^2 ... q^(K-2) q^(K-1)]
    % example for q=2 and K = 5
    % 1 2 4 8 16
    base = single(q.^(0:K-1)); 
    Y = permut(base*single(X)+1); 
    Y = int_to_tuple(Y,q,K);
end
