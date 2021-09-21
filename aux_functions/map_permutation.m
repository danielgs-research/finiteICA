function [Y] = map_permutation(X, permut,q,K)
    % X is the input data
    % permut is a permutation array of indexes from X
    % q is the alphabet size
    % K is the input dimensionality

    % The line below creates an array like this:
    % [q^(K-1) q^(K-2) ... q^2 q^1 1]
    % example for q=2 and K = 5
    % 16 8 4 2 1
    base = single(q.^(K-1:-1:0));
    Y = permut(base*single(X)+1)-1;
    Y = single((dec2base(Y,q,K)-48)');
end
