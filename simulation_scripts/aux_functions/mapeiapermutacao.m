function [Y] = mapeiapermutacao(S,X, permut,q,K)
    % S são as amostras originais
    % X são as observações, após serem misturadas
    % permut é um array (neste caso coluna (NO QICA está sendo linha..)) que indica
    % como se permutaram as posições
    % q é o primo de (q^K)
    % K é a quantidade de fontes.

    % The below line creates an array like this:
    % [q^(K-1) q^(K-2) ... q^2 q^1 1]
    % example for q=2 and K = 5
    % 16 8 4 2 1
    base = q.^(K-1:-1:0);



    Y = permut(base*X+1)-1;

    % 17/04/2019
    % pq 48? 48 is the '0', character zero in ASCII

    %THis is the original one
    % Y = double(dec2bin(Y,K)-48)';

    %this is my modification:
    % 25/04/2019
    Y = double(dec2base(Y,q,K)-48)';
    Nobs = size(Y,2);
end
