function [S] = generate_from_probs(probs,P,K,Nobs)		
    % P is the alphabet size
    % K is the dimension
    % Nobs = the number of observations to be generated

    full_size = P^K;

    assert(size(probs,2) == full_size);        

    % generating the samples here
    x = 1:full_size;
    idx = discretize(rand(1,Nobs),[0,cumsum(probs)]);
    idx(isnan(idx)) = randi(full_size,[1 sum(isnan(idx))]);
    S = x(idx)';

end