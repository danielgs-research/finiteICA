function [S] = gerainteiros_from_probs(probs,P,K,Nobs)
	% modification by Mateus Marcuzzo in 10/05/2019
	% Inspired in geravetorsinais from Daniel Guerreiro e Silva
    % K is the number of sources
    % Nobs = the number of signals to be generated

    full_size = P^K;

    assert(size(probs,2) == full_size);        

    % generating the samples here
%     S = discrete_rnd(range_of_values,probs,Nobs,1); %OCTAVE
%     S = randsample(full_size, Nobs, true, probs); %Statistics toolbox
    x = 1:full_size;
    idx = discretize(rand(1,Nobs),[0,cumsum(probs)]);
    idx(isnan(idx)) = randi(full_size,[1 sum(isnan(idx))]);
    S = x(idx)';

end