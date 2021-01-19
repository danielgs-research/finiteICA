function [S] = gerainteiros_from_probs(probs,P,K,Nobs)
	% modification by Mateus Marcuzzo in 10/05/2019
	% Inspired in geravetorsinais from Daniel Guerreiro e Silva

    % q is a prime
    % K is the number of sources
    % Nobs = the number of signals to be generated

    % probs = zeros(P,K);

    full_size = P^K;

    assert(size(probs,2) == full_size);    
    
          

    % generating the samples here


    range_of_values = 1:full_size;
    S = discrete_rnd(range_of_values,probs,Nobs,1);


end