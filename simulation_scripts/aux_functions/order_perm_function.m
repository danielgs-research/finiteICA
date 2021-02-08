function [Y, opt_perm] = order_perm_function(Px,X,P,K)
% function: Short description
%
% Extended description

[~, opt_perm] = sort(Px);

[Y,] = mapeiapermutacao([], X, opt_perm,P,K);

end  % function
