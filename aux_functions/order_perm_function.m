function [Y, opt_perm] = order_perm_function(Px,X,P,K)
% function: Short description
%
% Extended description

[~, opt_perm] = sort(Px);

[Y,] = map_permutation(X, opt_perm,P,K);

end  % function
