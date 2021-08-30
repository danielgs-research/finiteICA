function r = binornd(n,p,varargin)
% r = binornd(n,p)
% r = binornd(n,p,sz1)
% r = binornd(n,p,sz1,sz2,...)
% Generates random numbers from the binomial distribution specified by 
% the number of trials n and the probability of success for each trial p.
% sz1,sz2: size of the output matrix
% 
% Binlin Wu, wub1@southernct.edu
% 04/02/2021
%
nn = [varargin{:}];
r = zeros([nn,1]);
x = rand([n nn 1]);
for i = 1:numel(r(:))
    r(i) = nnz(x(:,i) <= p);
end