%Simulated Annealing for ICA over Galois Fields
%Daniel Guerreiro e Silva - 29/5/2018
%%%%INPUT:
%X: N x Nobs mixtures matrix over GF(q^m)
%field: list of GF(q^m) elements (non-prime fields)
%beta: Temperature decay rate
%k: number of iterations at each temperature level

%%%%OUTPUT:
%B: separation matrix
%Nfits: total number of cost function evaluations

function [B,Nfits] = sa4ica(X,q,m,field,beta,k)

epsilon = 1e-3;

N = size(X,1); %number of mixtures
B = eye(N,N); %initial solution
[h, Nfits] = entrp(X,q,m); %initial entropies
if m > 1
    B = B - 1;
end

T = 1; %initial temperature (Kirckpatrick)

while T>epsilon

    for it=1:k
        %generate a random move
        i =randi(N);
        j = randi(N);
        V = eye(N);
        if m > 1 %non-prime field            
            V = V - 1;
            c = randsrc(1,1,0:q^m-2);
        else %prime field
            c = randi(q-1);
        end
        V(i,j) = c;  %switch Xi by the combination Xi + c.Xj
        
        Xnew = product_GFmatrix(V, X, q, m, field); %extract candidate sources
        
        [hnew, count] = entrp(Xnew(i,:),q,m);  %cost function
        Nfits = Nfits + count;
        delta_H = hnew - h(i);
        %update candidate solution
        if(delta_H<0 || exp(-delta_H/T) > rand())                    
            B = product_GFmatrix(V,B,q,m,field); %update solution
            X = Xnew;
            h(i) = hnew;
        end       
    end

    T = beta * T;    
%     fprintf(1,'%.4f\n', T);
end

end
    
% cost function - entropy of a signal Y over GF(q^m)
function [v,count] = entrp(Y,q,m)

P = q^m; %field order

[n, Nobs]= size(Y);

lg_cte = log(P); %correction factor to calculate always logP entropies
v = ((P-1)/(2*Nobs)).*ones(n,1); %Miller-Madow correction

if(m>1)%non-prime field
    Py = histc(Y,-1:P-2,2)./Nobs; %pmf estimation for all q^m symbols
else
    Py = histc(Y,0:P-1,2)./Nobs; %pmf estimation for all q symbols
end
lgPy = log(Py)./lg_cte;
Pzero = (Py==0);
PlgP = Py.*lgPy;
PlgP(Pzero) = 0;
v = v - sum(PlgP,2);

count = n;

end





