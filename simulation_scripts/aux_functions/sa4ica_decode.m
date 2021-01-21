function [B] = sa4ica_decode(Px,parameters,beta,k)
	%SA4ICA Probabilities Tensor version - Prime fields only!
	%Good for large number of samples - 2^11 or larger!!!
	% global count;

	% histH = [];
	epsilon = 1e-3;

	N = parameters.K;
	q = parameters.P;
	% Nfits = 0;
	%flag_parada = 0;
	%Ntotal_comb = 0;
	B = eye(N,N); %initial solution
	[h, ~] = decode_(eye(N),Px,parameters); %initial entropies
	% if m > 1
	%     B = B - 1;
	% end

	T = 1; %Kirckpatrick

	while T>epsilon

	    for it=1:k
	        %generate a random move
	        i =randi(N);
	        j = randi(N);
	        V = eye(N);
	%         if m > 1
	%             V = V - 1;
	%             c = randsrc(1,1,0:q^m-2);
	%         else
	            c = randi(q-1);
	%         end
	        V(i,j) = c;  %switch Xi by the combination Xi + c.Xj

	%         Xnew = produtomatrizGF(V, X, q, m, field);
	        [hnew, Pxnew] = decode_(V, Px, parameters);
	        hnew = hnew(i);
	%         H = entrp([combined_signal; X(i,:)],q,m);
	%         [hnew, count] = entrp(Xnew(i,:),q,m);
	%         Nfits = Nfits + count;
	        delta_H = hnew - h(i);
	        %update candidate solution
	        if(delta_H<0 || exp(-delta_H/T) > rand())
	%             B = produtomatrizGF(V,B,q,m,field);
	            B = produtomatrizGF(V,B,q,1,[]);
	            Px = Pxnew;
	            h(i) = hnew;
	        end
	    end
	%     Nfits = Nfits + k;
	    T = beta * T;
	%     fprintf(1,'%d %.4f\n', Nit, T);
	end
	% Nfits = count;
end


function [h,PPy] = decode_(W,PPx,parameters)
	% decode function for sa4ica
	%PPx: probabilities tensor
    PPy = PPx;
    r = parameters.r;
    q = parameters.P;
    K = parameters.K;
    lex = parameters.Lex;
    % global r q K lex;
    % K = length(W);
    lg_cte = log2(q); %correction factor to calculate always logP entropies
    % r=q.^(0:K-1);

    %K-D fft for obtaining the characteristic tensor
    fPPyn=fftn(reshape(PPx,q*ones(1,K)));
    fPPy=fPPyn(:);

    %obtain the characteristic vectors of
    %the linear combinations
    qf=ones(q,q^K);
    qf(2,:)=fPPy;
    if q>2
        qf(q,:)=conj(fPPy);
        for m=2:q/2
            mLex=rem(m*lex,q);
            qf(m+1,:)=fPPy(r*mLex+1);
            qf(q+1-m,:)=conj(qf(m+1,:));
        end
    end

    %translate characteristic vectors into probabilities
    %vectors and then into entropies
    ffq=ifft(qf);
    ffq=max(ffq,eps);
    h=-sum(ffq.*log2(ffq+eps),1);
    h = h(W*r'+1)./lg_cte;

    WLex = rem(W*lex, q);
    PPy(r*WLex+1) = PPx;%new prob tensor after W transform.

    % lgP = log(p)./lg_cte;
    % Pzero = (p==0);
    % PlgP = p.*lgP;
    % PlgP(Pzero) = 0;
    % h = sum(PlgP,2);

    % h = -p.*log2(p) - (1-p).*log2(1-p); %binary marginal entropies
    % h = sum(h);

    % v(it) = 1-(h/M);


end
