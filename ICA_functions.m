% This a script file which compilates many discrete ICA functions in one place
% 31/07/2020
% Also, we did some indentation correction for reading purposes
% The original functions may have many lines of comments,
% These ones here will be more breafly written.

% This is a script file.
1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% FUNCTIONS IN THIS FILE:

% decode_ -- used in sa4ica_decode

% sa4ica_decode -- the sa4ica algorithm itself

% QICA_function -- QICA algorithm, made by Painsky's et al.

% produtomatrizGF -- utils function

% OD_search -- utils used in Painsky's et al functions

% mapeiapermutacao -- util function which does maps the permutation 
%   -- needs a better explanation, can't remember

% mapeiainteiro_to_tuple -- maps an integer form column_wise_matrix to a tuple
%   -- seems to have something wrong here

% is_linear_comb -- check if's is a linear combination from the giving matrix
%   -- not being used explicitly, may need to be incorpored in GLICA_function

% GLICA_function -- GLICA's Painsky's et al algorithm

% geravetorsinais -- utils which generate the signal given parameters

% geramatrizmistura -- generates a mixing matrix

% gerainteiros_from_probs -- generates integers given probabilities

% generate_random_pmf -- generates a random pmf, used in joint_pmf to ICA experiments

% add_to_counter -- utils for generate_pai_P

% generate_pai_P -- ??
%    -- don't remember what it does, but it's related to all possible tuples

% generate_pai -- similar to pai_P, but restricted to base 2

% estimate_marg_probs -- estimates marginal probs from the observations

% entropy_from_frequencies -- calculates entropy given frequencies

% calc_k_params -- utils for QICA_function or BICA_function

% add_to_counter_v2 -- renaming for add_to_counter, used in some of the
% BICA_function or QICA_function, made by Painsky's et al

% calc_k_params_P -- utils for QICA_function or BICA_function

% allVL1nonrecurs -- utils for BICA_function and/or QICA_function

% BICA_function -- the BICA algorithm by Painsky's et al

% assing_slopes -- utils for BICA_function and/or QICA_function

% america -- AMERICA algorithm (Yeredor)


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [opt_est_vals,opt_appox_lin_min2,opt_appox_ent_with_appox_vals2,overall_opt_perm,opt_v_vec]=QICA_function(n,P,p,exhaustive,min_k,max_k,I)
	sorted_p=sort(p);

	[pai]=generate_pai_P(n,P); 
	tot_ent=-sum(p(p>0).*log2(p(p>0)));

	ent_vec=tot_ent*ones(max_k-min_k+1,1);

	lin_min_vec=zeros(max_k-min_k+1,1);
	appox_lin_min_vec=zeros(max_k-min_k+1,1);
	ent_with_lin_params=zeros(max_k-min_k+1,1);
	appox_ent_with_lin_params=zeros(max_k-min_k+1,1);
	est_var=zeros(max_k-min_k+1,1);

	appox_lin_min_vec2=zeros(max_k-min_k+1,1);
	appox_ent_with_lin_params2=zeros(max_k-min_k+1,1);

	num_of_steps_to_find_solution=zeros(max_k-min_k+1,1);
	num_of_initializations_to_find_solution=zeros(max_k-min_k+1,1);

	for k=min_k:max_k
	    [slopes, consts]=calc_k_params_P(k,P);

	    if exhaustive

	        %exhustive algorithm:

	        % I'm using the non_mex form
	        % 25/04/2019
	        % I'm wondering if the k^(P-1) here is correct for just k. Like BICA_function
	        %
	        v = allVL1nonrecurs(k^(P-1), n); %we have n entropies to calculate, each one of them may fall in each of the k^(p-1) cells we defined. 
	        v=fliplr(v);

	        lin_min=inf;
	        opt_v_vec=NaN;

	        for iter=1:size(v,1)
	            broken=0;
	            coef=zeros(P^n,1);
	            calculated_k=0;
	            v_vec=v(iter,:);
	            l=find(v_vec, 1, 'first');
	            while calculated_k<n
	                if (max(slopes(l,:))==inf)
	                    broken=1;
	                    break;
	                end
	                a=assign_slopes(pai(:,calculated_k+1:calculated_k+v_vec(l)),slopes(l,:));
	                coef=coef+sum(a,2);
	                calculated_k=calculated_k+v_vec(l);
	                l=l+1;
	            end
	            if ~broken

	                [sorted_coef ind]=sort(coef,'descend');
	                lin_min_iter=sorted_coef'*sorted_p+v_vec*consts;


	                 if (lin_min_iter<lin_min)
	                    lin_min=lin_min_iter;
	                    if lin_min_iter<tot_ent
	                        problem=1;
	                    end
	                    opt_v_vec=v_vec;
	                 end
	            end
	         end


	        %calc est_vals
	        coef=zeros(P^n,1);
	        calculated_k=0;
	        for l=1:size(v,2)
	            a=assign_slopes(pai(:,calculated_k+1:calculated_k+opt_v_vec(l)),slopes(l,:));
	            coef=coef+sum(a,2);
	            calculated_k=calculated_k+opt_v_vec(l);
	        end


	        %rank the vector p according to the order of the coef vector 
	        vector=coef';
	        [sorted_values order]=sort(vector);
	        rank=zeros(1,length(vector));
	        rank(order) = 1:length(vector);


	        ind=max(rank)-rank+1;
	        opt_perm=sorted_p(ind);
	        est_vals=zeros(P,n);
	        for m1=1:P
	            for m2=1:n
	                est_vals(m1,m2)=sum((pai(:,m2)==(m1-1)).*opt_perm);
	            end
	        end
	        opt_est_vals=est_vals;
	        lin_min_vec(k-min_k+1)=lin_min;
	        ent_with_est_vals=sum(sum(-est_vals.*log2(est_vals)));
	        ent_with_lin_params(k-min_k+1)=ent_with_est_vals;
	        opt_appox_ent_with_appox_vals2=ent_with_est_vals;
	        overall_opt_perm=ind;
	        opt_appox_lin_min2=lin_min;
	    end
	    
	    
	    %%%%% appox. algorithm
	    if ~exhaustive
	        finish=0;
	        counter=1;
	        %opt_appox_lin_min - holds the value of the upper bound piecewise objective
	        %opt_appox_ent_with_appox_vals - holds the value of the true objective, when applying to it the params we found
	        
	        %For debugging purposes, we have two "minimizations" we perform - in the first we keep track 
	        % of both variables, according to the lowset value of
	        % opt_appox_ent_with_appox_vals we have seen so far:
	        opt_appox_lin_min=inf;    
	        opt_appox_ent_with_appox_vals=inf;

	        %In the second we keep track of both variables, according to the 
	        % lowset value of opt_appox_lin_min we have seen so far:
	        opt_appox_lin_min2=inf;
	        opt_appox_ent_with_appox_vals2=inf;

	        opt_iter=inf;
	        %we run the appox alg. enough times and keep the optimal vals we find
	        while ~finish
	            s = sort(randperm(size(consts,1)+n-1,size(consts,1)-1));
	            v_vec = diff([0 s size(consts,1)+n]) - 1;

	            [appox_lin_min appox_ent_with_appox_vals appox_est_vals appox_v_vec iter terminated opt_perm]=OD_search(P,n,p,slopes,consts,pai,v_vec);

	            if terminated==1
	                test_it=1;
	            end
	            if appox_ent_with_appox_vals<opt_appox_ent_with_appox_vals
	                opt_appox_lin_min=appox_lin_min;
	                opt_appox_ent_with_appox_vals=appox_ent_with_appox_vals;
	                opt_iter=iter;
	            end
	            if appox_lin_min<opt_appox_lin_min2
	                opt_appox_lin_min2=appox_lin_min;
	                opt_appox_ent_with_appox_vals2=appox_ent_with_appox_vals;
	                opt_v_vec=appox_v_vec;
	                opt_iter=iter;
	                opt_est_vals=appox_est_vals;
	                overall_opt_perm=opt_perm;
	            end

	            if counter==I-1
	                finish=1;
	            end
	            counter=counter+1;
	        end
	    

	        appox_lin_min_vec(k-min_k+1)=opt_appox_lin_min;
	        appox_ent_with_lin_params(k-min_k+1)=opt_appox_ent_with_appox_vals;

	        appox_lin_min_vec2(k-min_k+1)=opt_appox_lin_min2;
	        appox_ent_with_lin_params2(k-min_k+1)=opt_appox_ent_with_appox_vals2;

	        num_of_steps_to_find_solution(k-min_k+1)=opt_iter;
	        num_of_initializations_to_find_solution(k-min_k+1)=counter;
	    
	    end  
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [C] = produtomatrizGF(A,B,q,m,field) %#codegen
	%auxiliary routine to perform matrix multiplication over GF(q^m)
	%INPUT
	%A, B: input matrices, defined over GF(q^m)
	%field: list of GF(q^m) elements

	%OUTPUT
	%C: product matrix

	%Daniel Guerreiro e Silva - 12/01/2015
    if(m==1)%prime field
        C = rem(A*B,q);
    else %non-prime field
        lines = size(A,1);
        columns = size(B,2);
        C = zeros(lines,columns);   
        for j=1:columns
            for i=1:lines        
                list = gfmul(A(i,:),B(:,j)',field);
                C(i,j) = gfsum(list,field);
            end
        end                        
        C(C==-Inf) = -1;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lin_min ent_with_appox_vals est_vals opt_v_vec iter terminated ind]=OD_search(P,n,p,slopes,consts,pai,v_vec)


    lin_min=inf;
    ent_with_appox_vals=inf;
    est_vals=NaN;
    opt_v_vec=NaN;
    terminated=0;
    
    
    iter=0;
    flag=1;
    sorted_p=sort(p);
    while flag
        if iter>1
           working=1; 
        end
        %extract coef
        l=find(v_vec, 1, 'first');
        calculated_k=0;
        coef=zeros(P^n,1);
        while calculated_k<n
            a=assign_slopes(pai(:,calculated_k+1:calculated_k+v_vec(l)),slopes(l,:));
            coef=coef+sum(a,2);
            calculated_k=calculated_k+v_vec(l);
            l=l+1;
        end

        sorted_coef=sort(coef,'descend');
        lin_min_iter=sorted_coef'*sorted_p+v_vec*consts;

        %extact est vals
        
        %rank p according to the order of the coef vector 
        vector=coef';
        [sorted_values order]=sort(vector);
        rank=zeros(1,length(vector));
        rank(order) = 1:length(vector);


        ind=max(rank)-rank+1;        
        
        opt_perm=sorted_p(ind);
        
        est_vals=zeros(P,n);
        %est_vals_in_ranges=zeros(P,n);
        for m1=1:P
            for m2=1:n
                est_vals(m1,m2)=sum((pai(:,m2)==(m1-1)).*opt_perm);
                %est_vals_in_ranges(m1,m2)=size(range,1)-sum(est_vals(m1,m2)<range);
            end
        end
        
        %We want to find for each set of parameters at which cell it lies.
        %We use the concavity of the linear approximation to find the minmal objective which determines the minimizing cell 
        sorted_est_vals=sort(est_vals);
        est_vals_tag=est_vals(1:end-1,:);
        [minimum occupied_cells]=min(slopes*est_vals_tag+repmat(consts,1,size(est_vals_tag,2)));
        
        %extract est_v_vec
        %sorted_est_vals_in_ranges=sort(est_vals_in_ranges);
        %sorted_est_vals_in_ranges=sorted_est_vals_in_ranges(1:end-1,:);
        %weights=repmat(P.^(P-2:-1:0)',1,n);
        %occupied_cells=sum((sorted_est_vals_in_ranges-1).*weights,1)+1;
        est_v_vec=zeros(1,size(v_vec,2));
        for m=1:n
              est_v_vec(occupied_cells(1,m))=est_v_vec(occupied_cells(1,m))+1;
        end
        
        if var(est_v_vec-v_vec)>0
            %test if we are not in decent direction
            if lin_min_iter-10^-3>lin_min
                prob=1;
            end
            v_vec=est_v_vec;
            lin_min=lin_min_iter;
            ent_with_appox_vals=sum(sum(-est_vals.*log2(est_vals)));
            opt_v_vec=v_vec;
            iter=iter+1;
        else
            flag=0;
            %test if we are not in decent direction
            if lin_min_iter-10^-3>lin_min
                prob=1;
            end
            lin_min=lin_min_iter;
            ent_with_appox_vals=sum(sum(-est_vals.*log2(est_vals)));
            opt_v_vec=v_vec;
        end
        if iter>10000
            flag=0;
            terminated=1;
            ent_with_appox_vals=inf;
            lin_min=inf;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [S] = mapeiainteiro_to_tuple(column_wise_matrix,P,K)
	% Maps from range 1 to P^K it's respective integer to respective tuple
	% in component form to effective BSS or ICA task
	% The -1 is necessary, since 0 is the first number and full element
	% is the last: P^K - 1

	% 48 is 0 in ascii form...
	S = dec2base(column_wise_matrix-1,P,K) - 48;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = is_linear_comb(most_recent_W,row_candidate,P)

	% This function is not explicitly being used, but a 
	% construction is used in GLICA_function (10/05/2019)
	% It was a test setup 
    k = size(most_recent_W,1);
    value = any(sum(abs(mod(generate_pai_P(k,P)(1:end-1,:) * most_recent_W,P) - row_candidate),2)
    ==0);
    
	% mod(mod(generate_pai_P(k,P)(1:end-1,:) * most_recent_W,P) - row_candidate,P)
	% mod(generate_pai_P(k,P)(1:end-1,:) * most_recent_W,P)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = GLICA_function(X,param)
	% this fuction seems really closer to america
	% and was created by Painsky et.al
	% We'll do an implementation accordingly to
	% the original paper: 
	% "Linear Independent Component Analysis over
	% Finite Fields: Algorithms and Bounds"

	% X is the observations matrix. which is K x Nobs
	% output: the separation matrix W

	P = param.P;
	K = param.K;
	PK = param.PK;
	
	% Nothing is said about Lex in Painsky's work
	% Lex = param.Lex;

    
	eqepes = 1e-9;


	vT_matrix = generate_pai_P(K,P);
    
    %america does in this order and uses Lex to reference this elements
    % flipud(vT_matrix)
    
    
	U = produtomatrizGF(vT_matrix , X,P,1,[]);
    
    marg_probs = estimate_marg_probs(U,P)';
    U_entropies = entropy_from_frequencies(marg_probs)';
    U_entropies(end) = NaN;
    
    
    [U_entropies_sorted U_entropies_sorted_index] = sort(U_entropies);
    
    
    index_entropy = 1;
    W = [];
    k=1;
    while k<=K

        if isempty(W)
                W = [ vT_matrix(U_entropies_sorted_index,:)(index_entropy,:)];
                index_entropy += 1;
                k+=1;
                continue;
        end
        
        row_candidate = vT_matrix(U_entropies_sorted_index,:)(index_entropy,:);
        
        % checks if the row candidate is a combination of pre-existing lines
        existing_lines = size(W,1);
        
        
        vT_matrix_view = vT_matrix(1:P^existing_lines,1+K-existing_lines:end);
        
        %excluding zero
        vT_matrix_view = vT_matrix_view(1:end-1,:);

        %if the row_candidate is NOT a linear combination...
        % of the existing rows in W... add it to the W matrix
        if(~ any(sum(abs(mod(vT_matrix_view * W,P) - row_candidate),2) ==0) )
            W = [W ; row_candidate];
            k+=1;
        end
        index_entropy+=1;

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [probs,S] = geravetorsinais(q,m,K,thre,Nobs)
	%auxiliary routine to generate K non-uniform sources over GF(q^m)
	%INPUT
	%K: number of sources
	%thre: non-uniformity threshold
	%Nobs: number of observations

	%OUTPUT
	%probs: PxK probability matrix (pmf) for each source
	%S: KxNobs generated sources matrix

	%Daniel Guerreiro e Silva - 12/01/2015
	P = q^m;

	probs = zeros(P,K);
	for k=1:K
	    prob=rand(P,1);
	    sprob=sum(prob);
	    prob=prob./sprob(ones(P,1),:);
	    KLD = 1 + prob'*(log(prob)./log(P)); %Kullback-Leibler Divergence
	    while(KLD<thre || max(prob)>.98 || min(prob)==0) %non-uniformity and non-degenerate requirement
	        prob=rand(P,1);
	        sprob=sum(prob);
	        prob=prob./sprob(ones(P,1),:);
	        KLD = 1 + prob'*(log(prob)./log(P));
	    end
	    probs(:,k)=prob;
	end        
	cprobs=cumsum(probs);
	RND=rand(K,Nobs);
	S=zeros(K,Nobs);
	for k=1:K
	    for cp=1:P-1
	        S(k,:)=S(k,:)+(RND(k,:)>cprobs(cp,k));
	    end
	end

	if(m>1)
	    S = S - 1;
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function A = geramatrizmistura(q,m,field,K)
	%auxiliary routine to randomly generate a KxK mixing matrix over GF(q^m)
	%INPUT
	%field: non-uniformity threshold
	%Nobs: number of observations

	%OUTPUT
	%probs: PxK probability matrix (pmf) for each source
	%S: KxNobs generated sources matrix

	%Daniel Guerreiro e Silva - 12/01/2015
	P = q^m;
	A = randi(P,K) - 1;
	A = A + diag(diag(A)==0);
	AL = tril(A);
	AU = eye(K) + triu(A,1);

	if(m>1)%adjustment if non-prime field
	    AL = AL - 1;
	    AU = AU - 1;
	end

	A = produtomatrizGF(AU,AL,q,m,field);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [probs] = generate_random_pmf(Q)
	%This routine generates a random pmf with cardinality Q
	% This routine is used to give a joint_pmf
	% to ICA experiments

	% Mateus Marcuzzo da Rosa 15/07/2019

    probs =rand(Q,1);
    sprob=sum(probs);
    probs=probs./sprob(ones(Q,1),:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [counter]=add_to_counter(counter,P)
	% It's used in the following function.
    n=size(counter,2);
    flag=1;
    while flag
       t=counter(1,n)+1;
       if t<=P
           counter(1,n)=t;
           flag=0;
       else
           counter(1,n)=0;
           if n-1>0
                n=n-1;
           else
               flag=0;
           end
       end
    end
               
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pai]=generate_pai_P(n,P)
    counter=zeros(1,n);
    %pai=counter;
    number_of_words=P^n;
    pai=zeros(number_of_words,n);
    for m=1:number_of_words
        pai(m,:)=counter;
        counter=add_to_counter(counter,P-1);
    end
    %pai=pai(1:end-1,:);      
    pai=flipud(pai);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pai]=generate_pai(n)
    pai=zeros(2^n,n);
    for m=0:n-1
        len=2^(n-m-1);
        vec_instance=[ones(len,1); zeros(len,1)];
        vec=repmat(vec_instance,2^m,1);
        pai(:,m+1)=vec;
    end
                       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Not being used, but here just for register

% function [probs] = generate_binomial_pmf(N,P_succ)
% 	% This routine generates a binomial pmf with cardinality N
% 	% Given the probability of success P
% 	% This routine is used to give a joint_pmf
% 	% to ICA experiments

% 	% Mateus Marcuzzo da Rosa 05/08/2019
% 	% Which is actually not needed, since Octave gives support for
% 	% a function does the same, and faster.

%     probs = zeros(N,1);

%     N = N-1;
%     % Note that we must use (N-1 k), not the N itself, to
%     % do a binomial with support size N

%     for i=0:(N)
%     	j = i+1;
%     	probs(j) = nchoosek(N,i)*(P_succ^i)*(1-P_succ)^(N-i);
% 	end

       

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function marg_probs = estimate_marg_probs(Y,P)
	% This function estimates the marginal probabilities given the
	% Y matrix, of the observations. It must known which prime P its
	% working

	marg_probs= [];
	for symbol=0:(P-1)
		marg_probs = [marg_probs mean(Y==symbol,2)];
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = entropy_from_frequencies(ffq, log_base=2)
	% FFQ must be already in frequency form, integer values are not count


	% IMPORTANT: THE FREQUENCIES MUST BE COLUMN-WISE
	% I.E assert(all(sum(ffq,1) == 1))


	% since it's double values and we have float problems, 
	% we will assert the way below.

	assert(all(sum(ffq,1) >= 0.95));
	assert(all(sum(ffq,1) <= 1.1));


	h=-sum(ffq.*log2(ffq+eps),1)./log2(log_base);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [slopes, consts,range]=calc_k_params(k)
    range_length=0.5/k;
    range=[0:range_length:0.5]';
    mid_range_points=[range_length/2:range_length:0.5-range_length/2]';
    slopes=log2(1-mid_range_points)-log2(mid_range_points);
    ent_at_mid_range_points=-mid_range_points.*log2(mid_range_points)-(1-mid_range_points).*log2(1-mid_range_points);
    consts=ent_at_mid_range_points-slopes.*mid_range_points;
    
    %slopes=[0.8; 0.7; 0.6;];
    %consts=[0.2; 0.4; 0.5;];



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Originally, add_to_counter, but we already have a function
% here with the same name, which gives support to another function
% so, this is the version _2.
function [counter]=add_to_counter_v2(counter,k,P)
    n=size(counter,2);
    flag=1;
    while flag
       t=counter(1,n)+1;
       if t<=k
           counter(1,n)=t;
           flag=0;
       else
           counter(1,n)=1;
           if n-1>0
                n=n-1;
           else
               flag=0;
           end
       end
    end
               
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [slopes, consts]=calc_k_params_P(k,P)

    range_length=0.5/k;
    range=[0:range_length:0.5]';
    mid_range_points=[range_length/2:range_length:0.5-range_length/2]';
    mid_range_points_in_all_dims= repmat(mid_range_points,1,P-1);
    %number_of_points=size(mid_range_points,1);
    slopes=inf*ones(k^(P-1),P-1);   %in total we have k^(P-1) cells, each one contains P-1 slope parameters (gradients)
    consts=zeros(k^(P-1),1);   %in total we have k^(P-1) cells, each contains a single value.
    
    counter=ones(1,P-1);  %this counter indicates the desired point at each axis 
    for m1=1:size(slopes,1)
       p= mid_range_points_in_all_dims(counter); %p holds the mid range points at the current counter. These are actually the probabilities according to which we calculate the entropy. 
       if sum(p)<1
            all_p=[p 1-sum(p)];
            ent_at_p=sum(-all_p.*log2(all_p));
            for m2=1:P-1
                mid_range_point=p(m2);
                slopes(m1,m2)=log2(all_p(end))-log2(mid_range_point);
                
       
            end
            consts(m1)=ent_at_p-p*slopes(m1,:)';
       end
       counter=add_to_counter_v2(counter,k,P);
    end
    %calc correct range (where the gradient intersect)
%     slopes_of_single_cell=slopes(1:size(range,1)-1,end);
%     consts_of_single_cell=consts(1:size(range,1)-1,end);
%     slopes_difference=slopes_of_single_cell(1:end-1)-slopes_of_single_cell(2:end);
%     consts_difference=consts_of_single_cell(2:end)-consts_of_single_cell(1:end-1);
%     correct_range=consts_difference./slopes_difference;
%     range=[0; correct_range; 0.5;];
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = allVL1nonrecurs(n, L1)
	% function v=allVL1eq(n, L1);
	% INPUT
	%    n: length of the vector
	%    L1: desired L1 norm
	% OUTPUT:
	%    if head is not defined
	%      v: (m x n) array such as sum(v,2)==L1
	%         all elements of v is naturel numbers {0,1,...}
	%         v contains all (=m) possible combinations
	%         v is (dictionnary) sorted
	% Algorithm:
	%    NonRecursive

	% Chose (n-1) the splitting points of the array [0:(n+L1)]
	s = nchoosek(1:n+L1-1,n-1);
	m = size(s,1);

	s1 = zeros(m,1,class(L1));
	s2 = (n+L1)+s1;

	v = diff([s1 s s2],1,2); % m x n
	v = v-1;

end % allVL1nonrecurs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [permuted_p,opt_perm,est_vals,min_ent]=BICA_function(p,min_k,max_k)

	min_ent=inf;
	permuted_p=inf;
	opt_perm=inf;
	est_vals=inf;
	n=log2(length(p));


	[sorted_p sorted_p_ind]=sort(p);
	%sorted_p=p(sorted_p_ind);
	[pai]=generate_pai(n); 


	ent_with_lin_params=zeros(max_k-min_k+1,1);
	lin_min_vec=zeros(max_k-min_k+1,1);


	for k=min_k:max_k
	    k;
	    [slopes, consts,range]=calc_k_params(k);

	    lin_min=inf;
	    opt_v_vec=NaN;

	    %v = allVL1(k, n);
	    v = allVL1nonrecurs(k,n); %returns all combinations of n balls in k boxes
	    v=fliplr(v);

	    for iter=1:size(v,1)
	        coef=zeros(2^n,1);
	        calculated_k=0;
	        v_vec=v(iter,:);
	        for l=1:size(v,2)
	            a=slopes(l)*pai(:,calculated_k+1:calculated_k+v_vec(l));
	            coef=coef+sum(a,2);
	            calculated_k=calculated_k+v_vec(l);
	        end
	        sorted_coef=flipud(sort(coef));
	        lin_min_iter=sorted_coef'*sorted_p+v_vec*consts;
	        if lin_min_iter<lin_min
	            lin_min=lin_min_iter;
	            opt_v_vec=v_vec;
	        end
	    end


	    %calc est_vals
	    coef=zeros(2^n,1);
	    calculated_k=0;
	    for l=1:size(v,2)
	        a=slopes(l)*pai(:,calculated_k+1:calculated_k+opt_v_vec(l));
	        coef=coef+sum(a,2);
	        calculated_k=calculated_k+opt_v_vec(l);
	    end

	    
	    %rank p according to the order of the coef vector 
	    vector=coef';
	    [sorted_values order]=sort(vector);
	    rank=zeros(1,length(vector));
	    rank(order) = 1:length(vector);
	    
	    
	    ind=max(rank)-rank+1;
	    opt_p=sorted_p(ind);
	    est_vals=sum(pai.*repmat(opt_p,1,n));
	    lin_min_vec(k-min_k+1)=lin_min;
	    ent_with_lin_params(k-min_k+1)=sum(-est_vals.*log2(est_vals)-(1-est_vals).*log2(1-est_vals));
	    if k==min_k
	        opt_perm= sorted_p_ind(ind');
	        %sorted_p=p(sorted_p_ind);
	        permuted_p=opt_p;
	        min_ent=ent_with_lin_params(k-min_k+1);
	        
	    else
	        if  ent_with_lin_params(k-min_k+1)<min_ent
	            opt_perm= sorted_p_ind(ind');
	            permuted_p=opt_p;
	            min_ent=ent_with_lin_params(k-min_k+1);  
	        end
	    end
	end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M]=assign_slopes(A,probs)
    % tic
    % M=zeros(size(A,1),size(A,2)); 
    % for m1=1:size(A,1)
    %     for m2=1:size(A,2)
    %         if A(m1,m2)+1<=size(probs,2)
    %             M(m1,m2)=probs(1,A(m1,m2)+1);
    %         end
    %     end
    % end
    % toc
    % tic
    M=zeros(size(A,1),size(A,2));
    for m1=1:size(probs,2)
        M=M+probs(1,m1)*(A==(m1-1));
    end
    % toc
    
                 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There's a cleaner function reading of AMERICA, but we'll
% put it here later (31/07/2020)

function B=america(PPx,param)
	%The Ascending Minimization of EntRopies for ICA
	%(AMERICA) algorithm
	%input: PPx - the (estimated) probabilities tensor
	%output: B - the estimated separating matrix

	% global P K PK r Lex

	P = param.P;
	K = param.K;
	PK = param.PK;
	Lex = param.Lex;
	r = param.r;

	eqeps=1e-9; %a threshold for deciding
	            %equal entropies


	%K-D fft for obtaining the characteristic tensor
	fPPyn=fftn(reshape(PPx,P*ones(1,K)));
	fPPy=fPPyn(:);
	        
	%obtain the characteristic vectors of 
	%the linear combinations
	qf=ones(P,PK);
	qf(2,:)=fPPy;
	if P>2
	    qf(P,:)=conj(fPPy);
	    for m=2:P/2
	        mLex=mod(m*Lex,P);
	        qf(m+1,:)=fPPy(r*mLex+1);
	        qf(P+1-m,:)=conj(qf(m+1,:));
	    end
	end

	%translate characteristic vectors into probabilities
	%vectors and then into entropies
	ffq=ifft(qf);
	ffq=max(ffq,eps);
	h=-sum(ffq.*log2(ffq+eps),1);
	%mark irrelevant entropies (such as the one related
	%to the all-zeros (trivial) combination, and subsequent
	%"used" entropies - with a NaN
	h(1)=NaN;

	B=[];
	k=1;
	%sorted entropies (ascending order)
	[sh shix]=sort(h);
	inh=1;
	while k<=K
	    vh=sh(inh);
	    mix=shix(inh);
	    for itry=inh+1:PK
	        if abs(sh(itry)-vh)>eqeps, break; end
	    end
	    %randomized selection in case of a tie
	    neq=itry-inh;
	    if neq>1
	        ipick=floor(rand*neq);
	        pinh=inh+ipick;
	        tmph=sh(inh);
	        tmpi=shix(inh);
	        sh(inh)=sh(pinh);
	        shix(inh)=shix(pinh);
	        sh(pinh)=tmph;
	        shix(pinh)=tmpi;
	    end
	    %test if the selected is not a linear combination
	    %of the previous ones
	    mix=shix(inh);
	    b=Lex(:,mix);
	    Bb=[B b];
	    TLex=Lex(1:k,2:P^k);
	    test0=mod(Bb*TLex,P);
	    if ~any(sum(test0,1)==0)    %not a linear combination
	        B=Bb;   
	        k=k+1;
	    end
	    inh=inh+1;
	end

	B=B';
end;

pkg load linear-algebra
function new_pmf_conj = order_permutation_algorithm(pmf_conj)
    line_pmf_conj = pmf_conj(:)
    [line_pmf_conj, indices] = sort(line_pmf_conj)

    ndimensions = ndims(pmf_conj)
    new_pmf_conj = zeros(size(pmf_conj))

    % NO PYTHON: Checar o BSS_functions.py!
    % for index,p in enumerate(product(range(pmf_conj.shape[0]),repeat=ndim)):
    %   new_pmf_conj[p] = line_pmf_conj[index]
    % iterar sob o produto cartesiano.
    % pmf_conj.shape[0] é size(pmf_conj)(0), aqui.
    % para um pmf_conj de shape/size 2 2 2

    % queríamos iterar sobre o produto cartesiano de [1 2] x [1 2] x [1 2]
    % que nos daria, primeiro:
    %  1 1
    %  1 2
    %  2 1
    %  2 2

    % depois 
    % 1 1 1
    % 1 1 2
    % etc

    % o problema aqui, é que o octave só aceita como cartprod(arg1,arg2...)
    % e não como um argumento só, que é uma lista.

    % então devemos rodar alguma espécie de loop. no momento não consigo fazer.
    % TODO!!!

end