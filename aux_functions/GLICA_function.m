function W = GLICA_function(X,P,K)	
	% An implementation accordingly to
	% the original paper:
	% "Linear Independent Component Analysis over
	% Finite Fields: Algorithms and Bounds"

	% X is the observations matrix. which is K x Nobs
	% output: the separation matrix W

    PK = P^K;
    vT_matrix = single(dec2base(PK-1:-1:0,P,K) - 48);
    U_entropies = zeros(1,PK);
    
    for it=1:PK-1
        u = product_GFmatrix(vT_matrix(it,:),X,P,1,[]);        
        marg_prob = estimate_marg_probs(u,P)';
        U_entropies(it) = entropy_from_frequencies(marg_prob)';        
    end
    U_entropies(end) = NaN;
    [~, U_entropies_sorted_index] = sort(U_entropies);

    
    index_entropy = 1;
    W = [];
    k=1;
    while k<=K

        if isempty(W)
                tmp = vT_matrix(U_entropies_sorted_index,:);
                W =  tmp(index_entropy,:);
                index_entropy = index_entropy + 1;
                k= k + 1;
                continue;
        end

        tmp = vT_matrix(U_entropies_sorted_index,:);
        row_candidate = tmp(index_entropy,:);

        % checks if the row candidate is a combination of pre-existing lines
        existing_lines = size(W,1);


        vT_matrix_view = vT_matrix(1:P^existing_lines,1+K-existing_lines:end);

        %excluding zero
        vT_matrix_view = vT_matrix_view(1:end-1,:);

        %if the row_candidate is NOT a linear combination...
        % of the existing rows in W... add it to the W matrix
        if(~ any(sum(abs(mod(vT_matrix_view * W,P) - row_candidate),2) ==0) )
            W = [W ; row_candidate];
            k = k + 1;
        end
        index_entropy= index_entropy + 1;

    end

end
