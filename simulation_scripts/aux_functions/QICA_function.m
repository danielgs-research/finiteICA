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


%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%



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



%%%%%%%%%%%%%%%%


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


%%%%%%%%%%


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
