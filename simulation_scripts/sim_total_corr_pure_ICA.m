% NOT READY YET!!!! (30/08/2020)
% TO-DO: refactor this
% transform it to a script file
% This code is greatly inspired in Painsky's code of QICA.m
% not QICA_function.m !
% 16/09/2019
% adding these comments to prevent future problems
% This is the pure Zipf pmf simulation.

%% By now it englobes just the pure_zipf experiment
some_primes = [2];
n_sources = 3;
n_samples = [128, 256, 512, 1024];
n_trials = 10;
qica_min_k = 4;
qica_max_k = 8;
%the_distribution =

addpath('aux_functions/');
% sim_start_time = datetime(); %MATLAB
sim_start_time = localtime(time()); %OCTAVE

% Must be a cel array, so we can do a strfind...
algorithms_names = {'america';'sa4ica';'QICA';'GLICA';'order'};

the_algorithms = 1:length(algorithms_names);


space = [ length(some_primes), length(n_sources), length(n_samples),...
 length(the_algorithms) ];

n_cases = prod(space);

trial_time = zeros( [ size(zeros( [ space ] )), n_trials] );

total_corr_results = zeros( [size(trial_time) ] );

total_time = tic;

for p_i = 1:length(some_primes)
    for k_i = 1:length(n_sources)
%        The code below was for debugging
%        all_possible_tuples = generate_pai_P(n_sources(k_i),some_primes(p_i))'
%        apt = all_possible_tuples;
        for t_i = 1:length(n_samples)

            P = some_primes(p_i);
            K = n_sources(k_i);
            Nobs = n_samples(t_i);
            fprintf('\n');


            for trial=1:n_trials

                PK = P^K;
                %%%%%%% from super_comparativo.m %%%%%%
                %% Which was directly extracted from america code
                %% itself

                % construct a "Lexicon": Lex(:,n) is the (n-1)-th
                Lex=single(zeros(K,PK));
                nvec=single(0:PK-1);
                for k=1:K-1
                    Lex(k,:)=mod(nvec,P);
                    nvec=floor(nvec/P);
                end
                Lex(K,:)=nvec;
                %construct an "index-vector translation" vector:
                r=P.^(0:K-1); %AMERICA parameter



                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %First we are going Zipf like Painsky's code
                q=P;
                n=K; %dimension of the random vector
                %Zipf distribution
                s=1.05; %zipf distribution parameter
                zipf_p=[1:1:q^n];
                zipf_p=1./(zipf_p.^s);
                zipf_p=zipf_p'/sum(zipf_p);

                joint_pmf  = zipf_p';

                generated_integers = gerainteiros_from_probs(joint_pmf,P,K,Nobs);
                X = mapeiainteiro_to_tuple(generated_integers,P,K);
                X = X';


                Y = mapeiapermutacao([],X,1:PK,P,K);
                % A permutacao trivial em X deve dar o pr�prio X
                % Comentar este c�digo quando seguro de consist�ncia
                assert( all( all( X==Y ) ) );

                % THe -1 is necessary,since the first element is always zero.
                estimated_joint_from_samples = estimate_marg_probs(generated_integers'-1,PK);

                h_joint_from_samples = entropy_from_frequencies(estimated_joint_from_samples');

%                disp('gerou joint pmf');

                parameters.P = P;
                parameters.K = K;
                parameters.PK = PK;
                parameters.Lex = Lex;
                parameters.r = r;



                % calculo tensor de probabilidades - AMERICA
                idx=r*X;
                Px = single(zeros(1,PK));
                for t=1:Nobs
                    Px(idx(t)+1) = Px(idx(t)+1) + 1;
                end

                Px=Px/Nobs;

                h_joint=-sum(Px(Px>0).*log2(Px(Px>0)));

                % a toler�ncia aqui � de 0.2
                assert( h_joint <= h_joint_from_samples + 0.2);
                assert( h_joint >= h_joint_from_samples - 0.2);
%                disp('h_joint - h_joint_from_samples');
%                disp(h_joint - h_joint_from_samples);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% AMERICA ALGORITHM EVALUATION %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                algo_i = find(strcmp(algorithms_names,'america'));

                start_time = tic;
                [Wm] = america(Px,parameters);
                trial_time(p_i,k_i,t_i,algo_i,trial) = toc(start_time);

                % This is the demixing matrix, we are going to do U*A
                % and check the number of hits. If everyone hits, we count it
%		        U = produtomatrizGF(Wm,A,P,1,[]);
%                U_america = U

                % there's some legacy code happening here, we'll leave it here.
                % it resembles supercomparativo.m code
%		        if(1>1)
%		            Z = (U>-1); %null element in GF(P^1)
%		        else
%		            Z = (U>0);
%		        end

%		        hits = sum(sum(Z,2)==1);
%		        if(hits == K)
%		        	bss_succ_rate(p_i,k_i,t_i,algo_i) += 1;
%		        end

                Y = produtomatrizGF(Wm,X,P,1,[]);
%                Y_possible_tuples_america = produtomatrizGF(U_america,all_possible_tuples,P,1,[])

                % estimate_marg_probs(Y,P)' must be alike with 'the_pmfs'

                h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
                h_marg = sum(h_marg(:));
                total_corr_results(p_i,k_i,t_i,algo_i,trial) = h_marg-h_joint;


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% SA4ICA ALGORITHM EVALUATION %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                algo_i = find(strcmp(algorithms_names,'sa4ica'));

                start_time = tic;

                % the decode_ function, makes a confusion with a pre-existing
                % decode function. so I renamed it to decode_. It is part of sa4ica
                [Wsa] = sa4ica_decode(Px,parameters,0.995,5);
                trial_time(p_i,k_i,t_i,algo_i,trial) = toc(start_time);

                %once again, we are going to test the number of hits
%		        U = produtomatrizGF(Wsa,A,P,1,[]);
%                Usa4ica = U

%		        if(1>1)
%		            Z = (U>-1); %null element in GF(q^m)
%		        else
%		            Z = (U>0);
%		        end

%		        hits = sum(sum(Z,2)==1);
%		        if(hits == K)
%		        	bss_succ_rate(p_i,k_i,t_i,algo_i) += 1;
%		        end

                Y = produtomatrizGF(Wsa,X,P,1,[]);
%                Y_possible_tuples_sa4ica = produtomatrizGF(Usa4ica,all_possible_tuples,P,1,[])
                h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
                h_marg = sum(h_marg(:));
                total_corr_results(p_i,k_i,t_i,algo_i,trial) = h_marg-h_joint;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  QICA ALGORITHM EVALUATION  %%%%%
                %%%%% WE'VE SENT THIS TO EXPERIMENT 2 %
                %%%%%        WHICH IS PURE ICA     %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                algo_i = find(strcmp(algorithms_names,'QICA'));

                 start_time = tic;


                [opt_est_vals,opt_appox_lin_min2,opt_appox_ent_with_appox_vals2,opt_perm,opt_v_vec]=QICA_function(K,P,Px',0,qica_min_k,qica_max_k,1000);


                trial_time(p_i,k_i,t_i,algo_i,trial) = toc(start_time);

                [Yqica,] = mapeiapermutacao([], X, opt_perm,P,K);

                 h_marg=entropy_from_frequencies(estimate_marg_probs(Yqica,P)');
                 h_marg = sum(h_marg(:));

                 total_corr_results(p_i,k_i,t_i,algo_i,trial) = h_marg-h_joint;

                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %%%%%  Order Permutation ALGORITHM EVALUATION  %%%%%
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 algo_i = find(strcmp(algorithms_names,'order'));

                 start_time = tic;

                 [Yorder, opt_perm] = order_perm_function(Px,X,P,K);

                 trial_time(p_i,k_i,t_i,algo_i,trial) = toc(start_time);

                 h_marg=entropy_from_frequencies(estimate_marg_probs(Yorder,P)');
                 h_marg = sum(h_marg(:));

                 total_corr_results(p_i,k_i,t_i,algo_i,trial) = h_marg-h_joint;


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  GLICA ALGORITHM EVALUATION  %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                algo_i = find(strcmp(algorithms_names,'GLICA'));

                start_time = tic;
                [Wglica] = GLICA_function(X,parameters);
                trial_time(p_i,k_i,t_i,algo_i,trial) = toc(start_time);

                Y = produtomatrizGF(Wglica,X,P,1,[]);
                h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
                h_marg = sum(h_marg(:));
                total_corr_results(p_i,k_i,t_i,algo_i,trial) = h_marg-h_joint;


            end
        end
    end
end
toc(total_time)


mean_trial_time = mean(trial_time, 5);
mean_total_corr_results = mean(total_corr_results, 5);


% saves with the date (day/month/year) and the hour: hh:mm
% start and ending times
start_time_str = strftime('sim_data_start_%d_%m_%Y_%H_%M',sim_start_time);%OCTAVE
saved_sim_str = strftime('_end_%d_%m_%Y_%H_%M_sim_total_corr_pure_ICA',localtime(time()));%OCTAVE
% start_time_str = 'sim_data_start_' + string(sim_start_time,'yyyy_MM_dd_HH:mm'); %MATLAB
% saved_sim_str = '_end_' + string(datetime(),'yyyy_MM_dd_HH:mm') + '_sim_total_corr_pure_ICA';%MATLAB

%plot(1:4,ans(1,1,:,1),1:4,ans(1,1,:,2),1:4,ans(1,1,:,3),1:4,ans(1,1,:,4),1:4,ans(1,1,:,5));legend('america','sa4ica','QICA','GLICA','order');

saved_sim = sprintf('sim_data/%s%s',start_time_str,saved_sim_str);
save(saved_sim)
