clear;
close all;

some_primes = [5];
n_sources = [2 3 4];
n_samples = 10.^(2:6);
n_trials = 50;
qica_min_k = 2;
qica_max_k = 10;
sig_digits = 4;

addpath('aux_functions/');
if isempty(gcp('nocreate')) %MATLAB
    parpool(10);
end
sim_start_time = datetime(); %MATLAB
%sim_start_time = localtime(time()); %OCTAVE

% Must be a cel array, so we can do a strfind...
algorithms_names = {'america';'sa4ica';'QICA';'GLICA';'order'};
n_algorithms = length(algorithms_names);
the_algorithms = 1:n_algorithms;

% distributions_names = {'Zipf';'Binomial p=0.2';'Binomial p=0.5';'Random'};
distributions_names = {'Zipf', 'Binomial p=0.2', 'Random'};
the_distributions = 1:length(distributions_names);


space = [length(distributions_names), length(some_primes), length(n_sources), length(n_samples), n_trials];

n_cases = prod(space); %number of tested scenarios per algorithm

america_trial_time = zeros(n_cases,1);
america_total_corr_results = zeros(n_cases,1);

sa4ica_trial_time = zeros(n_cases,1);
sa4ica_total_corr_results = zeros(n_cases,1);

QICA_trial_time = zeros(n_cases,1);
QICA_total_corr_results = zeros(n_cases,1);

GLICA_trial_time = zeros(n_cases,1);
GLICA_total_corr_results = zeros(n_cases,1);

order_trial_time = zeros(n_cases,1);
order_total_corr_results = zeros(n_cases,1);

total_time = tic;
% s = RandStream('mt19937ar','Seed',trial_i);

parfor idcase = 1:n_cases    
    [dist_i, p_i, k_i, t_i, trial_i] = ind2sub(space,idcase);
    P = some_primes(p_i);
    K = n_sources(k_i);
    Nobs = n_samples(t_i);
    
    PK = P^K;
    %% Directly extracted from america code

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

    if(dist_i == find(strcmp(distributions_names,'Zipf')))
        %First we are going Zipf like Painsky's code
        %Zipf distribution
        s=1.05; %zipf distribution parameter
        zipf_p=[1:1:PK];
        zipf_p=1./(zipf_p.^s);
        zipf_p=zipf_p'/sum(zipf_p);
        joint_pmf  = zipf_p';
        generated_integers = gerainteiros_from_probs(joint_pmf,P,K,Nobs);
    elseif(dist_i == find(strcmp(distributions_names,'Binomial p=0.2')))
        generated_integers = binornd(PK-1,.2,Nobs,1)+1;
    elseif(dist_i == find(strcmp(distributions_names,'Binomial p=0.5')))
        generated_integers = binornd(PK-1,.5,Nobs,1)+1;
    elseif(dist_i == find(strcmp(distributions_names,'Random')))
        joint_pmf  = single(generate_random_pmf(PK)');
        generated_integers = gerainteiros_from_probs(joint_pmf,P,K,Nobs);
    end

    X = mapeiainteiro_to_tuple(generated_integers,P,K);
    X = X';
   
    % calculo tensor de probabilidades - AMERICA  
    Px = single(zeros(1,PK));
    for t=1:Nobs
        Px(generated_integers(t)) = Px(generated_integers(t)) + 1;
    end

    Px=Px/Nobs;

    h_joint=-sum(Px(Px>0).*log2(Px(Px>0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% AMERICA ALGORITHM EVALUATION %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    start_time = tic;
    [Wm] = america(Px,P,K,PK,Lex,r);
    america_trial_time(idcase) = toc(start_time);

    Y = produtomatrizGF(Wm,X,P,1,[]);

    h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
    h_marg = sum(h_marg(:));
    america_total_corr_results(idcase) = h_marg-h_joint;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% SA4ICA ALGORITHM EVALUATION %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = tic;

    % the decode_ function, makes a confusion with a pre-existing
    % decode function. so I renamed it to decode_. It is part of sa4ica
    [Wsa] = sa4ica_decode(Px,P,K,Lex,r,0.995,5);
    sa4ica_trial_time(idcase) = toc(start_time);

    Y = produtomatrizGF(Wsa,X,P,1,[]);
    h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
    h_marg = sum(h_marg(:));
    sa4ica_total_corr_results(idcase) = h_marg-h_joint;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  QICA ALGORITHM EVALUATION  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = tic;

    [opt_est_vals,opt_appox_lin_min2,opt_appox_ent_with_appox_vals2,opt_perm,opt_v_vec]=QICA_function(K,P,Px',0,qica_min_k,qica_max_k,1000);
    QICA_trial_time(idcase) = toc(start_time);

    [Yqica,] = mapeiapermutacao(X, opt_perm,P,K);

    h_marg=entropy_from_frequencies(estimate_marg_probs(Yqica,P)');
    h_marg = sum(h_marg(:));

    QICA_total_corr_results(idcase) = h_marg-h_joint;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Order Permutation ALGORITHM EVALUATION  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(P==2)
        start_time = tic;

        [Yorder, opt_perm] = order_perm_function(Px,X,P,K);

        order_trial_time(idcase) = toc(start_time);

        h_marg=entropy_from_frequencies(estimate_marg_probs(Yorder,P)');
        h_marg = sum(h_marg(:));

        order_total_corr_results(idcase) = h_marg-h_joint;
    else
        order_trial_time(idcase) = NaN;
        order_total_corr_results(idcase) = NaN;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  GLICA ALGORITHM EVALUATION  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = tic;
    [Wglica] = GLICA_function(X,P,K);
    GLICA_trial_time(idcase) = toc(start_time);

    Y = produtomatrizGF(Wglica,X,P,1,[]);
    h_marg=entropy_from_frequencies(estimate_marg_probs(Y,P)');
    h_marg = sum(h_marg(:));
    GLICA_total_corr_results(idcase) = h_marg-h_joint;


end
toc(total_time)



[america_mean_trial_time, america_trial_time] = summary_statistics(america_trial_time,space, sig_digits);
[america_mean_total_corr_results, america_total_corr_results] = summary_statistics(america_total_corr_results,space, sig_digits);

[sa4ica_mean_trial_time, sa4ica_trial_time] = summary_statistics(sa4ica_trial_time,space, sig_digits);
[sa4ica_mean_total_corr_results, sa4ica_total_corr_results] = summary_statistics(sa4ica_total_corr_results,space, sig_digits);

[QICA_mean_trial_time, QICA_trial_time] = summary_statistics(QICA_trial_time,space, sig_digits);
[QICA_mean_total_corr_results, QICA_total_corr_results] = summary_statistics(QICA_total_corr_results,space, sig_digits);

[GLICA_mean_trial_time, GLICA_trial_time] = summary_statistics(GLICA_trial_time,space, sig_digits);
[GLICA_mean_total_corr_results, GLICA_total_corr_results] = summary_statistics(GLICA_total_corr_results,space, sig_digits);

[order_mean_trial_time, order_trial_time] = summary_statistics(order_trial_time,space, sig_digits);
[order_mean_total_corr_results, order_total_corr_results] = summary_statistics(order_total_corr_results,space, sig_digits);


% saves with the date (day/month/year) and the hour: hh:mm
% start and ending times
% start_time_str = strftime('sim_data_start_%d_%m_%Y_%H_%M',sim_start_time);%OCTAVE
% saved_sim_str = strftime('_end_%d_%m_%Y_%H_%M_sim_total_corr_pure_ICA',localtime(time()));%OCTAVE
start_time_str = 'sim_data_start_' + string(sim_start_time,'yyyy_MM_dd_HH:mm'); %MATLAB
saved_sim_str = '_end_' + string(datetime(),'yyyy_MM_dd_HH:mm') + '_sim_total_corr_pure_ICA';%MATLAB

saved_sim = sprintf('sim_data/%s%s',start_time_str,saved_sim_str);
save(saved_sim)
