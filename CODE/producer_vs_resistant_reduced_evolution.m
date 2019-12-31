%% setup
scoring_type = 'loser_remains_winner_gets_rest'; % 'winner_gets_all', 'loser_dies_winner_gets_rest', 'loser_remains_winner_gets_rest'
type_resist = 'max'; % max , plus, resist
weight_dist = 'evenly'; %evenly, half_producer 
decay = 1; % antibiotic decay (1) quick decay; (0) no decay ;
start_option = 'zero'; % the starting phenotype zero or rand
Cost = [0.005 0.01] ; % resistance and production costs
%Cost = [0.01, 0.03]; %[0.005 0.01] ; % resistance and production costs
type_mutation = 'gamma'; %Type of distribution for the mutation. either nominal or gamma
switch type_mutation
    case 'normal'
        Mut_size = [-0.1 -0.1]; % average size of resistant and production mutations (typically should be <=0)
        Mut_size_std = [0.1 0.1]; % standard deviation of resistant and production mutations
    case 'gamma'
        a_gamma = [3 3]; % a in gamma distribution
        b_gamma = [0.2 0.2]; % b (=1/beta) in gamma distribution
        move_fig = 0;
end
Mut_0 = [0.01 0.01] ; % chance of null mutations causing complete loss of resistant(1) or production(2) 
minimal_phen = [0.01  0];
time_limit = 0.8; %0.8;
self_competition = 0;
competition_between_resistants = 0;
save_fig  = true;

N_P = 1;  % number of producers 
N_K = 1; %number of antibiotics
N_R = 1; % number of resistants   
maxit = 1e4; % max number of fixations   
num_eff_producers = 2; % the effective number of inhibiting producer
for i_num_comp = 2:20
    fig_num = 100+i_num_comp*5;

%%
if strcmp(start_option,'chosen')
    intial_production = 0.5*ones(N_P,1);
    intial_resistance = 0.4*ones(N_R,1);
    intial_resistance(1) = 0.4 ;%0.25; % 0.25 - resistant to 0.5
    intial_resistance(2) = 0.4 ;%0.25; % 0.25 - resistant to 0.5
    intial_resistance(3) = 0.4 ;%0.25; % 0.25 - resistant to 0.5
end
%%
dt_print = 10000;
dt_view = 100000;
%%
N= N_R + N_P;

switch start_option
    case 'rand'
        Phen = rand(N_K,2,N);
    case 'zero'
        Phen = zeros(N_K,2,N) ; %1:Res, 2:Production
        Phen(:,1,N_P+1:N) = minimal_phen(1); %1:Res, 2:Production
end

clear growth_curves
growth_curve_filename = ['CODE/saved_profiles' filesep 'growth_curves_log100.mat'] ; %path to database of growth curves
if exist(growth_curve_filename,'file')
    load(growth_curve_filename,'growth_curves') ;
else
    growth_rate_ratio = logspace(-4,0,100) ;
    epsilon = 1e-3 ;
    n0 = 1e-3;
    num_timepoints = 1000 ;
    growth_curves = create_growth_curve_database(growth_rate_ratio,num_timepoints,n0,epsilon) ;
    save(growth_curve_filename,'growth_curves') ;
end

t = 0 ; % number of cycles
it = 1 ; % number of fixation events
t_v = nan(maxit,1) ;
improvement = nan(maxit,1) ; %saves how beneficial was the mutation
Phen_v = nan(N_K,2,N,maxit) ; % keeps all phenotypes versus time
Phen_v(:,:,:,it) = Phen ; 
t_v(it,1) = t ;
max_rounds = 1e8;
i_round = 0; 

%% run
cost_matrix = nan(N,N);
weight_matrix = nan(1,N);

switch weight_dist
    case 'evenly'
        weight_matrix(1,1:N) = 1.0/(N);
    case 'half_producer'
        weight_matrix(1,1:N_P) = 0.5/(N_P);
        weight_matrix(1,N_P+1:N) = 0.5/(N_R);
end

mutation_dist = 'equal_to_pop_size';
switch mutation_dist
    case 'equal_to_pop_size'
        mutation_distribution= weight_matrix;
    case 'different_rates'
        mutation_distribution = nan(1,N);
        production_rate = 0.5;
        mutation_distribution(1,1:N_P) = production_rate/(N_P);
        mutation_distribution(1,N_P+1:N) = (1-production_rate)/(N_R);
end

for i = 1:N
    for j= 1:i
        g1 = growth_rate(Phen(:,:,i), Cost);
        g2 = growth_rate(Phen(:,:,j), Cost);
        [cur_growth] = interp_growth_curve(g1,g2, growth_curves);
        g = [g1 g2];
        [y_i,y_j] = single_droplet(Phen(:,:,i),Phen(:,:,j),g,type_resist, scoring_type, decay, time_limit, cur_growth);
        cost_matrix(i,j) = y_i;
        cost_matrix(j,i) = y_j;
    end
end

%%
while (it<maxit)&&(t<max_rounds)
    cur_n = randsample(N,1,1,mutation_distribution);
    %cur_n = mod(t,N)+1;
    if cur_n<=N_P
        p = 2; %production
        k = randi(N_K);
    else 
        p = 1; %resistance
        k = randi(N_K);
    end
    WT = Phen(:,:,cur_n) ;
    MT = WT ;
    % mutate
    P0 = MT(k,p) ;
    switch type_mutation
        case 'normal'
            t = t + 1;
            P0 = P0 + Mut_size(p) + Mut_size_std(p)*randn ; 
            P0 = P0 * (rand>Mut_0(p));
            P0 = max(P0,minimal_phen(p));
        case 'gamma'
            t = t + 1;
            P0 = P0 + (a_gamma(p)-1)*b_gamma(p) - gamrnd(a_gamma(p),b_gamma(p)) ; %to make sure mode is 0
            P0 = P0 * (rand>Mut_0(p));
            P0 = max(P0,minimal_phen(p));
    end
    
    MT(k,p) = P0; %max(P0,minimal_phen(p));
  
    if all(all(WT == MT))
        continue;
    end
    % calc average fitness
    fWT = 0 ;
    fMT = 0 ;
    M_cost_matrix = cost_matrix;
    for i = 1:N
        if i==cur_n && (~self_competition)
            continue
        elseif (cur_n)>N_P && i>N_P && (~competition_between_resistants)
            continue
        elseif (cur_n)<=N_P && i<=N_P && (~competition_between_resistants)
            continue
        else
            fWT = fWT + weight_matrix(i)*cost_matrix(cur_n,i);

            g1 = growth_rate(MT, Cost);
            g2 = growth_rate(Phen(:,:,i), Cost);
            g = [g1 g2];
            [cur_growth] = interp_growth_curve(g1,g2, growth_curves);
            [y_i,y_j] = single_droplet(MT,Phen(:,:,i),g,type_resist, scoring_type, decay, time_limit, cur_growth);
            %[y1,y2,~,~] = single_droplet_with_solver(MT, Phen(:,:,i), g, type_resist, scoring_type, decay, time_limit);

            M_cost_matrix(cur_n,i) = y_i;
            M_cost_matrix(i,cur_n) = y_j;

            fMT = fMT + weight_matrix(i)*y_i;
        end
    end
    % fixation
    if cur_n<=N_P
        threshold = fMT/fWT - 1;
    else
        f_MT_0 = 0.44-0.0074*MT(k,p); % MT fitness against non producing producer 
        f_WT_0 = 0.44-0.0074*WT(k,p); % WT fitness against non producing producer
        %this fourmala was taken from fitting the function fitness_r(0,r) and needs to be adjusted if competition parameter change
        threshold = (f_MT_0*(i_num_comp-num_eff_producers)+num_eff_producers*fMT)/(f_WT_0*(i_num_comp-num_eff_producers)+num_eff_producers*fWT) - 1;
    end
    if (rand < threshold)
        it = it + 1 ;
        Phen(:,:,cur_n) = MT ;
        Phen_v(:,:,:,it) = Phen ;
        t_v(it,1) = t ;
        improvement(it,1) = fMT/fWT ;
        cost_matrix = M_cost_matrix;
        
        if self_competition
            g1 = growth_rate(MT, Cost);
            g = [g1 g1];
            [cur_growth] = interp_growth_curve(g1,g1, growth_curves);
            [y_i,~] = single_droplet(MT,MT,g,type_resist, scoring_type, decay, time_limit, cur_growth);
            cost_matrix(cur_n,cur_n) = y_i;
        end
    else
        %threshold
    end
    
    if ~mod(t,dt_print), disp([t,it]); end
    if ~mod(t,dt_view)
        if N_K>1
            plotPhenMultipleAntibiotics(fig_num, N_K, N, Phen_v,t_v)
        else
            plotPhen(fig_num, N_P, N, Phen_v, t_v)
        end
    end
    
    if ~still_prducting && stop_flag_num_peak && any(Phen(:,2,1:N_P)>eps_production)
        already_produced = already_produced + 1;
        still_prducting = 1;
    end
    if still_prducting  && all(Phen(:,2,1:N_P)<eps_production)
        still_prducting = 0;
        if already_produced>=stop_flag_num_peak
            break
        end
    end
end
if t>=max_rounds
    it = it + 1 ;
    Phen(:,:,cur_n) = WT ;
    Phen_v(:,:,:,it) = Phen ;
    t_v(it,1) = t ;
    improvement(it,1) = fMT/fWT ;
end
if N_K>1
    plotPhenMultipleAntibiotics(fig_num, N_K, N, Phen_v,t_v)
else
    plotPhen(fig_num,N_P, N, Phen_v, t_v)
end
if save_fig
    figure(fig_num);
    savefig(['multiple_figs/simulated_N_p_' num2str(i_num_comp) '_N_r_' num2str(N_R) '_N_k_' num2str(N_K) '_eff_producers_' num2str(num_eff_producers) '.fig']);
end
% load chirp
% sound(y,Fs)
end