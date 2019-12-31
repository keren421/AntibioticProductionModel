Cost = [0.005, 0.01];
time_limit = 0.8;
decay = 1;
scoring_type = 'loser_remains_winner_gets_rest'; % different scoring methods 'winner_gets_all', 'loser_dies_winner_gets_rest', 'loser_remains_winner_gets_rest'
with_solver = 0; % 1- to use solver for each point, 0 - to use the saved profiles

production = 0.0:0.005:10; %logspace(log(0.01), log10(4), 400); %
resistance = 0.0:0.005:10; %logspace(log(0.01), log10(2), 400); %0.0:0.01:2;

if with_solver == 0
    clear growth_curves
    growth_curve_filename = ['CODE/saved_profiles' filesep 'growth_curves_log100.mat'] ;%growth curves database
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
end
save_file = 0;
map_name = 'loser_remains_time_limit_10'; % name of file if saving
fig_num = 1200;
%% calc fitness for each point
producer = zeros(1,2);
resistant = zeros(1,2);
cost_producer = nan(length(production),length(resistance));
cost_resistant = nan(length(production),length(resistance));

for i_p = 1:length(production)
    for i_r = 1:length(resistance)
        producer(2) = production(i_p);
        resistant(1) = resistance(i_r);
        
        g1 = growth_rate(producer, Cost);
        g2 = growth_rate(resistant, Cost);
        g = [g1,g2];
        
        if with_solver == 1
            [y_p,y_r,~,~] = single_droplet_with_solver(producer,resistant,g,'max', scoring_type, decay, time_limit);
        elseif with_solver == 0
            g1 = growth_rate(producer, Cost);
            g2 = growth_rate(resistant, Cost);
            cur_growth = interp_growth_curve(g1,g2, growth_curves);
            [y_p,y_r] = single_droplet(producer,resistant, g,'max', scoring_type, decay, time_limit, cur_growth);
        end
        cost_producer(i_p,i_r) = y_p;
        cost_resistant(i_p,i_r) = y_r;
    end
end
fitness_producer = cost_producer;
fitness_resistant= cost_resistant;

%% Draw maps of fitness
figure(fig_num); clf;

subplot(1,2,1);
imagesc(production, resistance, fitness_producer');
set(gca,'Ydir','Normal')
title('Producer Fitness');
colorbar;
xlabel('Production');
ylabel('Resistance');

subplot(1,2,2); 
imagesc(production, resistance, fitness_resistant');
set(gca,'Ydir','Normal')
title('Resistant Fitness');
colorbar;
xlabel('Production');
ylabel('Resistance');
%% Add minimal resistance line
minimal_resistance = nan(length(production),1);
for i_p = 1:length(production)
    [~,I] = max(cost_resistant(i_p,:));
    if (~isempty(I))
        minimal_resistance(i_p) = resistance(I);
    end
end
figure(); plot(production, minimal_resistance,'linewidth',2)
title('Minimal Resistance for Producing Value')
xlabel('Production');
ylabel('Minimal Resistance');
set(gca,'fontsize',14);
grid on

figure(fig_num);
for i = 1:2
    subplot(1,2,i); hold all;
    plot(production, minimal_resistance,'r','linewidth',2,'displayname','minimal resistance')
end
%% Add optimal production line
optimal_production = nan(length(resistance),1);
for i_r = 1:length(resistance)
    [~,I] = max(cost_producer(:,i_r));
    if (~isempty(I))
        optimal_production(i_r) = production(I);
    end
end
figure(); plot(optimal_production,resistance,'linewidth',2)
title('Minimal Resistance for Producing Value')
xlabel('Production');
ylabel('Minimal Resistance');
set(gca,'fontsize',14);
grid on
[~,i_max_production] = max(optimal_production);
figure(fig_num);
for i = 1:2
    subplot(1,2,i); hold all;
    plot(optimal_production(1:i_max_production), resistance(1:i_max_production),'k','linewidth',2,'displayname','optimal production')
end

%% Add optimal production line
loss_of_function = nan(length(production),1);

for i_p = 1:length(production)
    sensetive = (resistance<minimal_resistance(i_p));
    for i_r = sum(sensetive):-1:1
        if cost_producer(i_p,i_r)<cost_producer(1,i_r)
            loss_of_function(i_p) = resistance(i_r);
        end
    end
end
figure(); plot(production,loss_of_function,'.','linewidth',2)
title('Neutral Production line')
xlabel('Production');
ylabel('Minimal Resistance');
set(gca,'fontsize',14);
grid on
[~,i_max_production] = max(optimal_production);
figure(fig_num);
for i = 1:2
    subplot(1,2,i); hold all;
    plot(production, loss_of_function,'k','linewidth',2,'displayname','Neutral Production')
end
%% Draw lines instead of map
sample_points = 1:100:1000;
figure(); hold all
CM = jet(length(sample_points));
for i = 1:length(sample_points)
    
    cur_point = sample_points(i);
    cur_r = resistance(cur_point);
    plot(production, cost_producer(:,cur_point),'linewidth',2,...
        'color',CM(i,:),'displayname',['R = ' num2str(cur_r)]);
    slope = gradient(cost_producer(:,cur_point),production);

end
grid on;
box on
legend show
set(gca,'fontsize',14)
xlabel('Production')
ylabel('Score');
%% Draw lines instead of map
sample_points = 51:50:275;
figure(); hold all
CM = jet(length(sample_points));
for i = 1:length(sample_points)
    
    cur_point = sample_points(i);
    cur_p = production(cur_point);
    plot(resistance, cost_producer(cur_point,:),'linewidth',2,...
        'color',CM(i,:),'displayname',['P = ' num2str(cur_p)]);
    slope = gradient(cost_producer(:,cur_point),production);

end
grid on;
box on
legend show
set(gca,'fontsize',14)
title('Producer''s Fitness')
xlabel('Resistance')
ylabel('Fitness Score');
%% Draw lines instead of map
sample_points = 1:50:275;
figure(); hold all
CM = jet(length(sample_points));
for i = 1:length(sample_points)
    
    cur_point = sample_points(i);
    cur_p = resistance(cur_point);
    plot(production, cost_resistant(cur_point,:),'linewidth',2,...
        'color',CM(i,:),'displayname',['R = ' num2str(cur_p)]);
    slope = gradient(cost_producer(:,cur_point),production);

end
grid on;
box on
legend show
set(gca,'fontsize',14)
title('Producer''s Fitness')
xlabel('Resistance')
ylabel('Fitness Score');
%% save file
if save_file 
    folder_path = ['C:\Users\Keren\Documents\MATLAB\AntibioticProduction\fitness_maps\examples\' map_name '\'];
    mkdir(folder_path)
    savefig([folder_path 'Fitness_matrix.fig'])
    saveas(gcf,[folder_path 'Fitness_matrix.jpg'])
    save([folder_path 'vars.mat'],'Cost','run_name', 'scoring_type','decay', 'production','resistance','cost*','fitness*') 
end