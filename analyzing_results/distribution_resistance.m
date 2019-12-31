fig_num = 10;
dt = 10;
multiple_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\N_p_1_N_r_15_N_k_1.fig';
close all;
[t,producer, resistant] = draw_p_r_from_fig(multiple_file);
[~,producer, resistant] = resample_at_set_times(t, producer, resistant,dt);
multiple_draws_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\N_p_1_N_r_15_average.fig';
[t,producer_draws, resistant_draws] = draw_p_r_from_fig(multiple_draws_file);
% [~,producer_draws, resistant_draws] = resample_at_set_times(t, producer_draws, resistant_draws,dt);

%%
load('C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\fitness_map.mat')
minimal_resistance = nan(length(production),1);
for i_p = 1:length(production)
    [~,I] = max(cost_resistant(i_p,:));
    if (~isempty(I))
        minimal_resistance(i_p) = resistance(I);
    end
end
%%
p_bins = linspace(0,max(producer),40);

minimal_resistance_bins = interp1(production, minimal_resistance, p_bins);
minimal_resistance_mid_bin = interp1(production, minimal_resistance, 0.5*(p_bins(2:end)+p_bins(1:end-1)));

for i = 1:length(p_bins)-1
    going_up = producer(2:end)-producer(1:end-1)>=0.01;
    going_up = [going_up 0];
    i_in_bin = find(producer>p_bins(i)&producer<p_bins(i+1)&going_up);
    minimal_resistance_in_bin = interp1(production, minimal_resistance, producer(i_in_bin))';
    resistance_in_bin = resistant(i_in_bin,:) - minimal_resistance_in_bin;
    
    %X=randperm(numel(resistance_in_bin));
    %ShuffledData=reshape(resistance_in_bin(X),size(resistance_in_bin));
    %mean_resistance = mean(ShuffledData,2); %
    mean_resistance = mean(resistance_in_bin,2);
    if length(mean_resistance)<1e2
        continue
    end
%     i_in_bin_single = find(producer_single>p_bins(i)&producer_single<p_bins(i+1));
%     resistance_single_in_bin = resistant_single(i_in_bin_single);
    i_in_bin_draws = find(producer_draws>p_bins(i)&producer_draws<p_bins(i+1));
    resistance_draws_in_bin = resistant_draws(i_in_bin_draws);
    
    figure(); hold on
    r_bins = linspace(0.8*min(min(resistance_in_bin)), 0.8*max(max(resistance_in_bin)),50);
    a = histogram(mean_resistance,r_bins, 'Normalization', 'probability','displayname','Average Resistance');

    b=histogram(resistance_in_bin(:),r_bins, 'Normalization', 'probability','displayname','All Resistance');
%     histogram(resistance_single_in_bin(:),r_bins, 'Normalization', 'probability','displayname','Single Competitor');
%    histogram(resistance_draws_in_bin(:),r_bins, 'Normalization', 'probability','displayname','Multiple draws');
%     histogram(mean_resistance_ordered,r_bins, 'Normalization', 'probability','displayname','Average Resistance - ordered');
%     histogram(resistance_in_bin_ordered(:),r_bins, 'Normalization', 'probability','displayname','All Resistance - ordered');
legend([a,b],{'Average Resistance', 'All Resistance'});
    max_prob = max(max(a.Values), max(b.Values));
    plot([0, 0]+ minimal_resistance_bins(i) - minimal_resistance_mid_bin(i) , [0 1.2*max_prob], 'k','linewidth',3);
    plot([0, 0]+ minimal_resistance_bins(i+1) - minimal_resistance_mid_bin(i) , [0 1.2*max_prob], 'k','linewidth',3);
    ylim([0, 1.2*max_prob]);
    pd = fitdist(mean_resistance,'Normal');
    text(r_bins(3),0.8*max_prob,['\sigma_{average}/\sigma_{all} = ' num2str(var(resistance_in_bin(:))/var(mean_resistance))])
    pd = fitdist(resistance_in_bin(:),'Normal');
    title(['Production between ' num2str(p_bins(i)) ' to '  num2str(p_bins(i+1))])
end
%%


function [t, producer, resistant] = draw_p_r_from_fig(figure_file)
    fig_handle = open(figure_file);
    axesObjs = get(fig_handle, 'Children');
    axesObjs = axesObjs(length(axesObjs));
    lineObjs = get(axesObjs, 'Children');

    num_competitors = nan(length(lineObjs),1);

    t = lineObjs(1).XData(1:end);
    producer = nan(length(t),1);
    resistant = nan(length(t),length(lineObjs)-1);
    num_r = 1;
    for i = 1:length(lineObjs)

        if strcmp(lineObjs(i).DisplayName,'Producer')
            producer(:,1) = lineObjs(i).YData(1:end);
        elseif contains(lineObjs(i).DisplayName,'Competitor')
            resistant(:,num_r) = lineObjs(i).YData(1:end);
            num_r = num_r+1;
        end
    end
    close(fig_handle);
end
function [t_sampled, producer_sampled, resistant_sampled] = resample_at_set_times(t, producer, resistant,dt)
    t_sampled = min(t):dt:max(t);    
    producer_sampled = interp1(t,producer, t_sampled);
    resistant_sampled = interp1(t,resistant, t_sampled);
end