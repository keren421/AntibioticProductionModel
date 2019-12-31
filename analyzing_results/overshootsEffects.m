fig_num = 10;
dt = 10;

single_r_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\N_p_1_N_r_1_N_k_1.fig';
no_overshoots_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\N_p_1_N_r_1_no_overshoots.fig';

%single_r_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\multiple_producers\N_p_1_N_r_1_N_k_1_no_overshots.fig';
%single_r_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\no_competition_between_competitors\same_mutation_rate_loser_remains\time_correction\1_p_1_r.fig';

%close all;
[t,producer, resistant] = draw_p_r_from_fig(single_r_file);
[t_overshoots,producer_overshoots, resistant_overshoots] = draw_p_r_from_fig(no_overshoots_file);
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
lose_vs_gain_production = sign(producer(2:end)-producer(1:end-1));
went_up =true;
going_down= [];
going_up= [];
for i =1:length(producer)-1
    if went_up & (producer(i)-producer(i+1))>0 && producer(i)>=0.0
        going_down = [going_down i];
        went_up = false;
    end
    if (producer(i)-producer(i+1))<0 && producer(i)>=0.0
        going_up = [going_up i];
        went_up = true;
    end
end
minimal_line_t = interp1(production, minimal_resistance, producer(1:end-1));
distance_from_minimal = (resistant(1:end-1)-minimal_line_t);
figure(); plot(producer)
hold all; plot(going_down,producer(going_down),'.r')
hold all; plot(going_up,producer(going_up),'.g')
hold all; plot(minimal_line_t,'k')
hold all; plot(resistant,'r')

fig_handle = figure(); yyaxis left; hold all;
bins = linspace(-0.2,0.4, 50);
a = histogram(distance_from_minimal([going_up going_down]),bins,'displayname','All instances');
b = histogram(distance_from_minimal(going_down), bins,'displayname','Lose Production');
%plot([0 0],[0 1.2*max(a.Values)],'-','linewidth',2,'color','k')
ylim([0 1.2*max(a.Values)]);

xlabel('Distance from Minimal Resistance');
ylabel('Histogram');
set(gca,'fontsize',12)
box on;
%%
bins = linspace(-0.2,0.4, 50);
up = nan(1,length(bins)-1);
down = nan(1,length(bins)-1);

minimal_line_t = interp1(production, minimal_resistance, producer(1:end-1));
distance_from_minimal = (resistant(1:end-1)-minimal_line_t);
lose_vs_gain_production = sign(producer(2:end)-producer(1:end-1));

for i =1:length(bins)-1
    up(i) = sum(lose_vs_gain_production(distance_from_minimal>bins(i)&distance_from_minimal<bins(i+1))>0);
    down(i) = sum(lose_vs_gain_production(distance_from_minimal>bins(i)&distance_from_minimal<bins(i+1))<0);
end
figure(fig_handle); hold all; yyaxis right;
%a=area(0.5*(bins(1:end-1)+bins(2:end)),up./(up+down),'FaceAlpha',0.6);
b=plot(0.5*(bins(1:end-1)+bins(2:end)),smooth(down./(up+down)),'displayname','Probability to Lose Production','linewidth',2);
legend([a b],{'Gain Production','Lose Production'})
plot([0 0],[0 1.0],'-','linewidth',2,'color','k')
legend([a b],{'Gain Production','Lose Production'})
%xlabel('Distance from Minimal Resistance');
ylabel('Probability to Lose Production')
ylim([0 1.1]);
set(gca,'fontsize',12)
box on;
%%
dt = 1e4;
[average_p, error_p] = find_avereage_p(producer,t,dt);
[average_p_overshoots, error_p_overshoots] = find_avereage_p(producer_overshoots,t_overshoots,dt);
X = categorical({'No overshots','Actual Data'});
X = reordercats(X,{'No overshots','Actual Data'});
Y = [average_p_overshoots average_p];

figure(); bar(X,Y)
ylabel('Average Production Rate')
hold on;
errorbar(X ,Y, [error_p, error_p_overshoots],'linewidth',2,'displayname','average production'); %,'color',CM(j,:))
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

function [average_p, error_p] = find_avereage_p(production,time,dt)
    cum_production = cumtrapz(time,production);
    interped_cum_production = interp1(time,cum_production,0:dt:max(time));
    interped_cum_production(1) = 0;
    window_mean = diff(interped_cum_production)./dt;
    average_p = mean(window_mean);
    error_p = std(window_mean)/sqrt(length(window_mean));
end