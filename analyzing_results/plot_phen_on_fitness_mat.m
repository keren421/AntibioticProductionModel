
fig_num = 10;
figure_file = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\N_p_1_N_r_1_N_k_1.fig'; %figure from single simulation
t_min = 0.2e4; % starting time of burst
t_max = 0.6e4; % end time of burst
window_lim = [1e4,3e4]; %cuts the simulation

fig_handle = open(figure_file);
axesObjs = get(fig_handle, 'Children');
axesObjs = axesObjs(length(axesObjs));
lineObjs = get(axesObjs, 'Children');

data_struct = cell(length(lineObjs),2); % time, production
num_competitors = nan(length(lineObjs),1);
for i = 1:length(lineObjs)
    t = lineObjs(i).XData(1:end);
    if strcmp(lineObjs(i).DisplayName,'Producer')
        producer = lineObjs(i).YData(1:end);
    elseif contains(lineObjs(i).DisplayName,'Competitor')
        resistant = lineObjs(i).YData(1:end);
    end
end
t = t - min(window_lim);
    close(fig_handle);
%%
I_zoom = find(t>t_min&t<t_max);
I_zero = find(producer(I_zoom)==0);
I_zero(diff(I_zero)==1)=[];
I_zoom =I_zoom(1:max(I_zero));

figure1 = figure('position',[0,0,800,300]);clf; hold on;
plot(t,producer,'-','Marker','.','MarkerSize',15,'linewidth',2,'displayname','Producer','Color',[0.929411764705882 0.690196078431373 0.129411764705882]')
plot(t,resistant,'-','Marker','.','MarkerSize',15,'linewidth',2,'displayname','Competitor','Color',[0.0 0.45 0.74]')
legend1 = legend(gca,'show');
set(legend1,...
    'Location','northwest',...
    'FontWeight','normal','linewidth',1,'fontsize',12);
box on;
xlim(window_lim-min(window_lim))
ylim([0,1.3*max(producer)]);
rectangle('Position',[t_min 0, t_max-t_min, 1.1*max(producer(I_zoom))],'linewidth',2)
xlabel('Time [# mutations]');
ylabel('Production/Resistance');

set(gca,'FontSize',12)
axes1 = axes('Parent',figure1,...
    'Position',[0.63 0.50 0.25 0.40]);
hold(axes1,'on');
plot(t(I_zoom),producer(I_zoom),'--','Marker','.','MarkerSize',15,'linewidth',2,'displayname','Producer')
plot(t(I_zoom),resistant(I_zoom),'-','Marker','.','MarkerSize',15,'linewidth',2,'displayname','Competitor','Color',[0.929411764705882 0.690196078431373 0.129411764705882])
ylim([0  1.1*max(producer(I_zoom))])
xlim([t_min , t_max-100])
box(axes1,'on');
set(axes1,'FontSize',12);

%%
load('fitness_map.mat') 
%%
color_map = [0.149019613862038 0.0901960805058479 0.400000005960464;0.154691487550735 0.112919472157955 0.434324681758881;0.160363361239433 0.135642856359482 0.468649327754974;0.166035234928131 0.158366248011589 0.502974033355713;0.171707108616829 0.181089639663696 0.537298679351807;0.177378982305527 0.203813031315804 0.5716233253479;0.183050841093063 0.226536422967911 0.605948030948639;0.188722714781761 0.249259814620018 0.640272676944733;0.194394588470459 0.271983206272125 0.674597322940826;0.200066462159157 0.294706612825394 0.708922028541565;0.205738335847855 0.31742998957634 0.743246674537659;0.20009545981884 0.343279868364334 0.75608366727829;0.194452583789825 0.369129747152328 0.768920660018921;0.188809707760811 0.394979625940323 0.781757593154907;0.183166831731796 0.420829504728317 0.794594585895538;0.177523955702782 0.446679383516312 0.807431578636169;0.171881079673767 0.472529292106628 0.820268571376801;0.166238188743591 0.498379170894623 0.833105564117432;0.160595312714577 0.524229049682617 0.845942556858063;0.154952436685562 0.550078928470612 0.858779549598694;0.149309560656548 0.575928807258606 0.87161648273468;0.143666684627533 0.6017786860466 0.884453475475311;0.138023808598518 0.627628564834595 0.897290468215942;0.133563980460167 0.639872372150421 0.860082387924194;0.129104167222977 0.652116179466248 0.822874248027802;0.124644346535206 0.664359986782074 0.785666167736053;0.120184525847435 0.6766037940979 0.748458027839661;0.115724705159664 0.688847601413727 0.711249947547913;0.111264884471893 0.701091408729553 0.674041867256165;0.106805063784122 0.71333521604538 0.636833727359772;0.102345243096352 0.725579023361206 0.599625647068024;0.0978854149580002 0.737822771072388 0.562417507171631;0.0934255942702293 0.750066578388214 0.525209426879883;0.0889657735824585 0.762310385704041 0.488001316785812;0.0845059528946877 0.774554193019867 0.450793206691742;0.0800461322069168 0.786798000335693 0.413585096597672;0.075586311519146 0.79904180765152 0.376377016305923;0.0711264908313751 0.811285614967346 0.339168906211853;0.0666666701436043 0.823529422283173 0.301960796117783;0.144444450736046 0.838235318660736 0.282352954149246;0.222222223877907 0.852941155433655 0.26274511218071;0.300000011920929 0.867647051811218 0.243137270212173;0.37777778506279 0.882352948188782 0.223529428243637;0.455555558204651 0.897058844566345 0.20392157137394;0.533333361148834 0.911764740943909 0.184313729405403;0.611111104488373 0.926470577716827 0.164705887436867;0.688888907432556 0.941176474094391 0.14509804546833;0.766666650772095 0.955882370471954 0.125490203499794;0.844444453716278 0.970588207244873 0.105882361531258;0.922222197055817 0.985294103622437 0.0862745121121407;1 1 0.0666666701436043;1 0.942810475826263 0.0758169963955879;1 0.885620951652527 0.0849673226475716;1 0.828431367874146 0.0941176488995552;1 0.771241843700409 0.103267982602119;1 0.714052319526672 0.112418308854103;1 0.656862735748291 0.121568635106087;1 0.599673211574554 0.13071896135807;1 0.542483687400818 0.139869287610054;0.940392136573792 0.488888919353485 0.129934638738632;0.880784332752228 0.435294151306152 0.120000004768372;0.821176469326019 0.381699353456497 0.11006536334753;0.761568665504456 0.328104585409164 0.100130721926689;0.701960802078247 0.274509817361832 0.0901960805058479];
%%
fig_num2 = figure('position',[0,0,800,300]); clf;
subplot(1,2,1);
contourf(production,resistance,cost_producer', 150,'EdgeColor','none');
view(2)
set(gca,'Ydir','Normal')
title('Producer Fitness');
xlim([1e-2 2.0]);
ylim([1e-2 1.0]);
h = colorbar;
ylabel(h, 'Fitness')
xlabel('Production');
ylabel('Resistance');
colormap(color_map);

subplot(1,2,2); 
contourf(production,resistance,cost_resistant', 150,'EdgeColor','none');
view(2)
set(gca,'Ydir','Normal')
title('Competitor Fitness');
xlim([1e-2 4.5]);
ylim([1e-2 4.5]);
h = colorbar;
ylabel(h, 'Fitness')
xlabel('Production');
ylabel('Resistance');
%% Add minimal resistance line and optimal production line
minimal_resistance = nan(length(production),1);
for i_p = 1:length(production)
    [~,I] = max(cost_resistant(i_p,:));
    if (~isempty(I))
        minimal_resistance(i_p) = resistance(I);
    end
end
optimal_production = nan(length(resistance),1);
for i_r = 1:length(resistance)
    [~,I] = max(cost_producer(:,i_r));
    if (~isempty(I))
        optimal_production(i_r) = production(I);
    end
end
[~,i_max_production] = max(optimal_production);

figure(fig_num2); hold all
for i = 1:2
    subplot(1,2,i); hold all;
    g(i) = plot(production, smooth(minimal_resistance),'-','color',[0.64 0.08 0.18],'linewidth',2,'displayname','minimal resistance');
    h(i) = plot(smooth(optimal_production(1:i_max_production)), resistance(1:i_max_production),'--','color', [0.64 0.08 0.18],'linewidth',2,'displayname','optimal production');
end

%%
for k = 1:2
    subplot(1,2,k);hold all
    set(gca,'fontsize',10)
    xlim([1e-2 2.0]);
    ylim([1e-2 1.0]);
    l_zoom = length(I_zoom);
    start_points = [max(producer(I_zoom(1:l_zoom-1)),0.01); resistant(I_zoom(1:l_zoom-1))];
    end_points = [max(producer(I_zoom(2:l_zoom)),0.01); resistant(I_zoom(2:l_zoom))];
    %end_points = [max(producer(circshift(I_zoom,-1)),0.01); resistant(circshift(I_zoom,-1))];
    arrow(start_points,end_points,8,90,40,2)  
    %(Start,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
    %if not using the arrow function:
    %hold on; plot(producer(I_zoom),resistant(I_zoom),'.-k',...
    %'linewidth',2,'markersize',12)
    legend([g(k) h(k)],{'minimal resistance','optimal production'})
end

