
multiple_competiotors = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\producers.fig';
%multiple_competiotors = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\gamma_3\producers_new.fig';

average_competiotors  = 'C:\Users\Keren\Documents\MATLAB\AntibioticProduction\previous_r\producers_average.fig';
fig_path = {multiple_competiotors,average_competiotors};
dt = [0.01, 0.01];
color = [0 0 0; 0 0.450980392156863 0.741176470588235];
style ={'k','b'};
fig_handle = figure(); hold all
for  i= 1:length(fig_path)
    [data_struct,num_competitors] = drawDataFromFig(fig_path{i});
    [average_production, var_production] = calcAverageProduction(data_struct,num_competitors,dt(i));
    if i==1
        i_1 = find(num_competitors==1);
        figure();subplot(2,1,1);
        area(data_struct{i_1,1},data_struct{i_1,2},'facecolor','b','edgecolor','b','facealpha',0.5)
        xlim([0 0.5e6])
        ylim([0 5]);
        legend('1')
        ylabel('Production')
        set(gca,'fontsize',14);
        i_max = find(num_competitors==15);
        subplot(2,1,2); area(data_struct{i_max,1},data_struct{i_max,2},'facecolor','r','edgecolor','r','facealpha',0.5)
        xlim([0 0.5e6])
        ylim([0 5]);
        legend(num2str(num_competitors(i_max)))
        ylabel('Production')
        xlabel('Time [# mutations]')
        set(gca,'fontsize',14);
    end
    [~,I] = sort(num_competitors);
    figure(fig_handle); hold all
    errorbar(num_competitors(I),average_production(I),var_production(I),'color',color(i,:),'LineStyle','none','linewidth',2,'displayname','Average Production')
    
    [fitresult, gof] = fit_to_results(num_competitors, average_production);
    h(i) = plot(fitresult);
    h(i).LineWidth = 2;
    h(i).Color = color(i,:);
    if i ==1
        errorbar(num_competitors(i_1),average_production(i_1),var_production(i_1),'-b','linewidth',2,'displayname','Average Production')
        errorbar(num_competitors(i_max),average_production(i_max),var_production(i_max),'r','linewidth',2,'displayname','Average Production')
    end
end
legend([h(1), h(2)],'Original Model','Multiple Draws')
xlabel('Number of Competitors')
ylabel('Average Production')
set(gca,'fontsize',14);

%figure_file = 'C:\Users\Keren\Documents\MATLAB\Einat\temp_results\producers.fig';
%figure_file = 'C:\Users\Keren\Documents\MATLAB\Einat\no_competition_between_competitors\1_r_multiple_draws\t_correction\producers.fig';
%figure_file = 'C:\Users\Keren\Documents\MATLAB\Einat\no_competition_between_competitors\1_r_multiple_draws\producers_no_correction.fig';


%%
function [data_struct,num_competitors] = drawDataFromFig(figure_file)
    fig_handle = open(figure_file);
    axesObjs = get(fig_handle, 'Children');
    axesObjs = axesObjs(length(axesObjs));
    lineObjs = get(axesObjs, 'Children');
    
    data_struct = cell(length(lineObjs),2); % time, production
    num_competitors = nan(length(lineObjs),1);
    for i = 1:length(lineObjs)
        temp_t = lineObjs(i).XData(2:end);
        temp_p = lineObjs(i).YData(2:end);
        data_struct{i,1} = temp_t(~isnan(temp_p));
        data_struct{i,2} = temp_p(~isnan(temp_p));
        num_competitors(i) = str2double(lineObjs(i).DisplayName);
    end
    close(fig_handle);
end

function [average_production, var_production] = calcAverageProduction(data_struct,num_competitors,orig_dt)
    dt = orig_dt;
    average_production = nan(length(data_struct),1);
    var_production = nan(length(data_struct),1);

    for i = 1:length(data_struct)
        temp_t = data_struct{i,1}; %/num_competitors(i);
        if dt<1
            dt = ceil(orig_dt*max(temp_t));
        end
        temp_p = data_struct{i,2};
    
        I = find(diff(temp_t)==0);
        cum_production = cumtrapz(temp_t,temp_p);
        %cum_production = cumtrapz(temp_p(1:end-1).*diff(temp_t));
        %interped_cum_production = interp1(temp_t(1:end-1),cum_production,0:dt:max(temp_t));
        interped_cum_production = interp1(temp_t,cum_production,0:dt:max(temp_t));
        interped_cum_production(1) = 0;
        window_mean = diff(interped_cum_production)./dt;
        average_production(i) = mean(window_mean);
        var_production(i) = std(window_mean)/sqrt(length(window_mean));
    end
end
function [peak_height, t_peaks] = calcPeakHeightDuration(data_struct)
    p_minimal = 0.01;    
    num_lines = length(data_struct);
    peak_height = cell(num_lines,1);
    t_peaks = cell(num_lines,1);

    average_t_on = nan(num_lines,1);
    median_t_on = nan(num_lines,1);
    var_t_on = nan(num_lines,1);
    average_t_off = nan(num_lines,1);
    median_t_off = nan(num_lines,1);
    var_t_off = nan(num_lines,1);
    average_peak_height = nan(num_lines,1);
    median_peak_height = nan(num_lines,1);
    var_peak_height = nan(num_lines,1);

    for i = 1:num_lines
        t = data_struct{i,1};
        p = data_struct{i,2};
        on_off = p>p_minimal; %average_production;
        deriv_on_off = diff(on_off);
        i_on_to_off = find(deriv_on_off<0);
        i_off_to_on = find(deriv_on_off>0);
        if i_on_to_off(1)<i_off_to_on(1)
            i_off_to_on = [1, i_off_to_on];
        end

        on_to_off = t(i_on_to_off);
        off_to_on = t(i_off_to_on);
        trim_length = min(length(on_to_off),length(off_to_on));
        on_to_off = on_to_off(1:trim_length);
        off_to_on = off_to_on(1:trim_length);

        time_on = on_to_off - off_to_on;
        time_off = off_to_on(2:end) - on_to_off(1:length(off_to_on)-1);

        average_t_on(i) = mean(time_on);
        median_t_on(i) = median(time_on);
        var_t_on(i) = std(time_on)/sqrt(length(time_on)); 

        average_t_off(i) = mean(time_off);
        median_t_off(i) = median(time_off);
        var_t_off(i) = std(time_off)/sqrt(length(time_off)); 

        peaks = nan(length(off_to_on),1);
        for j = 1:length(off_to_on)
            peaks(j) = max(p(i_off_to_on(j):i_on_to_off(j)));
        end
        peak_height{i} = peaks;
        t_peaks{i} = time_on';
        
        average_peak_height(i) = mean(peaks);
        median_peak_height(i) = median(peaks);
        var_peak_height(i) = std(peaks)/sqrt(length(peaks));
    end
end
function [p_sampled, t_sampled] = calcProductionTime(data_struct) 
    num_lines = length(data_struct);
    p_sampled = cell(num_lines,1);
    t_sampled = cell(num_lines,1);
    for i = 1:num_lines
        t_sampled{i} = 0:10:max(data_struct{i,1});
        p_sampled{i} = interp1(data_struct{i,1},data_struct{i,2},t_sampled{i})';
    end
end
function [fitresult, gof] = fit_to_results(num_competitors, average_production)
    [xData, yData] = prepareCurveData( num_competitors, average_production );

    % Set up fittype and options.
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.19196847565222 0.454836407407937 0.719396528670862];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

