function [y1,y2,production_resistance_1,production_resistance_2] = single_droplet_with_solver(P1,P2, g, type_resist, type_scoring, antibiotic_decay, time_limit)
    eps = 1e-6; %integral is computed until population 1-eps 
    intial_pop = 1e-3;
    
    g1 = g(1) ;
    g2 = g(2) ;
    
    if antibiotic_decay
        production_resistance_1 = P1(:,2);
        production_resistance_2 = P2(:,2);
    else
        production_resistance_1 = 1/g1 * P1(:,2)*log(1/(2*eps));
        production_resistance_2 = 1/g2 * P2(:,2)*log(1/(2*eps));
    end
            
    switch type_resist
        case 'max'
            R1 = max(P1(:,1),production_resistance_1);
            R2 = max(P2(:,1),production_resistance_2);
        case 'plus'
            R1 = P1(:,1) + production_resistance_1;
            R2 = P2(:,1) + production_resistance_2;
        case 'resist'
            R1 = P1(:,1);
            R2 = P2(:,1);
    end
                    
    num_antibiotics = length(P1(:,1));
    
    if isnan(time_limit)
        find_steady_t = @(lambda, eps, n0) (1/lambda)*log((1+n0*eps-eps-n0)/(n0*eps));
        t1 = find_steady_t(g1, eps, intial_pop);
        t2 = find_steady_t(g2, eps, intial_pop);
        max_t = min(t1,t2) ;
    else
        delta = 1- time_limit;
        max_t = log((1+intial_pop*delta-delta-intial_pop)/(intial_pop*delta));
    end
    
    bacteria_growth = @(t,y) [g1*y(1)*(1-y(1)-y(2)); ...
                             g2*y(2)*(1-y(1)-y(2))];
    %try rewtitting as d[logy]/dt
    %[t,y] = ode45(bacteria_growth,[0 max_t],[intial_pop; intial_pop]);
    [t,y] = ode45(bacteria_growth,linspace(0,max_t,100),[intial_pop; intial_pop]);
    %t = growth_curve(:,1);
    %y = growth_curve(:,2:3);
    
    t_death = max_t;
    overall_losing_bacteria = nan;
    for i = 1:num_antibiotics
        if antibiotic_decay
            concentration = y(:,1)* P1(i,2) + y(:,2)*P2(i,2);
        else
            concentration = cumtrapz(t,y(:,1)* P1(i,2) + y(:,2)*P2(i,2));
        end
        [resistance, most_sensitive] = min([R1(i),R2(i)]);
        %I = find(concentration>resistance,1,'first');
        if concentration(end)>resistance
            t_first = interp1(concentration,t,resistance,'pchip',0);
            if t_first< t_death %~isempty(I) && t(I)<t_death
                t_death = t_first; %t(I);
                overall_losing_bacteria = most_sensitive;
            end
        end
    end
    
    pop_size(1) = interp1(t,y(:,1),t_death); 
    pop_size(2) = interp1(t,y(:,2),t_death); 
    
    if ~isnan(overall_losing_bacteria)
        overall_winning_bacteria = 3 - overall_losing_bacteria ;
        switch type_scoring
            case 'winner_gets_all'
                if isnan(time_limit)
                    pop_size(overall_winning_bacteria) = 1;
                else
                    g = [g1, g2];
                    g = g(overall_winning_bacteria);
                    start_pop = pop_size(overall_winning_bacteria);
                    pop_size(overall_winning_bacteria) = 1/(1+((1-start_pop)/start_pop)*exp(-g*(max_t-t_death)));
                end
                pop_size(overall_losing_bacteria) = 0;
            case 'loser_dies_winner_gets_rest'
                if isnan(time_limit)
                    pop_size(overall_winning_bacteria) = 1 - pop_size(overall_losing_bacteria);
                else
                    g = [g1, g2];
                    g = g(overall_winning_bacteria);
                    start_pop = pop_size(overall_winning_bacteria);
                    K = 1 - pop_size(overall_losing_bacteria);
                    pop_size(overall_winning_bacteria) = K/(1+((K-start_pop)/start_pop)*exp(-g*(max_t-t_death)));
                end
                pop_size(overall_losing_bacteria) = 0;
            case 'loser_remains_winner_gets_rest'
                if isnan(time_limit)
                    pop_size(overall_winning_bacteria) = 1 - pop_size(overall_losing_bacteria);
                else
                    g = [g1, g2];
                    g = g(overall_winning_bacteria);
                    start_pop = pop_size(overall_winning_bacteria);
                    K = 1 - pop_size(overall_losing_bacteria);
                    pop_size(overall_winning_bacteria) = K/(1+((K-start_pop)/start_pop)*exp(-g*(max_t-t_death)));
                end   
        end
    end
    plot_sim = 1;
    if plot_sim
        figure();
        
        if t_death<max_t
            winner_growth = @(t,y)g(overall_winning_bacteria)*y*(1-y);
            [t_after_death,y_after_death] = ode45(winner_growth,linspace(t_death,max_t,100),start_pop);
            t_winner = [t(t<t_death); t_after_death];
            pop_winner = [y(t<t_death,overall_winning_bacteria); y_after_death];
            t_loser = [t(t<t_death); t_after_death];
            pop_loser = [y(t<t_death,overall_losing_bacteria); pop_size(overall_losing_bacteria)*ones(size(t_after_death))];
        else
            t_winner = t;
            pop_winner = y(:,1);
            t_loser = t;
            pop_loser = y(:,2);
        end
        
        subplot(2,1,1); hold all;
        plot(t_winner,pop_winner,'-','linewidth',2,'color',[0.93 0.69 0.13])
        plot(t_loser,pop_loser,'-','linewidth',2,'color',[0 0.45 0.74])
        plot([max_t max_t],[0 10],'--','linewidth',1,'color','k')
        text(max_t,-0.05,'t_{final}','EdgeColor','none');
        if t_death<max_t
            plot(t(t>t_death),y(t>t_death,overall_winning_bacteria),':','linewidth',2,'color',[0.93 0.69 0.13])
            plot(t(t>t_death),y(t>t_death,overall_losing_bacteria),':','linewidth',2,'color',[0 0.45 0.74])
            plot(t(end),pop_winner(end),'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 0],'MarkerFaceColor',[0.93 0.69 0.13])
            plot(t(end),pop_loser(end),'MarkerSize',10,'Marker','o','LineStyle','none','Color',[0 0 0],'MarkerFaceColor',[0 0.45 0.74])
            plot([t_death t_death],[0 1],'--','linewidth',1,'color','k')
            ylim([0 max(max(y))])
            text(1.0*t_death,-0.05,'t_{kill}','EdgeColor','none');
        end
        xlabel('t')
        ylabel('Population Size')
        legend('Producer','Competitor','Location','northwest')
        set(gca,'fontsize',12)
        box on
        
        subplot(2,1,2); hold all;
        if overall_winning_bacteria ==2 
            concentration = pop_loser*P1(1,2) + pop_winner*P2(1,2);
        else
            concentration = pop_winner* P1(1,2) + pop_loser*P2(1,2);
        end
        plot(t_winner,concentration,'linewidth',2,'color',[0.93 0.69 0.13])
        plot([t_death t_death],[0 10],'--','linewidth',1,'color','k')
        text(1.00*t_death,-0.05,'t_{kill}','EdgeColor','none');
        plot([max_t max_t],[0 10],'--','linewidth',1,'color','k')
        text(max_t,-0.05,'t_{final}','EdgeColor','none');
        plot([0 max_t],[resistance resistance],'--','linewidth',1,'color','k')
        text(0.2,resistance+0.05,'Antibiotic Resistance','EdgeColor','none');
        ylim([0 max(concentration)])
        xlabel('t')
        ylabel('Antibiotic Concentration')
        set(gca,'fontsize',12)
        box on
        
        y1 = pop_size(1);
        y2 = pop_size(2);
    end
end

