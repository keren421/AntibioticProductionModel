function [y1,y2] = single_droplet(P1,P2, g, type_resist, type_scoring, antibiotic_decay, time_limit, growth_curve)
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
        case 'keren_wants'
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
    
    %bacteria_growth = @(t,y) [g1*y(1)*(1-y(1)-y(2)); ...
    %                          g2*y(2)*(1-y(1)-y(2))];
    % try rewtitting as d[logy]/dt
    %[t,y] = ode45(bacteria_growth,[0 max_t],[intial_pop; intial_pop]);
    
    t = growth_curve(:,1);
    y = growth_curve(:,2:3);
    
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
            t_first = interp1(concentration,t,resistance,'linear',0);
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
    
    y1 = pop_size(1);
    y2 = pop_size(2);
end

