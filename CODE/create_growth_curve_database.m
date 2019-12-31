function [growth_object] = create_growth_curve_database(growth_rate_ratio,num_timepoints,n0,epsilon)
lambda1 = 1; 
lambdas_2 = growth_rate_ratio * lambda1 ;

find_steady_t = @(lambda, eps, n0) (1/lambda)*log((1+n0*eps-eps-n0)/(n0*eps));
t_final = find_steady_t(lambda1, epsilon, n0);
t_sample = linspace(0,t_final,num_timepoints);
growth_curves = nan(length(lambdas_2),2,length(t_sample));

for i_2 = 1:length(lambdas_2)
    lambda2 = lambdas_2(i_2);

    bacteria_growth = @(t,y) [lambda1*y(1)*(1-y(1)-y(2)); ...
                              lambda2*y(2)*(1-y(1)-y(2))];
    [t,y] = ode45(bacteria_growth,t_sample,[n0, n0]);
    y1 = y(:,1);
    y2 = y(:,2);
    growth_curves(i_2,1,:) = y1;
    growth_curves(i_2,2,:) = y2;
end

growth_object.growth_rate_ratio = growth_rate_ratio ;
growth_object.epsilon = epsilon;
growth_object.n0 = n0;
growth_object.num_timepoints = num_timepoints;
growth_object.growth_curves = growth_curves;
growth_object.t_sample = t_sample;
end