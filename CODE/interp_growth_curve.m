function [cur_growth] = interp_growth_curve(g_mutated,g_competitor, growth_curves)
    curves = growth_curves.growth_curves;
    g_ratios = growth_curves.growth_rate_ratio;
    t_sample = growth_curves.t_sample;
    g = [g_mutated,g_competitor];

    [fastest_growth, faster_growing] = max(g);
    growth_ratio = g(3-faster_growing)/g(faster_growing);
    num_curve = find(g_ratios>=growth_ratio,1,'first');
    
    if growth_ratio == g_ratios(num_curve)
        cur_growth = [t_sample'*(1/fastest_growth), ...
                      squeeze(curves(num_curve,faster_growing,:)), ...
                      squeeze(curves(num_curve,3 - faster_growing,:))];
    else
        faster_after = squeeze(curves(num_curve,faster_growing,:));
        slower_after = squeeze(curves(num_curve,3 - faster_growing,:));
        faster_before = squeeze(curves(num_curve - 1,faster_growing,:));
        slower_before = squeeze(curves(num_curve - 1,3 - faster_growing,:));

        faster_interp = (faster_before*(g_ratios(num_curve)-growth_ratio) + ...
                        faster_after*(growth_ratio-g_ratios(num_curve-1)))/(g_ratios(num_curve) - g_ratios(num_curve-1));

        slower_interp = (slower_before*(g_ratios(num_curve)-growth_ratio) + ...
                        slower_after*(growth_ratio-g_ratios(num_curve-1)))/(g_ratios(num_curve) - g_ratios(num_curve-1));

        cur_growth = [t_sample'*(1/fastest_growth), ...
                      faster_interp, ...
                      slower_interp];
    end
end

