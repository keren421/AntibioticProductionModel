function [] = plotPhen(fig_num, N_p, N, Phen_v,t_v)
    figure(fig_num);clf;
    hold on
    for k = 1:N
        if k<=N_p
            plot(t_v,squeeze(Phen_v(1,2,k,:)),'--','displayname','Producer')
        else
            plot(t_v,squeeze(Phen_v(1,1,k,:)),'-','displayname',['Competitor ' num2str(k-N_p)])
        end
        set(gca,'fontsize',14)
    end
    drawnow
end