a = 2;
b = 0.3;
num_samples = 1e6;
r1 = gamrnd(a,b,num_samples,1); %random(pd,100000,1);
peak = (a-1)*b;
r1 = peak - r1;
figure(); histogram(r1,linspace(-2,1,100),'Normalization','pdf')
title('Mutation Distribution');
xlabel('Mutation Size');
ylabel('pdf')
set(gca,'fontsize',14)