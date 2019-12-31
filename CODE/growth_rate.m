function g = growth_rate(P,cost)
g = exp(-sum(P(:,1)*cost(1) + P(:,2)*cost(2))) ;
end