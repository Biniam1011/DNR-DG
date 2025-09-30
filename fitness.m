function Z=fitness(x,start)
[VSI,V_sys,PTloss,QTloss]=objective_fun(x,start);
Z=PTloss;
end