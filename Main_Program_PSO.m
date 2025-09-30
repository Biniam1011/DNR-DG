
clear
close all
global d d2 V_base nbus
format shortg

%% Uncomment this section to run base case load flow
[d, V_base] = data();      % Read Systen data for base case
% nbus = max(max(d(:,1)));  % Number of buses
[V_nodg,PL_nodg,QL_nodg] = LoadflowBase(); % Load flow without DG

%% Uncomment this section to run DNR case load flow
% Read Systen data for DNR case
[d2, V_base] = DNR_data();
nbus = max(max(d(:,1)));
[V_dnr,PL_dnr,QL_dnr] = LoadflowDNR(); % Load flow without DG with DNR

%% This section for load flow with DNR and DG case
% DG limit in VARs
Qmin = -sum(d(:,6))*1e3;
Qmax =  sum(d(:,6))*1e3;

n = 1;           % Size of the swarm - no of birds
bird_setp = 20;  % Maximum number of birds steps
dim = 3;         % Dimension of the problem
% PSO parameter C1
c1min = 1.5;   c1max = 2.2;
% PSO parameter C2
c2min = 1.5;   c2max = 2.2;
% PSO momentum or inertia
wmax = 1.2;    wmin = 0.3;
R1 = rand(dim,n); R2 = rand(dim,n);
Rw = 1; %rand;
current_fitness = inf*ones(n,1);
current_position(1,:) =  Qmin + (Qmax - Qmin) * rand(1,n);
current_position(2,:) =  67;  %randi(nbus,1,n);
current_position(3,:) =  110; %randi(nbus,1,n);
velocity = .3*randn(dim,n) ;
local_best_position  = current_position;
for i = 1:n
    current_fitness(i) = Objective_func(current_position(:,i));
end
local_best_fitness  = current_fitness;
[global_best_fitness,g] = min(local_best_fitness);
for i=1:n
    globl_best_position(:,i) = local_best_position(:,g); %#ok
end
iter = 0;
itermax = bird_setp;
tic;
while  ( iter < bird_setp )
    iter = iter + 1;
    wtemp = wmax - ((wmax -wmin)/itermax)*iter;
    w = wmin + wtemp.*Rw;
    c1 = c1max - ((c1max -c1min)/itermax)*iter;
    c2 = c2max - ((c2max -c2min)/itermax)*iter;
    for i = 1:n
        current_fitness(i) = Objective_func(current_position(:,i));
    end
    for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
            local_best_fitness(i)  = current_fitness(i);
            local_best_position(:,i) = current_position(:,i);
        end
    end
    [current_global_best_fitness,g] = min(local_best_fitness);
    if current_global_best_fitness < global_best_fitness
        global_best_fitness = current_global_best_fitness;
        for i=1:n
            globl_best_position(:,i) = local_best_position(:,g);
        end
    end
    velocity = w *velocity+...
        c1*(R1.*(local_best_position-current_position))+...
        c2*(R2.*(globl_best_position-current_position));
    current_position = current_position + velocity;
    for j = 1:n
        current_position(:,j) = check(current_position(:,j));
    end
    fprintf(' Iter - % 2d , Objective function: %3.4f \n',...
              iter,0.5*global_best_fitness);
    global_best_fitness_iter(iter) = 0.5*global_best_fitness;
end
toc;
% PSO Convergence curve
figure(1);
plot(1:iter,global_best_fitness_iter,'k','LineWidth',2)
xlabel('Iteration Number')
ylabel('Real Power loss (kW)')
title('Convergence property of PSO algorithm')
grid on; grid minor;

%% Load flow with DG
X = globl_best_position(:,1);
[V_dg,PL_dg,QL_dg] = LoadflowDG(X);
nb = length(V_dg);
nbr = length(PL_dg);

%% Figures for Comparison
figure(2)
plot(1:nb,V_nodg,'b','LineWidth',2)
% legend('Base case','With DNR only')
xlabel('Buses')
ylabel('Vbus(pu)')
title('Bus voltage profile ')
grid on; grid minor;

figure(3)
plot(1:nb,V_nodg,'b',1:nb,V_dnr,'k','LineWidth',2)
legend('Base case','With DNR only')
xlabel('Buses')
ylabel('Vbus(pu)')
title('Bus voltage profile Comparison ')
grid on; grid minor;

figure(4)
plot(1:nb,V_nodg,'b',1:nb,V_dg,'k','LineWidth',2)
legend('Base case','With DG only')
xlabel('Buses')
ylabel('Vbus(pu)')
title('Bus voltage profile Comparison ')
grid on; grid minor;

figure(5)
plot(1:nb,V_nodg,'b',1:nb,V_dnr,'k',1:nb,V_dg,'r','LineWidth',2)
legend('Base case','With DNR only','With DNR and DG')
xlabel('Buses')
ylabel('Vbus(pu)')
title('Bus voltage profile Comparison ')
grid on; grid minor;

figure(6)
plot(1:nbr,PL_nodg,'b',1:nbr,PL_dnr,'k',1:nbr,PL_dg,'r','LineWidth',2)
legend('Base case','With DNR only','With DNR and DG')
xlabel('Branches')
ylabel('Ploss(kW)')
title('Active power loss Comparison with DNR and DG')
grid on; grid minor;

% figure(4)
% plot(1:nbr,QL_dnr,'k',1:nbr,QL_dg,'r','LineWidth',2)
% legend('With DNR only','With DNR and DG')
% xlabel('Branches')
% ylabel('Qloss(kVAR)')
% title('Reactive power loss Comparison with DNR and DG')
% grid on; grid minor;

%% Overall comparison
% Voltage comparison
[BB1,RR1] = min(V_nodg); [BB2,RR2]=min(V_dnr); [BB3,RR3]=min(V_dg);
fprintf('    Minimum bus voltage magnitudes     \n')
fprintf(' Base case            DNR case           DNR+DG case\n')
fprintf(' %2.4fpu at bus %2.0f',BB1,RR1),
fprintf('  %6.4f at bus %2.0f',BB2,RR2),
fprintf('  %6.4f at bus %2.0f\n',BB3,RR3),
disp('-------------------------------------------------------------')
% Power loss comparison
fprintf('    Total Active power loss     \n')
fprintf('  Base case     DNR case    DNR+DG case\n')
fprintf('  %2.4fpu ',sum(PL_nodg)),
fprintf('  %6.4f ',sum(PL_dnr)),
fprintf('   %6.4f \n',sum(PL_dg)),
disp('-------------------------------------------------------------')



