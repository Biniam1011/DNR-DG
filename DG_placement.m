
clc
clear
close all
format short g


%%%%%%%%%%% Base case powr flow results %%%%%%%%%%%%%%%%
start=0;
x=0;
 [VSI,Vm,PTloss,QTloss]=objective_fun(x,start);
disp(' ')
disp('======================================================')
disp('Results of 114-bus system  before DG placment ')
disp(' ')
disp(['Total active loss is: ' num2str(PTloss) ' MV'])
disp(['Total reactive loss is: ' num2str(QTloss) ' kvar'])
[value_v,index_v]=sort(abs(Vm));
disp(['Minimum voltage is: ' num2str(value_v(1)) ', at bus ' num2str(index_v(1))])
disp(['Maximum voltage is: ' num2str(value_v(end-1)) ', at bus ' num2str(index_v(end-1))])
disp(' ')
pause(0.5)

%%%%%%%%%%% Compansation %%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Case 1 %%%%%%%%%%%%%%%%%%%%%%%%
start=1;
%% parameters setting
a=2;
%for DG placment
%lb=[VSI(1:a,1)' 1*ones(1,a)]; % lower bound
%ub=[VSI(1:a,1)' 4*ones(1,a)];  % upper bound
lb=[VSI(1:a,2)' 1*ones(1,a)]; % lower bound
ub=[VSI(1:a,2)' 6*ones(1,a)];  % upper bound
nvar=2*a; % number of variable


NP=25;              % number particle
T=15;               % max of iteration

W=1;
C1=2;
C2=2;

alpha=0.05;

%% initialization
tic
empty.pos=[];
empty.cost=[];
empty.velocity=[];

particle=repmat(empty,NP,1);

for i=1:NP
particle(i).pos=lb+rand(1,nvar).*(ub-lb);

[particle(i).cost]=fitness(particle(i).pos,start);
particle(i).velocity=0;
end

bparticle=particle;

[value,index]=min([particle.cost]);
gparticle.cost=inf;
gparticle.pos=particle(index).pos;


%% main loop

best=zeros(T,1);
AVR=zeros(T,1);

for t=1:T

     for i=1:NP
         
          particle(i).velocity=W*particle(i).velocity...
                              +C1*rand(1,nvar).*(bparticle(i).pos-particle(i).pos)...
                              +C2*rand(1,nvar).*(gparticle.pos-particle(i).pos);
          
         particle(i).pos=particle(i).pos+particle(i).velocity;
         
         
        particle(i).pos=min(particle(i).pos,ub);
        particle(i).pos=max(particle(i).pos,lb);
          
          
          if sum(particle(i).pos(3:end))>=4
         [particle(i).cost]=fitness(particle(i).pos,start);
          else
              [particle(i).cost]=inf;
          end
         
         if particle(i).cost<bparticle(i).cost
             bparticle(i)=particle(i);
             
             if bparticle(i).cost<gparticle.cost
                 gparticle=bparticle(i);
             end
         end
         
         

     end
     
     
     
 W=W*(1-alpha);
 
 best(t)=gparticle.cost;
 AVR(t)=mean([particle.cost]);
 
 disp([ ' t = ' num2str(t)   ' BEST = '  num2str(best(t))]);
 

 
end


%% results
disp('====================================================')
disp([' Buses are   =  '  num2str(gparticle.pos(1:a))])
disp([' sizes are   =  '  num2str(gparticle.pos(a+1:end))])

j=0;
for i=1:a
    j=j+1;
    Buses(j)=gparticle.pos(i);
    Sizes(j)=gparticle.pos(a+i);
    
end
    disp(' ')
    disp('Buses    Number of DG points')
    disp([Buses' Sizes'])   
    disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VSI,Vm,PTloss,QTloss]=objective_fun(gparticle.pos,start);

figure(1)
hold on
plot(abs(Vm))
grid on; grid minor;

disp(' ')
disp('======================================================')
disp('Results of 114-bus system after DG placment')
disp(' ')
disp(['Total active loss is: ' num2str(PTloss) ' kW'])
disp(['Total reactive loss is: ' num2str(QTloss) ' kvar'])
[value_v,index_v]=sort(abs(Vm));
disp(['Minimum voltage is: ' num2str(value_v(1)) ', at bus ' num2str(index_v(1))])
disp(['Maximum voltage is: ' num2str(value_v(end-1)) ', at bus ' num2str(index_v(end-1))])
disp(' ')
legend('before DG placment','after DG placment')
grid on; grid minor;
figure(2)
bar([PTloss_DG PTloss; QTloss_DG QTloss])
legend('Active Loss After DG', 'Active Loss Before DG', 'Reactive Loss After DG', 'Reactive Loss Before DG')
xlabel('Power Loss Type')
ylabel('Power Loss (MW/kvar)')
title('Power Loss Comparison')
