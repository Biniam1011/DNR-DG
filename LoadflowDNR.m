function [V_nodg,PL_nodg,QL_nodg] = LoadflowDNR()
global d V_base
from_bus = d(:,1);
to_bus = d(:,2);
r = d(:,3);
x = d(:,4);
z = complex(r,x);
no_branch = length(d(:,1));
no_bus = max(max(d(:,1:2)));
V_bus = V_base*ones(no_bus,1);
V_bus_max = 1.05*V_base;
V_bus_min = 0.95*V_base;
PD = zeros(no_bus,1);
QD = zeros(no_bus,1);
% PDG = zeros(no_bus,1);
S_base = 20000;
Z_base = V_base^2/S_base;
S = zeros(no_bus,1);
I_bus = zeros(no_bus,1);
for a = 1:no_branch
    b = to_bus(a);
    PD(b) = (d(a,5)*1e3);
    QD(b) = (d(a,6)*1e3);
    S(b) = complex(PD(b),QD(b));
end

tol = inf;
iter = 1;
while tol > 1e-10 && iter < 50
    
    for a = 1:no_branch
        b = to_bus(a);
        I_bus(b) = conj(S(b)/V_bus(b));
    end
    
    I_branch = zeros(no_branch,1);
    for a = no_branch:-1:1
        b = to_bus(a);
        for c = 1:no_branch
            if b == from_bus(c)
                I_branch(a) = I_branch(a) + I_branch(c);
            end
        end
        I_branch(a) = I_branch(a) + I_bus(b);
    end
    
    for a = 1:no_branch
        b = from_bus(a);
        c = to_bus(a);
        V_bus(c) = V_bus(b) - z(a)*I_branch(a);
        if V_bus(c) < 0.9*V_base
            V_bus(c) = V_bus(c) + 0.1*(V_bus_max-V_bus_min);
        else
        end
    end
    Scalc = V_bus.*conj(I_bus);
    Sflow  = 1e-3*(V_bus(1:end-1).*conj(I_branch(:,1:end)));
    Sloss  = conj(I_branch.^2.*z)*0.5e-3;
    tol = max(abs(real(Scalc) - real(S)));
    iter = iter + 1;
end
S_inj(1) = V_bus(1)*conj(I_branch(1));
P(1) = real(S_inj(1)); Q(1) = imag(S_inj(1));
Pkw = P;  Qkw = Q;
LP = sum(real(Sloss)); %Pkw(1)-sum(PD);
LQ = sum(imag(Sloss)); %Qkw(1)-sum(QD);
Ploss = LP; Qloss = LQ;
V = V_bus./V_base;
Vang = angle(V_bus);
V_nodg = abs(V)';
PL_nodg = real(Sloss); 
QL_nodg = imag(Sloss);

% % Loss sensitivity index
% % LSI1 = -2*((PD(2:end)+Ploss/S_base).^2+...
% %            (QD(2:end)+Qloss/S_base).^2).*(r*Z_base)...
% %              ./(V_nosvc(:,2:end).^3)';
% 
% LSI1 = -2*((PD(2:end)/S_base+Ploss).^2+...
%            (QD(2:end)/S_base+Qloss).^2).*(r)...
%              ./((V_nosvc(:,2:end).^3)'*V_base);
% [Value1,Index1]=sort(LSI1);
% Index1=Index1+1;
% LSI1=[Value1,Index1];
% 
% LSI2 = 2*((QD(2:end)/S_base+Qloss).^2).*(r.*Z_base)...
%         ./(V_base*(V_nosvc(:,2:end).^2)');
% [Value2,Index2]=sort(LSI2,'descend');
% Index2=Index2+1;
% LSI2=[Value2,Index2];



disp('***********************************')
disp('     With DNR Before DG placement  ')
disp('***********************************')
% disp(' ')
fprintf('Active Power loss wz DNR: %2.3fkW \n',Ploss)
fprintf('Reactive Power loss wz DNR: %2.3fkVAR \n',Qloss)

disp('-----------------------------------')
disp('       Bus Voltage Profile wz DNR  ')
disp('-----------------------------------')
disp('  Bus    Vbus(pu)  Angle(deg)')
disp('-----------------------------------')
for N = 1:no_bus
    fprintf(' %4.0f',N),
    fprintf('  %8.4f',V_nodg(N)),
    fprintf('  %8.4f\n',Vang(N)),
end

PL = real(Sloss); QL = abs(imag(Sloss));
% Ploss_tot = sum(PL); Qloss_tot = sum(QL);

disp('-------------------------------------------------------------')
disp('        Power Losses through each branch  wz DNR   ')
disp('-------------------------------------------------------------')
disp('Br  From   To  Pflow(kW)  Qflow(kVAR)  Ploss(kW)  Qloss(kVAR)')
disp('-------------------------------------------------------------')
for M = 1:no_branch
    fprintf(' %2.0f',M),
    fprintf(' %4.0f',from_bus(M)),
    fprintf(' %6.0f',to_bus(M)),
    fprintf('  %10.2f',real(Sflow(M))),
    fprintf('  %8.2f',imag(Sflow(M))),
    fprintf('  %8.2f',PL(M)),
    fprintf('  %8.2f\n',QL(M)),
end
disp('-------------------------------------------------------------')
fprintf('Total power loss                       ')
fprintf('%10.2f',Ploss)
fprintf('%11.2f\n',Qloss)
disp('-------------------------------------------------------------')

