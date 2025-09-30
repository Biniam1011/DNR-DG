function [V_dg,PL_dg,QL_dg] = LoadflowDG(X)
global d V_base 
dg_q = X(1);
% dg_loc = round(X(2));
dg_loc1 = 67; %round(X(2));
dg_loc2 = 110; %round(X(3));
dg_loc = [dg_loc1 dg_loc2];
from_bus = d(:,1);
to_bus = d(:,2);
r = d(:,3);
x = d(:,4);
z = complex(r,x);
no_branch = length(d(:,1));
no_bus = max(max(d(:,1:2)));
V_bus = V_base*ones(no_bus,1);
PD = zeros(no_bus,1);
QD = zeros(no_bus,1);
DG_q = zeros(no_bus,1);
DG_q(dg_loc) = dg_q;
S = zeros(no_bus,1);
I_bus = zeros(no_bus,1);
for a = 1:no_branch
    b = to_bus(a);
    PD(b) = (d(a,5)*1e3);
    QD(b) = (d(a,6)*1e3) - DG_q(b);
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
    if V_bus(c) < 0.95*V_base
        V_bus(c) = V_bus(c) + 0.02*V_base;
    else
    end
end
    Scalc = V_bus.*conj(I_bus);
    Sflow  = 1e-3*(V_bus(1:end-1).*conj(I_branch(:,1:end)));
    Sloss  = conj(I_branch.^2.*z)*1e-3/2;
    tol = max(abs(real(Scalc) - real(S)));
    iter = iter + 1;
end
S_inj(1) = V_bus(1)*conj(I_branch(1));
P(1) = real(S_inj(1)); Q(1) = imag(S_inj(1));
Pkw = P; Qkw = Q;

LP = sum(abs(real(Sloss))); %Pkw(1)-sum(PD);
LQ = sum(abs(imag(Sloss))); %Qkw(1)-sum(QD);
% LP = sum(real(Sloss)); %Pkw(1)-sum(PD);
% LQ = sum(imag(Sloss)); %Pkw(1)-sum(PD);

Ploss = LP; Qloss = LQ;
PL_dg = abs(real(Sloss));
QL_dg = abs(imag(Sloss));
% PL_dg = real(Sloss); 
% QL_dg = imag(Sloss);

disp('***********************************')
disp('         DG placement             ')
disp('***********************************')
fprintf('DG location is Bus: %2.0f and %2.0f\n',dg_loc1,dg_loc2)
fprintf('DG Rating is: %2.3fkVAR \n',dg_q*1e-3)

disp('***********************************')
disp('         After DG placement')
disp('***********************************')
fprintf('Active Power loss is: %2.3fkW \n',Ploss)
fprintf('Reactive Power loss is: %2.3fkVAR \n',Qloss)

V = V_bus./V_base;
Vang = angle(V_bus);
V_dg = abs(V)';
disp('-----------------------------------')
disp('       Bus Voltage Profile         ')
disp('-----------------------------------')
disp('  Bus    Vbus(pu)  Angle(deg)      ')
disp('-----------------------------------')
for N = 1:no_bus
    fprintf(' %4.0f',N),
    fprintf('  %8.4f',V_dg(N)),
    fprintf('  %8.4f\n',Vang(N)),
end
PL = abs(real(Sloss)); QL = abs(imag(Sloss));
% Ploss_tot = sum(PL); Qloss_tot = sum(QL);

disp('-------------------------------------------------------------')
disp('        Power Losses through each branch     ')
disp('-------------------------------------------------------------')
disp('Br  From   To  Pflow(kW)  Qflow(kVAR)  Ploss(kW)  Qloss(kVAR)')
disp('-------------------------------------------------------------')
for M = 1:no_branch
    fprintf(' %2.0f',M),
    fprintf(' %4.0f',from_bus(M)),
    fprintf(' %6.0f',to_bus(M)),
    fprintf('  %10.2f',real(Sflow(M))),
    fprintf('  %8.2f',imag(Sflow(M))),
    fprintf('  %8.4f',PL(M)),
    fprintf('  %8.4f\n',QL(M)),
end
disp('-------------------------------------------------------------')
fprintf('Total power loss                       ')
fprintf('%10.4f',Ploss)
fprintf('%11.4f\n',Qloss)
disp('-------------------------------------------------------------')



