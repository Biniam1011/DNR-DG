function LP = Objective_func(X)
global d V_base
svc_q = X(1);
svc_loc1 = round(X(2));
svc_loc2 = round(X(3));
svc_loc = [svc_loc1 svc_loc2];
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
SVC_q = zeros(no_bus,1);
SVC_q(svc_loc) = svc_q;
S = zeros(no_bus,1);
I_bus = zeros(no_bus,1);
for a = 1:no_branch
    b = to_bus(a);
    PD(b) = (d(a,5)*1e3);
    QD(b) = (d(a,6)*1e3) - SVC_q(b);
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
    Sloss  = conj(I_branch.^2.*z)*1e-3;
    tol = max(abs(real(Scalc) - real(S)));
    iter = iter + 1;
    
end
S_inj(1) = V_bus(1)*conj(I_branch(1));
P(1) = real(S_inj(1));
Q(1) = imag(S_inj(1));
Pkw = P;
% LP = sum(real(Sloss)); %Pkw(1)-sum(PD);
LP = sum(abs(real(Sloss))); %Pkw(1)-sum(PD);
% LP = Pkw(1)-sum(PD);
% LQ = sum(imag(Sloss)); %Qkw(1)-sum(QD);

