function X = check(X)
global d nbus
SVC = X(1);
% svc_loc = round(X(2));
svc_loc1 = 67; %round(X(2));
svc_loc2 = 110; %round(X(3));
svc_loc = [svc_loc1 svc_loc2];
Pmin = 0;
Pmax = sum(d(:,5))*1e3;
if SVC < Pmin || SVC > Pmax
    SVC = Pmin + (Pmax - Pmin) * rand;
end
if svc_loc < 1 | svc_loc > nbus
    svc_loc = randi(nbus);
end
X = [SVC svc_loc];