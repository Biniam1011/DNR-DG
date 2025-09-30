function [VSI,Vm,PTloss,QTloss]=objective_fun(x,start)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sodonormal;

S_base=10000;   %(KVA)
V_base=15; %(kV)
Z_base=1000*(V_base^2)/(S_base);
I_base=S_base/(sqrt(3)*V_base);
if start==1
    for i=1:2
        bus_data(x(i),2)=bus_data(x(i),2)-409.3366  
    end
   
end

demanded_P=bus_data(:,2)/S_base;
demanded_Q=bus_data(:,3)/S_base;

R=line_data(:,3)/Z_base;
X=line_data(:,4)/Z_base;
nbus=length(bus_data);
nline=length(line_data);
%% calculating BIBC matrix of network
BIBC=zeros(nbus-1,nline);
for i=1:nline
    BIBC(:,line_data(i,2))=BIBC(:,line_data(i,1));
    BIBC(line_data(i,5),line_data(i,2))=1;
end
BIBC(:,1)=[]; %BIBC Matrix

%% Initilize voltage 
Vm=ones(1,nbus);   %voltage magnitude(pu)

%% start iteration
delta=1;eps=0.00001;iter=0;
MAXiter=1000;
while delta > eps
    iter=iter+1;
    if iter>MAXiter
        break;
    end
    
    for k=1:nbus
        Ibus(k,1)=(demanded_P(k)-sqrt(-1).*demanded_Q(k))./(conj(Vm(k)));
    end
    Inode=BIBC*Ibus(2:end); % It gives out branch current

    V(1)=1;
    for k=1:length(line_data)
       V(line_data(k,2))=Vm(line_data(k,1))-(R(k,1)+sqrt(-1)*X(k,1))*(Inode(line_data(k,5)));
    end
    delta=max(abs(V-Vm));
    Vm=V;
end
%% results 
if start==0
    figure(1)
    plot(abs(Vm))
    %,'--*b'
    ylabel('Voltage Magnitude (p.u.)')
    xlabel('Bus Number')
    title('Vlotage profile')
    grid on; grid minor;
end
%% loss calculation
    IL=Inode;
for k=1:length(line_data)
%     IL(k)=(V(line_data(k,1))-V(line_data(k,2)))/(R(k)+sqrt(-1).*X(k));
    S_line(k)=V(line_data(k,1))*conj(IL(k));
    PLoss(k)=R(k)*(abs(IL(k))^2)*S_base;
    QLoss(k)=X(k)*(abs(IL(k))^2)*S_base;
end

power_f_active=real(S_line)*S_base;
power_f_reactive=imag(S_line)*S_base;

PTloss=sum(PLoss); % total active or reactive loss
QTloss=sum(QLoss);% total reactive loss

%The strongest buss has high value of VSI
VSI= abs(Vm(1,2:end)').^4-4*(demanded_P(2:end,1).*X+demanded_Q(2:end,1).*R).^2-4*(demanded_P(2:end,1).*R+demanded_Q(2:end,1).*X).*abs(Vm(1,2:end)').^2;
%good case 
[Value3,Index3]=sort(VSI,'ascend');
Index3=Index3+1;
VSI=[Value3,Index3];




    