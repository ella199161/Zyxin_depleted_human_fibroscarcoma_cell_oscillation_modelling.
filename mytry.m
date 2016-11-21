clear
L=50;

MT=2;
Vm=5;
dt=0.001;
T=300;
nT=T/dt;
Ctotal=10;
Dx(1,1)=0.01;

for i=1:1:10
Dx(i,1)=0.01;
end

for i=10:1:15 
Dx(i,1)=0.02;
end
 for i=15:1:19 
Dx(i,1)=0.05;
end
 for i=19:1:22 
Dx(i,1)=0.2;
end
 for i=22:1:118
Dx(i,1)=0.5;
end
 for i=118:1:121
Dx(i,1)=0.2;
end
 for i=121:1:125
Dx(i,1)=0.05;
end
 for i=125:1:130 
Dx(i,1)=0.02;
end
 for i=130:1:140 
Dx(i,1)=0.01;
 end
ndg = length(Dx);
 N = length(Dx);
 Ngrid = N+1;
X = zeros(Ngrid,1);
for in = 1:1:N
    X(in+1) = X(in)+Dx(in);
end
C = zeros(Ngrid,1);
ns = ndg/5 - mod(ndg/5,1);
C(2:ns+1) = Ctotal / (sum(Dx(2:ns))+Dx(ns+1)/2+Dx(1)/2);
dCdt = zeros(Ngrid,1);

subplot(2,2,1);
plot(X,C,'bo-');
ylabel('C');
title('t = 0');
subplot(2,2,2);
plot(X,dCdt,'bo-');
ylabel('dC / dt');
title('t = 0');
for it = 1:1:nT
    Cava = (MT-C)/MT;

    ix = Ngrid;
    dCdt(ix) = Vm*Cava(ix)*C(ix-1)/Dx(ix-1);

    for ix = Ngrid-1:-1:2;
        dCdt(ix) = Vm * (Cava(ix)*C(ix-1)/(Dx(ix-1)+Dx(ix)) -Cava(ix+1)*C(ix)/(Dx(ix-1)+Dx(ix)));
            
    end

    ix = 1;
    dCdt(ix) = -Vm*Cava(ix+1)*C(ix)/Dx(ix);
    
    C = C + dCdt*dt;
    end

Ctest = sum((C(1:Ngrid-1)+C(2:Ngrid)).*Dx/2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ctest=0;
% for i=1;1;Ngrid
%     
subplot(2,2,3);
plot(X,C,'bo-');
xlabel('Time, sec');
ylabel('C');
title(['t = ', num2str(T)]);
subplot(2,2,4);
plot(X,dCdt,'bo-');
xlabel('Time, sec');
ylabel('dC / dt');
title(['t = ', num2str(T)]);