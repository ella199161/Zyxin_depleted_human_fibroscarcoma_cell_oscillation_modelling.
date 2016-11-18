clear

L = 50;
Ctotal = 10;
MT = 10;
Vm = 0.5;

dt = 0.001;
T =100.0;
nT = T/dt;

Dx(1,1) = 0.03;
for i=1:1:100
    Dx(i+1,1) = Dx(i,1) + Dx(i,1)*Dx(i,1)/2;
    if(Dx(i+1,1)>=0.5)
        break;
    end
end
ndg = length(Dx);
nad = (50-2*sum(Dx))/Dx(ndg) - mod((50-2*sum(Dx))/Dx(ndg),1);
nad = nad/2 - mod(nad/2,1);
Dx(ndg:ndg+nad) = Dx(ndg);
for i=ndg+nad+1:1:2*(ndg+nad)
    Dx(i) = Dx(2*(ndg+nad)-i+1);
end
N = length(Dx);

Ngrid = N+1;
X = zeros(Ngrid,1);
for in = 1:1:N
    X(in+1) = X(in)+Dx(in);
end

C = zeros(Ngrid,1);
ns = ndg/5 - mod(ndg/5,1);
C(5:ns+5) = Ctotal / (sum(Dx(5:ns+4))+Dx(ns+5)/2+Dx(4)/2);
% Ctest = sum((C(1:Ngrid-1)+C(2:Ngrid)).*Dx/2)


dCdt = zeros(Ngrid,1);

subplot(2,2,1);
plot(X,C,'bo-');
ylabel('C');
title('t = 0');
subplot(2,2,2);
plot(X,dCdt,'bo-');
ylabel('dC / dt');
title('t = 0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%option 1

for it = 1:1:nT
    Cava = (MT-C)/MT;

    ix = Ngrid;
    dCdt(ix) = Vm*Cava(ix)*C(ix-1)/Dx(ix-1);

    for ix = Ngrid-1:-1:2;
        dCdt(ix) = Vm * (Cava(ix)*C(ix-1)/(Dx(ix-1)+Dx(ix)) - ...
            Cava(ix+1)*C(ix)/(Dx(ix-1)+Dx(ix)));
    end

    ix = 1;
    dCdt(ix) = -Vm*Cava(ix+1)*C(ix)/Dx(ix);
    
    C = C + dCdt*dt;
end

Ctest = sum((C(1:Ngrid-1)+C(2:Ngrid)).*Dx/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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