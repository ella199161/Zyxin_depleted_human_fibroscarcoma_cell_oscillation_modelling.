 l=length(x);
 y=6*ones(1,l);
% y(1,:)=6;
% plot(x(end),MT1u(end),'x',x(end),MT1d(end),'>',x(end),y(end),'-','LineWidth',2,'MarkerSize',10);
% hold on
x=x/3600;
plot(x,MT1u,x,MT1d,x,y,'LineWidth',2);% 2 pixels
h=legend('MTu','MTd','KMT');
xlabel('time [h]');
ylabel('MT');
title_handle = title('MT and KMT');%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
% hold off;