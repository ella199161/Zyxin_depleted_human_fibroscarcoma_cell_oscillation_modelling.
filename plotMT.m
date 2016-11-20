l=length(x);
y=zeros(l);
y=KMT*3.1415926*5*5;
plot(x(end),MTu(end),'x',x(end),MTd(end),'>',x(end),MT0(end),'<',x(end),y(end),'-','LineWidth',2,'MarkerSize',10);
hold on
plot(x,MTu,x,MTd,x,MT0,x,y,'LineWidth',2);% 2 pixels
h=legend('MTu','MTd','MT0','KMT');
xlabel('time (t)');
ylabel('Content');
title_handle = title(parameter);%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
hold off;

