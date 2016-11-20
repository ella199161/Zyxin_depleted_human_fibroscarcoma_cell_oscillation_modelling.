%(* ::Package:: *)

% Ccs=sum(Cc') ;
% Cas=sum(Ca') ;
% Cmus=sum(Cmu');
% Cmds=sum(Cmd');

plot(x(end),Cc(end),'o',x(end),Ca(end),'*',x(end),Cmu(end),'^',x(end),Cmd(end),'<','LineWidth',2,'MarkerSize',10);
hold on
plot(x,Cc,x,Ca,x,Cmu,x,Cmd,'LineWidth',2);% 2 pixels
h=legend('Cc','Ca','Cmu','Cmd');
xlabel('time (t)');
ylabel('Concerntration (N/um)');
title_handle = title(parameter);%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
hold off;


