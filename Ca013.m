l=length(x);
y=zeros(l);
y=Kca;
plot(x(end),Ca0u(end),'o',x(end),Ca0d(end),'x',x(end),Ca1u(end),'>',x(end),Ca1d(end),'<',x(end),Ca3u(end),'*',x(end),Ca3d(end),'+',x(end),y(end),'-',x(end),MT0u(end),'o',x(end),MT0d(end),'x',x(end),MT1u(end),'>',x(end),MT1d(end),'<',x(end),MT3u(end),'*',x(end),MT3d(end),'+');
hold on
plot(x,Ca0u,x,Ca0d,x,Ca1u,x,Ca1d,x,Ca3u,x,Ca3d,x,y,x,MT0u,x,MT0d,x,MT1u,x,MT1d,x,MT3u,x,MT3d,'LineWidth',2);

h=legend('Ca0u','Ca0d','Ca1u','Ca1d','Ca3u','Ca3d','Kca','MT0u','MT0d','MT1u','MT1d','MT3u','MT3d');
hold off
title_handle = title('This is the original title');
set(title_handle,'String','Cau vs Cad at 0 1 3');