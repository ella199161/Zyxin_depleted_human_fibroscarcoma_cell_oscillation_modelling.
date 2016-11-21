%x=0:0.2:30.1;
l=length(x);
y=zeros(l);
y=Kca*3.1415926*5*5;
plot(x(end),Cau(end),'o',x(end),Cad(end),'>',x(end),y(end),'-','LineWidth',2);
hold on;
plot(x,Cau,x,Cad,x,y,'LineWidth',2);

h=legend('Cau','Cad','Kca');%,'Cmu','Cmd','MTu','MTd','MT0');
hold off;
title_handle = title('This is the original title');
set(title_handle,'String','Cau vs Cad');