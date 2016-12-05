loadall
figure(1);
% plot(x(end),Cc(end),'o',x(end),Ca(end),'*',x(end),Cmu(end),'^',x(end),Cmd(end),'<','LineWidth',2,'MarkerSize',10);
% hold on
plot(x,Cc,x,Ca,x,Cmu,x,Cmd,'LineWidth',3);% 2 pixels
set(gca,'FontSize',15);
h=legend('Cc','Ca','Cmu','Cmd');
xlabel('Time [t]','FontSize',30);
ylabel('Content','FontSize',30);
%title_handle = title(parameter);%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
% hold off;
saveas(1,'all_clean.eps','epsc')
saveas(1,'all_clean.fig')

figure(2);
l=length(x);
y=zeros(l);
y=Kca*3.1415926*5*5;
% plot(x(end),Cau(end),'o',x(end),Cad(end),'>',x(end),y(end),'-','LineWidth',3);
% hold on;
plot(x,Cau,x,Cad,'LineWidth',3);
set(gca,'FontSize',15);
h=legend('Cau','Cad');%,'Cmu','Cmd','MTu','MTd','MT0');
xlabel('Time [h]','FontSize',30);
ylabel('Content','FontSize',30);
% hold off;
% title_handle = title('This is the original title');
% set(title_handle,'String','Cau vs Cad');
saveas(2,'ca_clean.eps','epsc')
saveas(2,'ca_clean.fig')

figure(4);
Cmu1=transpose(Cmuk);
[~,h]=contourf(x(1:10:end),X,Cmu1(:,1:10:end));
set(h,'linestyle','none');
set(gca,'FontSize',15);
c=colorbar;
ylabel(c,'Period[h]');
xlabel('Time [h]','FontSize',30);
ylabel('L[\mum]','FontSize',30);
% title_handle = title('Cmu kymo');
saveas(4,'Cmu_clean.eps','epsc')
saveas(4,'Cmu_clean.fig')

figure(5);
Cmd1=transpose(Cmdk);
[~,h]=contourf(x(1:10:end),X,Cmd1(:,1:10:end));
set(h,'linestyle','none');
set(gca,'FontSize',15);
c=colorbar;
ylabel(c,'Period[h]');
xlabel('Time [h]','FontSize',30);
ylabel('L[\mum]','FontSize',30);
% title_handle = title('Cmd kymo');
saveas(5,'Cmd_clean.eps','epsc')
saveas(5,'Cmd_clean.fig')

figure(3);
l=length(x);
% y=zeros(l);
% y=KMT*3.1415926*5*5;
% plot(x(end),MTu(end),'x',x(end),MTd(end),'>',x(end),MT0(end),'<',x(end),y(end),'-','LineWidth',2,'MarkerSize',10);
% hold on
plot(x,MTu,x,MTd,x,MT0,'LineWidth',3);% 2 pixels
set(gca,'FontSize',15);
h=legend('MTu','MTd','MT0');
xlabel('Time [h]','FontSize',30);
ylabel('Content','FontSize',30);
% title_handle = title(parameter);%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
% hold off;
saveas(3,'mt_clean.eps','epsc')
saveas(3,'mt_clean.fig')

figure(6);
Ca1k=Cak(:,1:end-1);
 Ccsum=Ca1k+Cmuk+Cmdk;
 [~,h]=contourf(x(1:10:end),X,Ccsum(1:10:end,:)');
 
 set(h,'linestyle','none');
 set(gca,'FontSize',15);
 colorbar;
 xlabel('Time [h]','FontSize',30);
ylabel('L[\mum]','FontSize',30);
saveas(6,'kymoclean.eps','epsc');
saveas(6,'kymoclean.fig');



% plot(x,MTu,x,MTd,x,MT0,'LineWidth',3);% 2 pixels
% h=legend('MTu','MTd','MT0');
% xlabel('time [h]','FontSize',30);
% ylabel('Content','FontSize',30);
% % title_handle = title(parameter);%set(title_handle,'String','Dc10,Vm0.5,km0.01,kf0.01,ka0.5,ka*5,kd0.5,kp1.5,kp*15,KMT4,b2,ko0.025,kc500,a5,Kca3');% need to implement rates for each graph
% % hold off;
% saveas(5,'mt_clean.eps','epsc')
% saveas(5,'mt_clean.fig')

ax=zeros(4,1);
for i = 1:4
    ax(i)=subplot(4,1,i);
end
for i = 1:4
    figure(i)
    h = get(gcf,'Children');
    newh = copyobj(h,5)
    for j = 1:length(newh)
posnewh = get(newh(j),'Position');
possub  = get(ax(i),'Position');
set(newh(j),'Position',...
[posnewh(1) possub(2) posnewh(3) possub(4)])
    end
    delete(ax(i));
end




