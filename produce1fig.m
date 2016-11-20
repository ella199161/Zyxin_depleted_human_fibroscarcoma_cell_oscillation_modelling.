loadall
figure(1);
plot(x,Cc,x,Ca,x,Cmu,x,Cmd,'LineWidth',3);% 2 pixels
set(gca,'FontSize',15);
h=legend('Cc','Ca','Cmu','Cmd');
xlabel('time [t]','FontSize',30);
ylabel('Content','FontSize',30);
% saveas(1,'all_clean.eps','epsc')
% saveas(1,'all_clean.fig')

figure(2);
l=length(x);
y=Kca*3.1415926*5*5;
plot(x,Cau,x,Cad,'LineWidth',3);
set(gca,'FontSize',15);
h=legend('Cau','Cad');%,'Cmu','Cmd','MTu','MTd','MT0');
xlabel('time [h]','FontSize',30);
ylabel('Content','FontSize',30);
% saveas(2,'ca_clean.eps','epsc')
% saveas(2,'ca_clean.fig')

figure(3);
l=length(x);
plot(x,MTu,x,MTd,x,MT0,'LineWidth',3);% 2 pixels
set(gca,'FontSize',15);
h=legend('MTu','MTd','MT0');
xlabel('time [h]','FontSize',30);
ylabel('Content','FontSize',30);
% saveas(3,'mt_clean.eps','epsc')
% saveas(3,'mt_clean.fig')

figure(4);
Cmu1=transpose(Cmuk);
[~,h]=contourf(x(1:10:end),X,Cmu1(:,1:10:end));
set(h,'linestyle','none');
set(gca,'FontSize',15);
c=colorbar;
c.Label.String='Period[h]';
xlabel('time [h]','FontSize',30);
ylabel('L','FontSize',30);
% title_handle = title('Cmu kymo');
% saveas(4,'Cmu_clean.eps','epsc')
% saveas(4,'Cmu_clean.fig')
