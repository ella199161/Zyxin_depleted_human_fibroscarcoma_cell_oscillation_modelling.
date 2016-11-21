figure(4)
position=0;


% dx=x(2)-x(1);
rad=0:0.4:20;
 dx1 = 0.01;
 L=130;
dx2 = 0.02;
 dx3 = 0.05;
 dx4 = 0.2;
 dx5 = 1;
x1 = 10;
 x2 = 5;
x3 = 4;
 x4 = 3;
 x5 = L - 2;
len1=0:dx1:(x1*dx1);
len2=(x1*dx1+dx2):dx2:(x1*dx1+x2*dx2);
len3=(x1*dx1+x2*dx2+dx3):dx3:(x1*dx1+x2*dx2+x3*dx3);
len4=(x1*dx1+x2*dx2+x3*dx3+dx4):dx4:(x1*dx1+x2*dx2+x3*dx3+x4*dx4);
len5=(x1*dx1+x2*dx2+x3*dx3+x4*dx4+dx5):dx5:(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5);
len6=(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4):dx4:(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4);
len7=(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3):dx3:(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3*x3);
len8= (x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3*x3+dx2):dx2:(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3*x3+dx2*x2); 
len9=(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3*x3+dx2*x2+dx1):dx1:(x1*dx1+x2*dx2+x3*dx3+x4*dx4+x5*dx5+dx4*x4+dx3*x3+dx2*x2+dx1*x1);
len=[len1 len2 len3 len4 len5 len6 len7 len8 len9];
Cc_concentration=zeros(173,51);


for i=1:1:300
    Ccdata=load(['Z' num2str(j) '.dat']);

 position=position+(Cau(i*6)-Cad(i*6))/100;
% for j=1:length(rad)
%     Cc_concentration(:,j*6)=Ccdata(:,j)./(len(j)+len(j+1))'/2;
% end
[~,h2]=contourf(len+position,rad,Ccdata');
set(h2,'linestyle','none');
 %caxis([0 0.03]); 
colorbar;
xlabel('r');
ylabel('L');
range=300;
width=20;

ylim([-30,width+30])
xlim([-100,range])
    title_handle = title('Cell oscillation');
        
    text(-80,-25,['time=' num2str(i*60) 'second'],'color','r');
    %make gif
    f = getframe(gcf);
    Im = frame2im(f);
    [Ima,map]=rgb2ind(Im,256);
    if i==1
        imwrite(Ima,map,'Cellmovie1.gif','gif','DelayTime',0.02,'LoopCount',inf);
    else
        imwrite(Ima,map,'Cellmovie1.gif','gif','WriteMode','append','DelayTime',0.02);
    end
    title_handle = title('Cell oscillation');
        
   
end
   