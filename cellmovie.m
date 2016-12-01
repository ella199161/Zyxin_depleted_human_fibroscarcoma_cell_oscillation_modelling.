loadall
figure(1);
position=0;
v=200;
dx=x(2)-x(1);
maxv=85;
for i=1:5:length(x)
    i
    Cm=[Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:)];
   position=position+(Cau(i)-Cad(i))/maxv*5*200*dx;
    [~,h]=contourf(X-65+position, (1:2:20),Cm,[0:0.01:max(max(Cm))]);
    set(h,'linestyle','none');   
    caxis([0 1.1]); 
    colorbar;
    xlabel('r');
    ylabel('L');
    range=500;
width=20;
ylim([-30,width+30])
xlim([-100,range+50])
    title_handle = title('Cell oscillation');
        
    text(-45,-25,['time=' num2str(i*dx) 'h'],'color','r');
    %make gif
    f = getframe(gcf);
    Im = frame2im(f);
    [Ima,map]=rgb2ind(Im,256);
    if i==1
        imwrite(Ima,map,'cellmove3.gif','gif','DelayTime',0.02,'LoopCount',inf);
    else
        imwrite(Ima,map,'cellmove3.gif','gif','WriteMode','append','DelayTime',0.02);
    end
    
end
   