figure(2);

%Cm=zeros(10,91);
dx=x(2)-x(1);
for i=1:5:length(x)
    i
    Cm=[Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:); Cmuk(i,:)+Cmdk(i,:);Cmuk(i,:)+Cmdk(i,:)];
   
    [~,h]=contourf(X, (1:2:10),Cm(1:2:end,:),[0:0.05:max(max(Cm))]);
    set(h,'linestyle','none');   
    caxis([0 1.1]); 
    c=colorbar;
    ylabel(c,'polarity factors content')
    xlabel('r');
    ylabel('L');
    %title_handle = title('polarity factors');
        
    text(X(5),2,['time=' num2str(i*dx) 'h'],'color','r');
    %make gif
    f = getframe(gcf);
    Im = frame2im(f);
    [Ima,map]=rgb2ind(Im,256);
    if i==1
        imwrite(Ima,map,'Cm_small.gif','gif','DelayTime',0.1,'LoopCount',inf);
    else
        imwrite(Ima,map,'Cm_small.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    
end
   