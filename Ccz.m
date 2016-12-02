figure;
xn=load('xdir.dat');
for i=0:1:100;
    
    fz=load(['Z' num2str(i) '.dat']);
    
    plot(xn,fz);
%    ylim([0, 750]);
    %colormap(gray);
    clear fz;       
     
    %make gif
    f = getframe(gcf);
    Im = frame2im(f);
    [Ima,map]=rgb2ind(Im,256);
    if i==0
    imwrite(Ima,map,'metroplis.gif','gif','DelayTime',5,'LoopCount',inf);
    else
    imwrite(Ima,map,'metroplis.gif','gif','WriteMode','append','DelayTime',5);
    end
    

end
