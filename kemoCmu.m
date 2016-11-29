Cmu1=transpose(Cmuk);
[~,h]=contourf(x,X,Cmu1);
set(h,'linestyle','none');
colorbar;
xlabel('time (t)');
ylabel('L');
title_handle = title('Cmu kymo');