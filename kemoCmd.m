Cmd1=transpose(Cmdk);
[~,h]=contourf(x,X,Cmd1);
set(h,'linestyle','none');
colorbar;
xlabel('time (t)');
ylabel('L');
title_handle = title('Cmd kymo');
