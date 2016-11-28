Ca1=transpose(Cak);

[~,h]=contourf(Ca1);
set(h,'linestyle','none');
colorbar;
xlabel('time (t)');
ylabel('L(grid)');

title_handle = title('Ca kymo');
