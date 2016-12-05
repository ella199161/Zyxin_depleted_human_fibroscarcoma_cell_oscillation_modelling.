% loadall
% 
% 
% figure(1);
% plotsum
% figure(2);
% plotCa
% figure(3);
% kemoCmu
% figure(4);
% kemoCmd
% figure(5);
% kymoCa
% figure(6);
% plotMT
% figure(7);
% Ca013
% 
% %figure(8);
% %findCorrelation
loadall

figure(1);
plotsum
saveas(1,'all.eps','epsc')
saveas(1,'all.fig')

figure(2);
plotCa
saveas(2,'ca.eps','epsc')
saveas(2,'ca.fig')

figure(3);
kemoCmu
saveas(3,'Cmu.eps','epsc')
saveas(3,'Cmu.fig')

figure(4);
kemoCmd
saveas(4,'Cmd.eps','epsc')
saveas(4,'Cmd.fig')

figure(5);
kymoCa
saveas(5,'Cak.eps','epsc')
saveas(5,'Cak.fig')


figure(6);
plotMT
saveas(6,'mt.eps','epsc')
saveas(6,'mt.fig')

% figure(7);
% Ca013
% saveas(7,['.\L' num2str(L) '\Cai' num2str(L) '.eps'],'epsc')
% saveas(7,['.\L' num2str(L) '\Cai' num2str(L) '.fig'])

% figure(8);

% plot(w,psd);
% saveas(8,['.\L' num2str(L) '\corr' num2str(L) '.eps'],'epsc')
% saveas(8,['.\L' num2str(L) '\corr' num2str(L) '.fig'])

figure(7);
fdca
saveas(7,'fdca.eps','epsc')
saveas(7,'fdca.fig')

figure(8);
fdmt
saveas(8,'fdmt.eps','epsc')
saveas(8, 'fdmt.fig')


findtrans
peakpic
newfindtrans
% figure(11);
% kymoCc
% saveas(11,['Cck' num2str(L) '.eps'],'epsc')
[tau,af,w,psd]=compute_af_psd(Cau,dt);
if max(psd)>0
[maxtab, mintab] = peakdet(psd, 0.1*max(psd),w);
maxtab1=maxtab(1:end/2,:);
[a,b]=max(psd(1:end/2));
period=2*pi/w(b);
end