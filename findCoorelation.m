
xu=length(Cmu);
corrud=zeros(xu);
temp=corrcoef(Cmu,Cmd);
corrud(1)=temp(2,1);
for i=1:1:(xu-100)
    temp=corrcoef(Cmu(i:end),Cmd(1:(end-i+1)));
    corrud(i+1)=temp(2,1);
end
lco=length(corrud);
plot(x(1:lco),corrud);
h=legend('corrolation between Cmu and Cmd');

title_handle = title('This is the original title');
set(title_handle,'String','Cross corrlation');