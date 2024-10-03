%beta=5 TDS VS 明文幅值(bata=5，4G，不考虑噪声)




clc
clear

dataTT=[400,80,40,26,20,16,13,11,10,9,8];

acf=[];



for r=1:length(dataTT)

dataT=dataTT(r);


sim('beta5al1_5_v');
plot(x);
hold on;
plot(B1);

xx=x(end-30e4:end);
yy=B1(end-30e4:end);

T1=0;T2=3000;
%T1=0;T2=2000;


%_______________________________数据长_____________________________________
n=1;m=length(xx)-T2;

%_______________________________自相关_____________________________________
for i=T1:T2
    %k2((i-T1)+1)=mean((xx(n:m)-mean(xx(n:m))).*(xx(n+i:m+i)-mean(xx(n:m))));
    %k3((i-T1)+1)=(mean((xx(n:m)-mean(xx(n:m))).^2)*mean((xx(n+i:m+i)-mean(xx(n:m))).^2))^0.5;
    k2((i-T1)+1)=mean((xx(n:m)-mean(xx(n:m))).*(yy(n+i:m+i)-mean(xx(n:m))));
    k3((i-T1)+1)=(mean((xx(n:m)-mean(xx(n:m))).^2)*mean((yy(n+i:m+i)-mean(xx(n:m))).^2))^0.5;
    k4((i-T1)+1)=k2((i-T1)+1)/k3((i-T1)+1);
end




acf=[acf;k4(1)]


end
figure
V=40./dataTT;
plot(V,acf);








