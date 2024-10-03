%beta=5 TDS VS 明文幅值(bata=5，4G，不考虑噪声)




clc
clear

dataTT=[400,80,40,26,20,16,13,11,10,9,8];

acf=[];
ACF_up=[];
ACF_down=[];
ACF_peak=[];

dmi=[];
DMI_up=[];
DMI_down=[];
DMI_peak=[];
for r=1:length(dataTT)

dataT=dataTT(r);


sim('beta5al1_5');
plot(x);
hold on;
plot(B1);

xx=x(end-30e4:end);
yy=x(end-30e4:end);

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


% 
% %_____________________________互信息（mi.m文件）___________________________
% for i1=T1:T2
%     k5((i1-T1)+1)=mi(xx(n:m),yy(n+i1:m+i1));
% end





M1(r)=mean(k4);
D1=sqrt((sum((k4-M1(r)).^2))/length(k4));
AACF1(r)=M1(r)+D1;AACF2(r)=M1(r)-D1;
peak(r)=max(abs(k4(1200-100:1200+100)));
acf=[acf;k4];
ACF_up=[ACF_up;AACF1(r)];
ACF_down=[ACF_down;AACF2(r)];
ACF_peak=[ACF_peak;peak(r)];





end
V=40./dataTT;
plot(V,ACF_down);
hold on;
plot(V,ACF_up);
hold on;
plot(V,ACF_peak,'.');






