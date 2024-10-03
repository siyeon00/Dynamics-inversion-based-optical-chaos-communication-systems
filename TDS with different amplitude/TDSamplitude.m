%beta=5 TDS VS 明文幅值(bata=5，4G，不考虑噪声)

clc
clear

AAL=0:0.15:3;
acf=[];
ACF_up=[];
ACF_down=[];
ACF_peak=[];

dmi=[];
DMI_up=[];
DMI_down=[];
DMI_peak=[];
for r=1:length(AAL)

AL=AAL(r);

dataT=20;

sim('beta5AL');
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




M1(r)=mean(k4);
D1=sqrt((sum((k4-M1(r)).^2))/length(k4));
AACF1(r)=M1(r)+D1;AACF2(r)=M1(r)-D1;
peak(r)=max(abs(k4(1200-100:1200+100)));
acf=[acf;k4];
ACF_up=[ACF_up;AACF1(r)];
ACF_down=[ACF_down;AACF2(r)];
ACF_peak=[ACF_peak;peak(r)];



end


plot(AAL,ACF_down);
hold on;
plot(AAL,ACF_up);
hold on;
plot(AAL,ACF_peak,'.');

plot(AAL,DMI_down);
hold on;
plot(AAL,DMI_up);
hold on;
plot(AAL,DMI_peak,'.');


