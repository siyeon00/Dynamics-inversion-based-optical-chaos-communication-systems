%���ʧ��


clc
clear

dataT=20;

TT1=100;
TT2=100;
move=-TT1:2:TT2;
ber=[];
for jj=1:length(move)
    


 sim('beta5al_dismatch');

trainLen=5e4;
initLen=100;

% tau=B1(1e5*dataT:5e6);
% xin=x(1e5*dataT:5e6);
% xt=x(1e5*dataT+1200:5e6+1200);



tau=B1(2e6:5e6);
xtrain=x(2e6:5e6);
xtrain_x=1:1:length(xtrain);
xtrain_gap=1:1/25:xtrain_x(end);
xtrain_up=interp1(xtrain_x,xtrain,xtrain_gap,'linear');%��ֵ���x
    
    
%����ps

xin_up=xtrain_up(1:25:25*(trainLen-1)+1)';
%xtr_up=rec_up(1+25*(1+1200-1):25:1+25*(testLen+1200-1))';
xt_up=xtrain_up(25*(1+1200-1)+1:25:1+25*(trainLen+1200-1))';
%xtr_up2=rec_up(25*(1+1200-1)+1:25:1+25*(testLen+1200-1));
data_up=[xin_up,xt_up];




xin=x(2e6:5e6);
xt=x(2e6+1200:5e6+1200);
data=[xin(1:trainLen,:),xt(1:trainLen,:)];

    
    
 %����ESN������
inSize = 2;    %����ά��k
outSize = 1;   
%�������м�����Խڵ����
resSize = 1000;   %������ά��ΪN
a = 0.7; % leaking rate  ���Կ����Ǵ����صĸ����ٶ�  k��ʱ������0.5
rand( 'seed', 42 );
Win = (rand(resSize,1+inSize)-0.5) .* 1;    %N*(k+1)
W = rand(resSize,resSize)-0.5;    %N*N

disp 'Computing spectral radius...';
opt.disp = 0;
rhoW = abs(eigs(W,1,'LM',opt));  %�װ뾶�ǣ��������ֵ�ľ���ֵ��  
%rhoW=abs(max(abs(eig(W))));
disp 'done.' 
W = W .* (0.9/rhoW);    %��һ���������װ뾶

%Ϊ״̬�ռ���������ڴ�ռ�
X = zeros(1+inSize+resSize,trainLen-initLen);
%�������յ�Ŀ�����
Yt=tau(initLen+1:trainLen)';



% data2=[xin(5e5:2e6),xt(5e5:2e6)];
%ѵ��
xx = zeros(resSize,1);
for t = 1:trainLen
   u = data_up(t,:)'; %����
    %u = data1(end-trainLen+1-T1+t-1); %����
    xx=(1-a)*xx+ a*tanh( Win*[1;u] + W*xx ); %�м�״̬�ĸ���
    if t > initLen
		X(:,t-initLen) = [1;u;xx];   %����initLen֮�󣬿�ʼ��¼������״̬
	end
end

reg = 1e-8;
Wout = (     (X*X' + reg*eye(1+inSize+resSize)) \ (X*Yt')  )'; 



ytrain=Wout*X;

% 
% figure(1);
% plot(Yt(1:end));
% hold on;
% plot(ytrain(1:end));
% title('ѵ��������-rc���');

%% ����
testLen=5000*40;
rec=x(2e6+5e5:2e6+5e5+testLen+1200-1+25);
rec_x=1:1:length(rec);
gap=1:1/25:rec_x(end);
rec_up= interp1(rec_x,rec,gap,'linear');   
    
    
xinr_up=rec_up(1:25:25*(testLen-1)+1)';
xtr_up=rec_up(move(jj)+25*(1+1200-1)+1:25:1+25*(testLen+1200-1)+move(jj))';
data2_up=[xinr_up,xtr_up];



Y2 = zeros(outSize,testLen);

for t = 1:testLen 
   % u = data(trainLen+t); %���Լ�����
     u = data2_up(t,:)'; %���Լ�����
	xx = (1-a)*xx + a*tanh( Win*[1;u] + W*xx );
	y = Wout*[1;u;xx];
	Y2(:,t) = y;
    
end


errorLen = testLen;



e=Y2;
%�˲�
%��ͨ�˲�
Fp=1.5*10^6;  %ͨ����ֹƵ��/kHz
Fs=2.5*10^6;  %�����ֹƵ��/kHz
Ft=40*10^6;%����Ƶ��/kHz
Rp=1;  %ͨ��˥��/dB
Rs=30;   %���˥��/dB

wp=2*pi*Fp;%ģ�����
ws=2*pi*Fs;

Y=fft(e,Ft);                    %���п��ٸ���Ҷ�任

[nt,wn]=buttord(wp,ws,Rp,Rs,'s');     %���ͨ�˲����Ľ����ͽ�ֹƵ�ʣ�ntΪ�˲�����С������wnΪ���ֹƵ�ʣ�
                             %wp��ws�ֱ�Ϊͨ��Ƶ�ʺͽ�ֹƵ�ʣ�RpΪͨ�����˥����Rs�����С˥����

[bl,al]=butter(nt,wn,'s');                %��S���Ƶ����Ӧ�Ĳ��� 
[num,den]=bilinear(bl,al,Ft);           %˫���Ա任ʵ��S��Z��ı任 

e1=filter(num,den,e);               %���˲�����ź�


%%
%�˲�����λ
shift=23;

e2=zeros(length(e1),1);
len=length(e1)-shift;
e2(1:end-shift)=e1(end-len+1:end);




%% ����
bitlen=20;
Re=e2(1:end-shift); 
NN=fix(length(Re)/bitlen); 
Rm=[];
Re1=mean(Re);

for ii=1:NN
     if  Re(round((bitlen)/2)+(ii-1)*bitlen)>=Re1
  
        Rm=[Rm 1.5];
    else
        Rm=[Rm 0];
    end
end

Bit_real=tau(5e5:5e5+testLen);

Bit=[];
for j=1:NN
    if Bit_real((j-1)*bitlen+bitlen/2)== 1.5
        Bit=[Bit; 1.5];
    else
        Bit=[Bit; 0];
    end
end


diff=0;
for j=1:NN
    if Rm(j)==Bit(j) 
        diff=diff+1;
    end
end

BER=(NN-diff)/NN;
ber=[ber;BER]

    
clearvars -except ber dataT move jj



    
end


plot(move(1:2:end),ber(1:2:end))