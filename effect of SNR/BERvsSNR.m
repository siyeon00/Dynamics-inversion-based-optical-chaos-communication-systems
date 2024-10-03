%% 
clc
clear
dataT=40;%data transmission speed is set 1Gbps here

sim('system_reverse');


%% training
trainLen=5e4;
initLen=100;
tau=B1(1e5*dataT:5e6);%plaintext
xin=x(1e5*dataT:5e6);%chaos signal
xt=x(1e5*dataT+1200:5e6+1200); % 30ns delayed chaos signal
data=[xin(1:trainLen,:),xt(1:trainLen,:)];%input to RC
Yt=tau(initLen+1:trainLen)';%desired output
inSize = 2;    %Input dimension
outSize = 1;    %Output dimension


resSize = 1000;   % the internal units of RC
a = 0.7; % leaking rate  
rand( 'seed', 42 );
Win = (rand(resSize,1+inSize)-0.5) .* 1;    
W = rand(resSize,resSize)-0.5;    

disp 'Computing spectral radius...';
opt.disp = 0;
rhoW = abs(eigs(W,1,'LM',opt));  %Spectral radius  
disp 'done.' 
W = W .* (0.9/rhoW);   

X = zeros(1+inSize+resSize,trainLen-initLen);



xx = zeros(resSize,1);
for t = 1:trainLen
   u = data(t,:)'; %ÊäÈë
    xx=(1-a)*xx+ a*tanh( Win*[1;u] + W*xx ); %RC state
    if t > initLen
		X(:,t-initLen) = [1;u;xx];   
	end
end

reg = 1e-8;
Wout = (     (X*X' + reg*eye(1+inSize+resSize)) \ (X*Yt')  )'; 

ytrain=Wout*X; %RC output




figure(1);
plot(Yt(1:end));
hold on;
plot(ytrain(1:end));
title('Training');

%% testing

ber=[];
snr=2:30;
for jj=1:length(snr)
testLen=5000*40;
rec=x(1e5*dataT+19e3*dataT:1e5*dataT+19e3*dataT+testLen+1200-1);
 noi=awgn(rec,snr(jj),'measured');%signal with noise

xinr=noi(1:testLen);
xtr=noi(1+1200:testLen+1200);
data2=[xinr,xtr]; %received signal



Y2 = zeros(outSize,testLen);

for t = 1:testLen 
   
     u = data2(t,:)'; 
	xx = (1-a)*xx + a*tanh( Win*[1;u] + W*xx );
	y = Wout*[1;u;xx];
	Y2(:,t) = y;
  
end
figure(2);
plot(Y2);
hold on;
plot(tau(19e3*dataT:19e3*dataT+testLen));
title('testing');

errorLen = testLen;
corr(tau(19e3*dataT:19e3*dataT+testLen-1),Y2')


e=Y2;
%low-pass filter
Fp=0.5*10^6;  
Fs=1.5*10^6;  
Ft=40*10^6;%sampling rate/kHz
Rp=1;  %dB
Rs=30;   %dB

wp=2*pi*Fp;
ws=2*pi*Fs;

Y=fft(e,Ft);                   

[nt,wn]=buttord(wp,ws,Rp,Rs,'s');    
                             

[bl,al]=butter(nt,wn,'s');                
[num,den]=bilinear(bl,al,Ft);          

e1=filter(num,den,e);               %The filtered signal

shift=23;

e2=zeros(length(e1),1);
len=length(e1)-shift;
e2(1:end-shift)=e1(end-len+1:end);
figure(3);
plot(e);
hold on
plot(e2);
title('Effect of filtering');

%% decoding

bitlen=40;
Re=e2(1:end-shift); 
NN=fix(length(Re)/bitlen); 
Rm=[];
Re1=mean(Re);

for ii=1:NN
     if  Re(round((bitlen)/2)+(ii-1)*bitlen+1)>=Re1
        Rm=[Rm 1.5];
    else
        Rm=[Rm 0];
    end
end

Bit_real=tau(19e3*dataT+2:19e3*dataT+testLen);

Bit=[];
for j=1:NN
    if Bit_real(1+(j-1)*bitlen+round(bitlen/2))== 1.5
        Bit=[Bit;1.5];
    else
        Bit=[Bit; 0];
    end
end


BER=(NN-diff)/NN;
ber=[ber;BER]
end

