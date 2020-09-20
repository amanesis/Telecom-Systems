clear all;
close all;
clc;

% A.1
%given data
N=100;
bit_seq = (sign(randn(4*N,1))+1)/2;

%-----------------------------------------------------------------
% A.2
A=1;
X = bits_to_4PAM(bit_seq, A);

%-----------------------------------------------------------------
% A.3
XI = X(1:N);
XQ = X(N+1:2*N);

%-----------------------------------------------------------------
% A.4
%given data
T=1;
over=10;
Ts=T/over;
a=1;
Nf=1024;
[phi, t] = srrc_pulse(T,Ts,A,a);

%create Xd signal for XI, XdI
XdI = zeros(1,N*over);
XdI_t = linspace(0,N*T,N*over);
i=1;

for p=1:over:length(XdI)
XdI(p) = XI(i);
i=i+1;
end

%create Xd signal for XQ, XdQ
XdQ = zeros(1,N*over);
XdQ_t = linspace(0,N*T,N*over);
i=1;

for p=1:over:length(XdQ)
XdQ(p) = XQ(i);
i=i+1;
end

%convolution for XI(t)
tmin=XdI_t(1)+t(1);
tmax=XdI_t(end)+t(end);
XI_t = conv(phi,XdI);
time_axis_XI=linspace(tmin,tmax,length(XI_t));

%convolution for XQ(t)
tmin2=XdQ_t(1)+t(1);
tmax2=XdQ_t(end)+t(end);
XQ_t = conv(phi,XdQ);
time_axis_XQ=linspace(tmin2,tmax2,length(XQ_t));

%plot convolution XI(t)
figure;
plot(time_axis_XI,XI_t)
ylabel('XI(t)');
xlabel('Time t');

%plot convolution XQ(t)
figure;
plot(time_axis_XQ,XQ_t)
ylabel('XQ(t)');
xlabel('Time t');
Fs=1/Ts;

%plot periodogram of XI(t)
fft_XIt=fftshift(fft(XI_t,Nf)*Ts);
spec_XIt = abs(fft_XIt).^2;
TtotalXIt=length(XI_t)*Ts;
PxXIt=spec_XIt/TtotalXIt;
F_PxXIt=linspace(-Fs/2,Fs/2,length(PxXIt));

figure;
plot(F_PxXIt,PxXIt)
xlabel('Frequency (Hz)');
ylabel('Periodogram of XI(t)');

%plot periodogram of XQ(t)
fft_XQt=fftshift(fft(XQ_t,Nf)*Ts);
spec_XQt = abs(fft_XQt).^2;
TtotalXQt=length(XQ_t)*Ts;
PxXQt=spec_XQt/TtotalXQt;
F_PxXQt=linspace(-Fs/2,Fs/2,length(PxXQt));

figure;
plot(F_PxXQt,PxXQt)
xlabel('Frequency (Hz)');
ylabel('Periodogram of XQ(t)');

%-----------------------------------------------------------------
% A.5
%given data
F0=2;

%create XImod signal for XI
XImod_time=time_axis_XI;
XImod = XI_t.*(2*cos(2*pi*F0*XImod_time));

%create XQmod signal for XQ
XQmod_time=time_axis_XQ;
XQmod = XQ_t.*((-2)*sin(2*pi*F0*XQmod_time));

%plot XImod(t)
figure;
plot(XImod_time,XImod)
xlabel('Time t');
ylabel('XImod(t)');

%plot XQmod(t)
figure;
plot(XQmod_time,XQmod)
xlabel('Time t');
ylabel('XQmod(t)');

%plot periodogram of XImod(t)
fft_XImod=fftshift(fft(XImod,Nf)*Ts);
spec_XImod = abs(fft_XImod).^2;
TtotalXImod=length(XImod)*Ts;
PxXImod=spec_XImod/TtotalXImod;
F_PxXImod=linspace(-Fs/2,Fs/2,length(PxXImod));

figure;
plot(F_PxXImod,PxXImod)
xlabel('Frequency (Hz)');
ylabel('Periodogram of XImod(t)');

%plot periodogram of XQmod(t)
fft_XQmod=fftshift(fft(XQmod,Nf)*Ts);
spec_XQmod = abs(fft_XQmod).^2;
TtotalXQmod=length(XQmod)*Ts;
PxXQmod=spec_XQmod/TtotalXQmod;
F_PxXQmod=linspace(-Fs/2,Fs/2,length(PxXQmod));

figure;
plot(F_PxXQmod,PxXQmod)
xlabel('Frequency (Hz)');
ylabel('Periodogram of XQmod(t)');

%-----------------------------------------------------------------
% A.6
Xmod = XImod + XQmod;

%plot Xmod(t)
figure;
plot(XImod_time,Xmod)
xlabel('Time t');
ylabel('Xmod(t)');

%plot periodogram of Xmod(t)
fft_Xmod=fftshift(fft(Xmod,Nf)*Ts);
spec_Xmod = abs(fft_Xmod).^2;
TtotalXmod=length(Xmod)*Ts;
PxXmod=spec_Xmod/TtotalXmod;
F_PxXmod=linspace(-Fs/2,Fs/2,length(PxXmod));

figure;
plot(F_PxXmod,PxXmod)
xlabel('Frequency (Hz)');
ylabel('Periodogram of Xmod(t)');

%-----------------------------------------------------------------
% A.8
%given data
SNR=20;

%variance
s2w=(10*(A^2))/(Ts*(10^(SNR/10)));

%white Gaussian noise
Wt=sqrt(s2w).*randn(1,length(Xmod));
Zt = Xmod + Wt;

%plot Z(t)
figure;
hold on;
plot(XImod_time,Xmod)
plot(XImod_time,Zt,'red')
hold off;
legend('before','after adding Gaussian noise')

%-----------------------------------------------------------------
% A.9
ZI=Zt.*cos(2*pi*F0*XImod_time);
ZQ=Zt.*(-sin(2*pi*F0*XQmod_time));

%plot ZI(t)
figure;
plot(XImod_time,ZI)
xlabel('Time t');
ylabel('ZI(t)');

%plot ZQ(t)
figure;
plot(XQmod_time,ZQ)
xlabel('Time t');
ylabel('ZQ(t)');

%plot periodogram of ZI(t)
fft_ZI=fftshift(fft(ZI,Nf)*Ts);
spec_ZI = abs(fft_ZI).^2;
TtotalZI=length(ZI)*Ts;
PxZI=spec_ZI/TtotalZI;
F_PxZI=linspace(-Fs/2,Fs/2,length(PxZI));

figure;
plot(F_PxZI,PxZI)
xlabel('Frequency (Hz)');
ylabel('Periodogram of ZI(t)');

%plot periodogram of ZQ(t)
fft_ZQ=fftshift(fft(ZQ,Nf)*Ts);
spec_ZQ = abs(fft_ZQ).^2;
TtotalZQ=length(ZQ)*Ts;
PxZQ=spec_ZQ/TtotalZQ;
F_PxZQ=linspace(-Fs/2,Fs/2,length(PxZQ));

figure;
plot(F_PxZQ,PxZQ)
xlabel('Frequency (Hz)');
ylabel('Periodogram of ZQ(t)');

%-----------------------------------------------------------------
% A.10
%convolution for ZI(t)
tmin3=XImod_time(1)+t(1);
tmax3=XImod_time(end)+t(end);
ZI_t = conv(phi,ZI)*Ts;
time_axis_ZI=linspace(tmin3,tmax3,length(ZI_t));

%convolution for ZQ(t)
tmin4=XQmod_time(1)+t(1);
tmax4=XQmod_time(end)+t(end);
ZQ_t = conv(phi,ZQ)*Ts;
time_axis_ZQ=linspace(tmin4,tmax4,length(ZQ_t));

%plot convolution ZI(t)
figure;
plot(time_axis_ZI,ZI_t)
ylabel('ZI(t)');
xlabel('Time t');

%plot convolution XQ(t)
figure;
plot(time_axis_ZQ,ZQ_t)
ylabel('ZQ(t)');
xlabel('Time t');

%plot periodogram of ZI(t)
fft_ZIt=fftshift(fft(ZI_t,Nf)*Ts);
spec_ZIt = abs(fft_ZIt).^2;
TtotalZIt=length(ZI_t)*Ts;
PxZIt=spec_ZIt/TtotalZIt;
F_PxZIt=linspace(-Fs/2,Fs/2,length(PxZIt));

figure;
plot(F_PxZIt,PxZIt)
xlabel('Frequency (Hz)');
ylabel('Periodogram of ZI(t)');

%plot periodogram of ZQ(t)
fft_ZQt=fftshift(fft(ZQ_t,Nf)*Ts);
spec_ZQt = abs(fft_ZQt).^2;
TtotalZQt=length(ZQ_t)*Ts;
PxZQt=spec_ZQt/TtotalZQt;
F_PxZQt=linspace(-Fs/2,Fs/2,length(PxZQt));

figure;
plot(F_PxZQt,PxZQt)
xlabel('Frequency (Hz)');
ylabel('Periodogram of ZQ(t)');

%-----------------------------------------------------------------
% A.11
Y = zeros(2,round((length(ZI_t)-4*A*Fs)/over));
i=1;

for p=2*A*Fs:Fs:(length(ZI_t)-1)-2*A*Fs
Y(1,i)=ZI_t(p);
Y(2,i)=ZQ_t(p);
i=i+1;
end

scatterplot(Y');

%-----------------------------------------------------------------
% A.12
est_XI = detect_4_PAM((Y(1,:)), A);
est_XQ = detect_4_PAM((Y(2,:)), A);

%-----------------------------------------------------------------
% A.13
%symbol errors
SER=0;

for p=1:length(XI)
    
if(est_XI(p)~=XI(p) | est_XQ(p)~=XQ(p))
SER = SER + 1;
end

end

%-----------------------------------------------------------------
% A.14
est_bit = PAM_4_to_bits([est_XI,est_XQ],A);

%-----------------------------------------------------------------
% A.15
%bit errors
BER=0;

for p=1:length(bit_seq)
    
if(est_bit(p)~=bit_seq(p))
BER = BER + 1;
end

end

%-----------------------------------------------------------------
% B.1
%given data
K=200;
SNRdb_vector=[0:2:16];
M=16;
bps=log2(M); % bits per symbol
symbols=N*K;
bits=4*N*K;
Pser=ones(1,length(SNRdb_vector));
Pber=ones(1,length(SNRdb_vector));
theorySer = ones(1,length(SNRdb_vector));
theoryBer = ones(1,length(SNRdb_vector));

%experimental
Kser=0;
Kber=0;
c=1;

for SNRdb=0:2:16
   
for p=1:K
 
[SER, BER] = QAM_16(N, SNRdb);
Kser=Kser+SER;
Kber=Kber+BER;

end

Pser(c)=Kser/symbols;
Pber(c)=Kber/bits;

%theoretical
theorySer(c)=3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))))-(1/
bps)*3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))));
theoryBer(c)=(1/bps)*3/2.*theorySer(c);
Kser=0;
Kber=0;
c=c+1;

end

%-----------------------------------------------------------------
% B.2
figure;
semilogy(SNRdb_vector,Pser);
hold on;
semilogy(SNRdb_vector,theorySer,'green');
legend('Experimental SER','Theoretical SER');
xlabel('SNR in dB');
ylabel('Symbol Error Rate for 16-QAM');

%-----------------------------------------------------------------
% B.3
figure;
semilogy(SNRdb_vector,Pber);
hold on;
semilogy(SNRdb_vector,theoryBer,'green');
legend('Experimental BER','Theoretical BER');
xlabel('SNR in dB');
ylabel('Bit Error Rate for 16-QAM');