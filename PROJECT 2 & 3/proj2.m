
clc;
clear all; 
close all;

%Erwtima A.1-------------------------------------------------------

%mesw ekfwnisis
T=0.001;
over=5;
Ts=T/over;
A=3;
a=0.5;

%dimiourgia palmwn srrc
[phi1, t1] = srrc_pulse(T,Ts,A,a);
Nf=2048;
fft_phi1=fftshift(fft(phi1,Nf)*Ts);
%apoluti timi fft
abs_fft=abs(fft_phi1);
Fs=1/Ts;
F=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
spec_phi1 = abs_fft.^2;

%dimiourgia semilogy
figure();
semilogy(F,spec_phi1)
title('Energy Spectrum of SRRC Pulses');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum');

%Erwtima A.2-------------------------------------------------------

%mesw ekfwnisis
N=50;

%dimiourgia akolouthias N bits
b=(sign(randn(N,1))+1)/2;
X = bits_to_2PAM(b);

%Xd sima
Xd=zeros(1,N*over);
Xdt=linspace(0,N*T,N*over);
i=1;
for t=1:over:length(Xd)
Xd(t)=X(i);
i=i+1;
end

%Suneliksi
tmin=Xd(1)+t1(1);
tmax=Xd(end)+t1(end);
Xt=conv(phi1,Xd);
dt=linspace(tmin,tmax,length(Xt));

%dimiourgia plot
figure();
plot(dt,Xt)
title('Convolution');
ylabel('X(t)');
xlabel('Time (t)');

%Ypologismos fasmatikis puknotitas isxuos 
%median
sum=0;
for i=1:length(Xt)
sum=sum+Xt(i);
end
median=sum/length(Xt);

%variance
var=0;
for i=1:length(Xt)
var=var+(Xt(i)-median)^2;
end
variance=var/length(Xt);

%Fasmatiki puknotita isxuos
Sx=(variance/T)*spec_phi1;

%Erwtima A.3-------------------------------------------------------

fft_Xt=fftshift(fft(Xt,Nf)*Ts);
spec_Xt = abs(fft_Xt).^2;
Ttotal=length(dt)*Ts;
Px=spec_Xt/Ttotal;
F_Px=linspace(-Fs/2,Fs/2,length(Px));

%dimiourgia periodiagramatos plot
figure();
plot(F_Px,Px)
title('Periodogram of SRRC pulses')
xlabel('Frequency (Hz)');
ylabel('Periodogram');

%dimiourgia semilogy
figure();
semilogy(F_Px,Px)
title('Periodogram of SRRC pulses')
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses');

K=100;

figure();
for p=1:K
    b=(sign(randn(N,1))+1)/2;
    X = bits_to_2PAM(b);
    %create Xd signal
    Xd=zeros(1,N*over);
    Xdt=linspace(0,N*T,N*over);
    i=1;

    for t=1:over:length(Xd)
        Xd(t)=X(i);
        i=i+1;
    end

%convolution
    tmin=Xd(1)+t1(1);
    tmax=Xd(end)+t1(end);
    Xt=conv(phi1,Xd);
    dt=linspace(tmin,tmax,length(Xt));
    fft_Xt=fftshift(fft(Xt,Nf)*Ts);
    spec_Xt = abs(fft_Xt).^2;
    Ttotal=length(dt)*Ts;
    Px(p,:)=spec_Xt/Ttotal;
    F_Px=linspace(-Fs/2,Fs/2,length(Px));

    semilogy(F_Px,Px(p,:));
    xlabel('Frequency (Hz)');
    ylabel('100 Periodograms of SRRC pulses');
    hold on;
end

hold off;

%semilogy
figure;
semilogy(F,Sx, 'r');
hold on;
semilogy(F, mean(Px), 'b');
title('Theoretical and experimental approach PSD : 2-PAM');
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses');
legend('Theoritical','Experimental');
hold off;

K2=1000;
N2=100;

figure();
for p=1:K2
    b=(sign(randn(N2,1))+1)/2;
    X = bits_to_2PAM(b);
    %create Xd signal
    Xd=zeros(1,N2*over);
    Xdt=linspace(0,N2*T,N2*over);
    i=1;

    for t=1:over:length(Xd)
        Xd(t)=X(i);
        i=i+1;
    end

%convolution
    tmin=Xd(1)+t1(1);
    tmax=Xd(end)+t1(end);
    Xt=conv(phi1,Xd);
    dt=linspace(tmin,tmax,length(Xt));
    fft_Xt=fftshift(fft(Xt,Nf)*Ts);
    spec_Xt = abs(fft_Xt).^2;
    Ttotal=length(dt)*Ts;
    Px(p,:)=spec_Xt/Ttotal;
    F_Px=linspace(-Fs/2,Fs/2,length(Px));

    semilogy(F_Px,Px(p,:));
    xlabel('Frequency (Hz)');
    ylabel('1000 Periodograms of SRRC pulses');
    hold on;
end

hold off;
%semilogarithmic plot
figure();
semilogy(F, Sx, 'r');
hold on;

semilogy(F, mean(Px), 'b');
title('Theoretical and experimental approach PSD - increasing K,N');
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses');
legend('Theoritical','Experimental');
hold off;

%Erwtima A.4-------------------------------------------------------

%mesw ekfwnisis
N3=N/2;

X2 = bits_to_4PAM(b);

%dimiourgia Xd simatos
Xd2=zeros(1,N3*over);
Xdt2=linspace(0,N3*T,N3*over);
i=1;

for t=1:over:length(Xd2)
Xd2(t)=X2(i);
i=i+1;
end

%convolution
tmin2=Xd2(1)+t1(1);
tmax2=Xd2(end)+t1(end);
Xt2=conv(phi1,Xd2);
dt2=linspace(tmin2,tmax2,length(Xt2));

%dimiourgia plot suneliksis
figure();
plot(dt2,Xt2)
title('Convolution for 4-PAM');
ylabel('X(t)');
xlabel('Time t');

%upologismos fasmatikis puknotitas isxuos
%median
sum2=0;
for i=1:length(Xt2)
    sum2=sum2+Xt2(i);
end
median2=sum2/length(Xt2);

%variance
var2=0;
for i=1:length(Xt2)
    var2=var2+(Xt2(i)-median2)^2;
end
variance2=var2/length(Xt2);
Sx2=(variance2/T)*spec_phi1;

for p=1:K
    b=(sign(randn(N,1))+1)/2;
    X2 = bits_to_4PAM(b);
%dimiourgia Xd simatos
    Xd2=zeros(1,N3*over);
    Xdt2=linspace(0,N3*T,N3*over);
    i=1;
    for t=1:over:length(Xd2)
        Xd2(t)=X2(i);
        i=i+1;
    end
%convolution
    tmin2=Xd2(1)+t1(1);
    tmax2=Xd2(end)+t1(end);
    Xt2=conv(phi1,Xd2);
    dt2=linspace(tmin2,tmax2,length(Xt2));
    fft_Xt2=fftshift(fft(Xt2,Nf)*Ts);
    spec_Xt2 = abs(fft_Xt2).^2;
    Ttotal2=length(dt2)*Ts;
    Px2(p,:)=spec_Xt2/Ttotal2;
    F_Px2=linspace(-Fs/2,Fs/2,length(Px2));
end

%semilogy
figure();
semilogy(F,Sx2, 'r');
hold on;
semilogy(F, mean(Px2), 'b');
title('Theoretical and experimental approach PSD : 4-PAM');
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses');
legend('Theoritical','Experimental');

hold off;

figure();
plot(F, Sx, 'r');
hold on;
plot(F, Sx2, 'b');
xlabel('Frequency (Hz)');
ylabel('Bandwidth for 2-PAM and 4-PAM');
legend('2-PAM','4-PAM');

hold off;

%semilogy plot
figure();
semilogy(F, Sx, 'r');
hold on;
semilogy(F, Sx2, 'b');
xlabel('Frequency (Hz)');
ylabel('Bandwidth for 2-PAM and 4-PAM');
legend('2-PAM','4-PAM');

hold off;

%Erwtima A.5-------------------------------------------------------

%mesw ekfwnisis
T2=2*T;
over2=2*over;
Ts2=T2/over2;

%dimiourgia SRRC palmou
[phi2, t2] = srrc_pulse(T2,Ts2,A,a);
Xd3=zeros(1,N*over2);
Xdt3=linspace(0,N*T2,N*over2);

i=1;
for t=1:over2:length(Xd3)
    Xd3(t)=X(i);
    i=i+1;
end

%convolution
tmin3=Xd3(1)+t2(1);
tmax3=Xd3(end)+t2(end);
Xt3=conv(phi2,Xd3);
dt3=linspace(tmin3,tmax3,length(Xt3));
fft_Xt3=fftshift(fft(Xt3,Nf)*Ts);
spec_Xt3 = abs(fft_Xt3).^2;
Ttotal3=length(dt3)*Ts;
Px3(p,:)=spec_Xt3/Ttotal3;
F_Px3=linspace(-Fs/2,Fs/2,length(Px3));

%plot periodogram
figure();
plot(F_Px3,Px3)
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses - T2=2T');

%semilogy
figure();
semilogy(F_Px3,Px3)
xlabel('Frequency (Hz)');
ylabel('Periodogram of SRRC pulses - T2=2T');

K3=100;
N3=50;

figure();
for p=1:K3
    b=(sign(randn(N3,1))+1)/2;
    X3 = bits_to_2PAM(b);
%dimiourgia Xd simatos
    Xd3=zeros(1,N3*over2);
    Xdt3=linspace(0,N3*T2,N3*over2);
    i=1;
for t=1:over2:length(Xd3)
    Xd3(t)=X3(i);
    i=i+1;
end

%convolution
    tmin3=Xd3(1)+t2(1);
    tmax3=Xd(end)+t2(end);
    Xt3=conv(phi2,Xd3);
    dt3=linspace(tmin3,tmax3,length(Xt3));
    fft_Xt3=fftshift(fft(Xt3,Nf)*Ts2);
    spec_Xt3 = abs(fft_Xt3).^2;
    Ttotal3=length(dt3)*Ts2;
    Px3(p,:)=spec_Xt3/Ttotal;
    F_Px3=linspace(-Fs/2,Fs/2,length(Px3));
    semilogy(F_Px3,Px3(i,:));
    xlabel('Frequency (Hz)');
    ylabel('100 Periodograms of SRRC pulses - T2=2T');
    hold on;
end

hold off;
K4=1000;
N4=100;

figure();
for p=1:K4
    b=(sign(randn(N4,1))+1)/2;
    X = bits_to_2PAM(b);
%dimiourgia Xd simatos
    Xd=zeros(1,N4*over2);
    Xdt=linspace(0,N4*T2,N4*over2);
    i=1;
for t=1:over2:length(Xd)
    Xd(t)=X(i);
    i=i+1;
end

%convolution
    tmin=Xd(1)+t1(1);
    tmax=Xd(end)+t1(end);
    Xt=conv(phi1,Xd);
    dt=linspace(tmin,tmax,length(Xt));
    fft_Xt=fftshift(fft(Xt,Nf)*Ts2);
    spec_Xt = abs(fft_Xt).^2;
    Ttotal=length(dt)*Ts2;
    Px(p,:)=spec_Xt/Ttotal;
    F_Px=linspace(-Fs/2,Fs/2,length(Px));
    semilogy(F_Px,Px(i,:));
    xlabel('Frequency (Hz)');
    ylabel('1000 Periodograms of SRRC pulses - T2=2T');
    hold on;
end

hold off;