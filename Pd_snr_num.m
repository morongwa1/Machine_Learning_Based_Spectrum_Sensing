clear all,close all;
clc;
tic
n1=128;%
n2=256;
n3=512;
%n4=1024;
Pf1=0.1;%
%Pf2=0.01;
Gamma_dB=linspace(-8,4,20); %
L=length(Gamma_dB);
for i = 1:L      % 
    Gamma= power(10,Gamma_dB(i)/10);%
    Nd1=0;
    Nd2=0;
    Nd3=0;
    for j=1:500
        t=(0:n1-1)/10;  %
        Apt(i)=sqrt(Gamma);
        x=Apt(i).*cos(2*300*pi*t).*cos(2*pi*1000*t); %
        noise=randn(1,n1);%，，
        y=x+noise;%=
        Y=fft(y,n1);
        %power_avg=mean(abs(x).^2);%
        signalPower = (sum(x).^2)/length(x);
        noisePower = sum(abs(noise).^2)/length(noise);%=
        %SNR_dB(i)=10*log10(sigPower/noisePower); 
        Lambda1= (n1+sqrt(2*n1)*qfuncinv(Pf1));%——门限_TO Pd_theory
        T_simu(i)=sum(abs(Y).^2)/n1;
        if T_simu(i) > Lambda1
            Nd1=Nd1+1;
        end
    end
    Pd1_simu(i)=Nd1/j;
    Pd1_theory(i) = 0.5*erfc((Lambda1-n1-n1*Gamma)./sqrt(2*(2*n1+4*n1*Gamma)));%，SNR(i)=signalPower(i)
    
    for m=1:500
        t=(0:n2-1)/10;  %
        Apt(i)=sqrt(Gamma);
        x=Apt(i).*cos(2*300*pi*t).*cos(2*pi*1000*t); %
        noise=randn(1,n2);%，
        y=x+noise;%
        Y=fft(y,n2);
        %power_avg=mean(abs(x).^2);%
        signalPower = (sum(x).^2)/length(x);
        noisePower = sum(abs(noise).^2)/length(noise);%=mean(noise.^2)=std(noise).^2;
        %SNR_dB(i)=10*log10(sigPower/noisePower); 
        Lambda2= (n2+sqrt(2*n2)*qfuncinv(Pf1));%
        T_simu(i)=sum(abs(Y).^2)/n2;
        if T_simu(i) > Lambda2
            Nd2=Nd2+1;
        end
    end
    Pd2_simu(i)=Nd2/m;
    Pd2_theory(i) = 0.5*erfc((Lambda2-n2-n2*Gamma)./sqrt(2*(2*n2+4*n2*Gamma)));%，SNR(i)=signalPower(i)
    
    for p=1:500
        t=(0:n3-1)/10;  %
        Apt(i)=sqrt(Gamma);
        x=Apt(i).*cos(2*300*pi*t).*cos(2*pi*1000*t); %
        noise=randn(1,n3);%，，
        y=x+noise;%
        Y=fft(y,n3);
        %power_avg=mean(abs(x).^2);%
        signalPower = (sum(x).^2)/length(x);
        noisePower = sum(abs(noise).^2)/length(noise);%=mean(noise.^2)=std(noise).^2;
        %SNR_dB(i)=10*log10(sigPower/noisePower); 
        Lambda3= (n3+sqrt(2*n3)*qfuncinv(Pf1));%——门限_TO Pd_theory
        T_simu(i)=sum(abs(Y).^2)/n3;
        if T_simu(i) > Lambda3
            Nd3=Nd3+1;
        end
    end
    Pd3_simu(i)=Nd3/p;
    Pd3_theory(i) = 0.5*erfc((Lambda3-n3-n3*Gamma)./sqrt(2*(2*n3+4*n3*Gamma)));%，SNR(i)=signalPower(i)
    
    
    
end

figure
plot(Gamma_dB,Pd1_simu,'dk',Gamma_dB,Pd1_theory,'--k','LineWidth',2);
hold on

plot(Gamma_dB,Pd2_simu,'+k',Gamma_dB,Pd2_theory,'--k','LineWidth',1.5,'MarkerSize',7);
plot(Gamma_dB,Pd3_simu,'ok',Gamma_dB,Pd3_theory,'--k','LineWidth',1.5)
xlabel('SNR');
ylabel('Pd');
title ('Curve of detection probability and signal-to-noise ratio for the same warning probability and different sampling points')
legend('N=128 Pdsimu','N=128 Pdtheory','N=256 Pdsimu','N=256 Pdtheory','N=512 Pdsimu','N=512 Pdtheory')
grid on
toc