clear all,close all;
clc;
tic
N=128;%
%N=256;
%n3=512;
%n4=1024;
Pf1=0.1;% 
Pf2=0.01;
Pf3=0.001;
Gamma_dB=linspace(-8,8,20); %
%SNR_dB=-8:0.4:8;
L=length(Gamma_dB);
 for i = 1:L      % 
        Gamma= power(10,Gamma_dB(i)/10);%
        Nd1=0; 
        Nd2=0;
        Nd3=0;
      for j=1:500
        t=(0:N-1)/10;  %
        Apt(i)=sqrt(Gamma);
        x=Apt(i).*cos(2*300*pi*t).*cos(2*pi*1000*t); %     
        noise=randn(1,N);%£¬£¬
        y=x+noise;%=+
        Y=fft(y,N); 
        %power_avg=mean(abs(x).^2);%
        signalPower = (sum(x).^2)/length(x);
        noisePower = sum(abs(noise).^2)/length(noise);%=mean(noise.^2)=std(noise).^2;
        %SNR_dB(i)=10*log10(sigPower/noisePower);  
        Lambda1= (N+sqrt(2*N)*qfuncinv(Pf1));%¡ª¡ª Pd_theory
        Lambda2=(N+sqrt(2*N)*qfuncinv(Pf2));%¡ª¡ª TO Pd_theory
        Lambda3=(N+sqrt(2*N)*qfuncinv(Pf3));
        T_simu(i)=sum(abs(Y).^2)/N;
        if T_simu(i) > Lambda1
            Nd1=Nd1+1;
        end
        if T_simu(i) > Lambda2
            Nd2=Nd2+1;
        end
         if T_simu(i) > Lambda3
            Nd3=Nd3+1;
        end
     end
    Pd1_simu(i)=Nd1/j;
    Pd2_simu(i)=Nd2/j;
    Pd3_simu(i)=Nd3/j;
Pd1_theory(i) = 0.5*erfc((Lambda1-N-N*Gamma)./sqrt(2*(2*N+4*N*Gamma)));%£¬SNR(i)=signalPower(i) 
Pd2_theory(i) = 0.5*erfc((Lambda2-N-N*Gamma)./sqrt(2*(2*N+4*N*Gamma)));
Pd3_theory(i) = 0.5*erfc((Lambda3-N-N*Gamma)./sqrt(2*(2*N+4*N*Gamma)));
 end

figure
plot(Gamma_dB,Pd1_simu,'dk',Gamma_dB,Pd1_theory,'--k','LineWidth',2);
hold on
plot(Gamma_dB,Pd2_simu,'+k',Gamma_dB,Pd2_theory,'--k','LineWidth',1.5,'MarkerSize',7);
hold on
plot(Gamma_dB,Pd3_simu,'ok',Gamma_dB,Pd3_theory,'--k','LineWidth',1.5);
xlabel('SNR');
ylabel('Pd');
title ('Curves of detection probability and signal-to-noise ratio under different false alarms')
legend('Pf=0.1Pdsimu','Pf=0.1Pdtheory','Pf=0.01Pdsimu','Pf=0.01Pdtheory','Pf=0.001Pdsimu','Pf=0.001Pdtheory')
grid on 
toc