clear all;
clc;
close all;

n = 5;
ps = 1;
SNR1 = -5;
SNR2 = -8;
SNR3 = -10;
% Sim_Times=10000; %Monter-Carlo times
% m=5;
T=0.001;
%  
W=5*10^4;
% 
Fs = 2*W;
m = T*W;
n = 2*T*W;
% F0=W;
% Fs=2;
% Sig=sqrt(2)*sin(2*pi*F0/Fs*t); %single tone samples, Fs=2F0
% 
snr1 = 10.^(SNR1/10);
snr2 = 10.^(SNR2/10);
snr3 = 10.^(SNR3/10);
pn = (1/snr1)*ps;
mu0 = n*pn;
sigma0 = sqrt(2*n)*pn;
mu = n*(pn+ps);
sigma = sqrt(2*n*(pn^2+2*pn*ps));
% [noi,x0,mu0,sigma0,m0] = cnoi( n,pn );
% sig = randn(n,1);
sig = 1;
% 
count = 5000;
% 
lambda = [200:20:600];
lambda1 = [500:20:900];
lambda2 = [700:30:1300];
% 
% tt = [-5:0.4:3];
% cc = 10.^tt;
% tt1 = [-1:0.1:1];
% cc1 = 10.^tt;
% cc2 = [-0.01:0.001:0.01];
for kk = 1:1:length(lambda);
ff = 0;
dd = 0;
ff1 = 0;
dd1 = 0;
ff2 = 0;
dd2 = 0;
for ii=1:1:count;
t = (kk-1)*n+1:kk*n;
init_phase = 1/6*pi;
sig=2*sin(2*pi*W/Fs*t+init_phase);   %£¬
%sig=1;
%dot(sig,sig)/n
noi1 = randn(1,n);
noi1 = sqrt(1/snr1)*noi1;
rec1 = noi1 + sig;  %+
noi2 = randn(1,n);
noi2 = sqrt(1/snr2)*noi2;
rec2 = noi2 + sig;
noi3 = randn(1,n);
noi3 = sqrt(1/snr3)*noi3;
rec3 = noi3 + sig;
sum0_1 = dot(noi1,noi1);
sum0_2 = dot(noi2,noi2);
sum0_3 = dot(noi3,noi3);
sum1_1 = dot(rec1,rec1);
sum1_2 = dot(rec2,rec2);
sum1_3 = dot(rec3,rec3);
%
if (sum0_1 > lambda(kk));
ff = ff+1;
end
if (sum1_1 > lambda(kk));
dd = dd+1;
end
if (sum0_2 > lambda1(kk));
ff1 = ff1+1;
end
if (sum1_2 > lambda1(kk));
dd1 = dd1+1;
end 
if (sum0_3 > lambda2(kk));
ff2 = ff2+1;
end
if (sum1_3 > lambda2(kk));
dd2 = dd2+1;
end
end
Pd_1(kk) = dd/count;
Pf_1(kk) = ff/count;
Pd_2(kk) = dd1/count;
Pf_2(kk) = ff1/count;
Pd_3(kk) = dd2/count;
Pf_3(kk) = ff2/count;
end
figure(1);
plot(Pf_1,Pd_1,'-o',Pf_2,Pd_2,'-s',Pf_3,Pd_3,'-*');
hold on;
%grid on;
xlabel('False alarm probability');
ylabel('Detection probability');
legend('SNR=-5dB','SNR=-8dB','SNR=-10dB');
title ('Energy detection comparison under different SNR');
%grid on;
hold on;