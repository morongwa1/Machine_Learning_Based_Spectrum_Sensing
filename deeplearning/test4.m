clear all;
clc;
%close all;

n =300;                     
%m=ones(15,100);
pfa1=0.025;              
pfa2=0.05;
pfa3=0.1;
snr1=-14:1:0;               
k  = length(snr1);
 for i = 1:k                 
    s1 = 0;                  
    s2 = 0; 
    s3 = 0; 
    snr(i) = power(10,snr1(i)/10);
    for kk = 1:99         
       t = 1/n:1/n:4;             
      
       x = sin(5*pi*t);                
        x = ( x > 0);                
        x  =2*x-1;                   
         
        noise = randn(1,length(t));  
        pn = (std(noise)).^2;        
        a = sqrt(2*snr(i));       
        xx = a.*cos(2*pi*4*t);         
        y  = x.*xx;                  
       ps = mean(abs(y).^2);        
        yy = y+noise;               
        
        r1 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa1));    
        r2 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa2))
        r3 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa3))
        
        sum1 = sum(abs(yy).^2);      
        if sum1>r1
            s1 = s1+1;                 
        end
        if sum1>r2
            s2=s2+1;
        end
        if sum1>r3
            s3=s3+1;
        end    
    end
p1(i) = s1/100;             
p2(i) = s2/100;             
p3(i) = s3/100;             


pd(i) = 0.4*erfc((erfcinv(2*pfa1)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i)));
pd2(i) = 0.4*erfc((erfcinv(2*pfa2)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i)));
pd3(i) = 0.4*erfc((erfcinv(2*pfa3)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i)));
 end

figure(2);

plot(snr1,pd,'o-b',snr1,pd2,'-+r',snr1,pd3,'-*m');
%grid on
legend('pf=0.025','pf=0.05','pf=0.1',4);
title ('Based on energy detection-comparison of detection probabilities under different false alarm probabilities')
xlabel('SNR/dB');
ylabel('Detection probability');