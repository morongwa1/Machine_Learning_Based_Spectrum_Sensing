clear all;
clc;
%close all;

n =300;                     
ps=ones(100,15);             
yz=ones(100,15);            
db=ones(100,15);             
qu=ones(100,15);             
pfa1=0.10;              
snr1=-14:1:0;               
k= length(snr1);

 for i = 1:k                 
    s1 = 0;                   
    snr(i) = power(10,snr1(i)/10);
        
    for kk = 1:99          
       t = 1/n:1/n:4;              
       db(kk,i)=snr(i);
       
       x = sin(5*pi*t);                
        x = ( x > 0);                
        x  =2*x-1;                   
         
        noise = randn(1,length(t));  
        pn = (std(noise)).^2;        
        a = sqrt(2*snr(i));       
        xx = a.*cos(2*pi*4*t);        
        y  = x.*xx;                  
               
        yy = y+noise;                
        
        r1 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa1));    
        yz(kk,i)=r1;
        
        sum1 = sum(abs(yy).^2);      
        sum2=sum1/3.80;               
        ps(kk,i)=sum2;
        if sum2>r1
            s1 = s1+1;                 
            qu(kk,i)=1;
        else
            qu(kk,i)=-1;
        end
       
    end
    
    
p1(i) = s1/100;             

pd(i) = 0.40*erfc((erfcinv(2*pfa1)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i))); 

 end

ps_train=ps(:);       
yz_train=yz(:);
db_train=db(:);
qu_train=qu(:);
A=[ps_train,db_train,yz_train,qu_train];  
fid=fopen('CNNtrain.txt','wt');     
[h,l]=size(A);                                        
for ii=1:h
    for jj=1:l
        if jj==l
            fprintf(fid,'%g\n',A(ii,jj));      
        else
            fprintf(fid,'%g\t',A(ii,jj));      
        end
    end
end
fclose(fid);

            
figure(3);
plot(snr1,pd,'o-b',snr1,p1,'-*r');
grid on
legend('pf=0.10 Energy detection ','pf=0.10 Actual detection',4);
title ('Comparison of energy detection and actual detection probability')
xlabel('SNR/dB');
ylabel('pd Detection probability');