clc;
clear all;
close all;


n =300;                     %
ps_2=ones(100,15);             %,
yz_2=ones(100,15);               %,
db_2=ones(100,15);             %,
pfa1=0.10;              %=0.1
snr1=-14:1:0;               %
k= length(snr1);

 for i = 1:k                 %£¬
    s1 = 0;                  %£¬ 
    snr(i) = power(10,snr1(i)/10);%
        
    for kk = 1:99         %£¬£¬ 
       t = 1/n:1/n:4;              %=1200,4*300£¬
       db_2(kk,i)=snr(i);
       
       x = sin(5*pi*t);                %
        x = ( x > 0);                %£¬
        x  =2*x-1;                   %£¬1Óë-1£¬
         
        noise = randn(1,length(t));  %
        pn = (std(noise)).^2;        %£¬
        a = sqrt(2*snr(i));       %
        xx = a.*cos(2*pi*4*t);       %    
        y  = x.*xx;                  %£¬
               
        yy = y+noise;                %
        
        r1 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa1));    %,£¬n=2TW>
        yz_2(kk,i)=r1;
        
        sum1 = sum(abs(yy).^2);      %,
        sum2=sum1/3.80;               %
        ps_2(kk,i)=sum2;
        if sum2>r1
            s1 = s1+1;                 %£¬£¬
        end
       
    end
    
    
p1(i) = s1/100;             %£º/kk

pd(i) = 0.40*erfc((erfcinv(2*pfa1)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i))); %

 end
 
A=load('svmtrain.txt');     %
training0=A([1:1500],[1:2]);        %
training=training0';        %
[train,quality]=mapstd(training);        %

ps_3=A([1:1500],1);
yz_3=A([1:1500],3);

gind=find(A(:,4)==1);           %
bind=find(A(:,4)==-1);
group(gind)=1;                  %
group(bind)=-1;
group=group';  %

 ps_test=ps_2(:);
 %yz_test=yz_2(:);
 db_test=db_2(:);
 B0=[ps_test,db_test];    %
 B=B0';
 B=mapstd('apply',B,quality);

s=svmtrain(training',group,'Method','SMO','Kernel_Function','polynomial','showplot',true);       %£¬£¬
sv_index=s.SupportVectorIndices';   %
beta=s.Alpha';      %
b=s.Bias;           %
mean_and_std_trans=s.ScaleData;     %£¬
check=svmclassify(s,training','showplot',true);     %,£¬£¬£¬
err_rate=1-sum(group==check)/length(group)      %
solution=svmclassify(s,B');            %
solution=solution';
st=find(solution==1);                   %
sf=find(solution==-1);                   %


for i=1:k
    s2=0;
    for kk=1:99
        if yz_3(100*(i-1)+kk,1)<ps_3(100*(i-1)+kk,1)
        s2=s2+1;
        end
    end
    p3(i)=s2/100;        %
end

            
figure(5);
plot(snr1,pd,'o-b',snr1,p3,'-.+m');
grid on
legend('pf=0.10 Energy detection','pf=0.10SVM Testing');
title ('Comparison of energy detection and SVM detection probability')
xlabel('SNR/dB');
ylabel('pd Detection probability');