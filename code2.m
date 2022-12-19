
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

 for i = 1:k                 %，
    s1 = 0;                  %， 
    snr(i) = power(10,snr1(i)/10);%
        
    for kk = 1:99         %，， 
       t = 1/n:1/n:4;              %=1200,4*300，
       db_2(kk,i)=snr(i);
       
       x = sin(5*pi*t);                %
        x = ( x > 0);                %，
        x  =2*x-1;                   %，-1，
         
        noise = randn(1,length(t));  %
        pn = (std(noise)).^2;        %，
        a = sqrt(2*snr(i));       %
        xx = a.*cos(2*pi*4*t);       %    
        y  = x.*xx;                  %，
               
        yy = y+noise;                %
        
        r1 = (n+sqrt(n/2)*2*sqrt(2)*erfcinv(2*pfa1));    %,，n=2TW>
        yz_2(kk,i)=r1;
        
        sum1 = sum(abs(yy).^2);      %,
        sum2=sum1/3.80;               %
        ps_2(kk,i)=sum2;
        if sum2>r1
            s1 = s1+1;                 %，，
        end
       
    end
    
    
p1(i) = s1/100;             %：/kk

pd(i) = 0.40*erfc((erfcinv(2*pfa1)-snr(i)*sqrt(n/2))/sqrt(2+4*snr(i))); %

 end
 
A=load('svmtrain.txt');     %
training0=A([1:1500],[1:2]);        %
training=training0';        %
[train,quality]=mapstd(training);        %

ps_3=A([1:1500],1);
yz_3=A([1:1500],3);
length(ps_3);
length(yz_3);

gind=find(A(:,4)==1);           %-1
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
tic
s_linear = svmtrain(training',group,'Method','SMO','Kernel_Function','linear','showplot',true);nhhhhhhhhhhhhb 
toc
s_polynomial= svmtrain(training',group,'Method','SMO','Kernel_Function','polynomial','showplot',true);
toc
s_rbf = svmtrain(training',group,'Method','SMO','Kernel_Function','rbf','showplot',true); 
toc
s_quadratic = svmtrain(training',group,'Method','SMO','Kernel_Function','quadratic','showplot',true);
toc
s_mlp = svmtrain(training',group,'Method','SMO','Kernel_Function','mlp','showplot',true); 
toc


sv_index_linear = s_linear.SupportVectorIndices';   %
sv_index_polynomial = s_polynomial.SupportVectorIndices';   %
sv_index_rbf = s_rbf.SupportVectorIndices';   %
sv_index_quadratic = s_quadratic.SupportVectorIndices';   %
sv_index_mlp = s_mlp.SupportVectorIndices';   %


check_linear = svmclassify(s_linear,training','showplot',true);     %,，，，
check_polynomial = svmclassify(s_polynomial,training','showplot',true);  
check_rbf = svmclassify(s_rbf,training','showplot',true);  
check_quadratic = svmclassify(s_quadratic,training','showplot',true);  
check_mlp = svmclassify(s_mlp,training','showplot',true);  


Pd_linear = sum(group==check_linear)/length(group)*100    %
Pd_polynomial = sum(group==check_polynomial)/length(group)*100    %
Pd_rbf = sum(group==check_rbf)/length(group)*100     %
Pd_quadratic = sum(group==check_quadratic)/length(group)*100     %
Pd_mlp = sum(group==check_mlp)/length(group)*100     %

Pm_linear = (1-(sum(group==check_linear)/length(group)))*100    %
Pmpolynomial = (1-(sum(group==check_polynomial)/length(group)))*100      %
Pm_rbf = (1-(sum(group==check_rbf)/length(group)))*100      %
Pm_quadratic = (1-(sum(group==check_quadratic)/length(group)))*100      %
Pm_mlp = (1-(sum(group==check_mlp)/length(group)))*100      %

Pf_linear = sum(group~=check_linear)/length(group)*100    %
Pf_polynomial = sum(group~=check_polynomial)/length(group)*100      %
Pf_rbf = sum(group~=check_rbf)/length(group)*100     %
Pf_quadratic = sum(group~=check_quadratic)/length(group)*100      %
Pf_mlp = sum(group~=check_mlp)/length(group)*100     %

AC_linear = (sum(group==check_linear) - sum(group~=check_linear))/length(group)*100   %
AC_polynomial = (sum(group==check_polynomial) - sum(group~=check_polynomial))/length(group)*100    %
AC_rbf = (sum(group==check_rbf) - sum(group~=check_rbf))/length(group)*100  
AC_quadratic = (sum(group==check_quadratic) - sum(group~=check_quadratic))/length(group)*100      %
AC_mlp = (sum(group==check_mlp) - sum(group~=check_mlp))/length(group)*100      %

solution_linear = svmclassify(s_linear,B');            %
solution_polynomial = svmclassify(s_polynomial,B');            %
solution_rbf = svmclassify(s_rbf,B');        
solution_quadratic = svmclassify(s_quadratic,B');   
solution_mlp = svmclassify(s_mlp,B');            %


% solution = solution';


% st = find(check==1);
% 
% 
% sf = find(check==-1);    


[X_linear,Y_linear] = perfcurve(group,check_linear,1);
[X_polynomial,Y_polynomial] = perfcurve(group,check_polynomial,1);
[X_rbf,Y_rbf] = perfcurve(group,check_rbf,1);
[X_quadratic,Y_quadratic] = perfcurve(group,check_quadratic,1);
[X_mlp,Y_mlp] = perfcurve(group,check_mlp,1);


figure;
% plot(snr1,pd,'o-b',snr1,p3,'-.+m');
plot(X_linear,Y_linear,'m');
hold on
plot(X_polynomial,Y_polynomial,'b');
hold on
plot(X_rbf,Y_rbf,'c');
hold on
plot(X_quadratic,Y_quadratic,'r');
hold on
plot(X_mlp,Y_mlp,'g');
hold on
title ('Comparison of different SVM Kernel Functions')
legend('Linear','Polynomial','RBF','Quadratic','MLP');
xlabel('Probability of False Alarm');
ylabel('Probability of Detection');
