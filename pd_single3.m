function y=pd_single3(x1,x2,snr)
n1=256;%��������
n2=20;%��ͼ����
Pf1=linspace(x1,x2,n2);
dB=snr; 
L=length(Pf1);
for i = 1:L      % ����
    dB_binary= power(10,dB/10);%�����dbת��ʮ����
    Lambda(i)= (n1+sqrt(2*n1)*qfuncinv(Pf1(i)));%Ԥ���龯��������_TO Pd_theory
       Nd=0;  
       for j=1:1000  %������
            t=(0:n1-1)/100;  %
            Apt=sqrt(dB_binary);
            x=Apt*cos(2*300*pi*t).*cos(2*pi*1000*t); %�����ź�     
            noise=randn(1,n1);%��˹�����������ֵ������1
            y=x+noise;%�����ź�=�ź�+����
            Y=fft(y,n1);       
            signalPower = (sum(x).^2)/length(x);
            T_simu(i)=sum(abs(Y).^2)/n1;
                if  T_simu(i) > Lambda(i)
                    Nd=Nd+1;
                end
          end
   Pd_simu(i)=Nd/j;
   y=Pd_simu;
   Pd_theory(i) = 0.5*erfc((Lambda(i)-n1-n1*dB_binary)./sqrt(2*(2*n1+4*n1*dB_binary)));
end
%figure(1)
%plot(Pf1,Pd_simu,'dk',Pf1,Pd_theory,'-k');
%xlabel('Pf')
%ylabel('Pd_single')
%grid on