function [x,y]=pd_selfAdpt(x1,x2,snr)
n1=256;                 %��������
n2=20;                  %��ͼ����
time=1000;              %������
dB=snr;
Pf1=linspace(x1,x2,n2); %�龯��Χ
w=ones(1,20);        %Ȩ��
L=length(Pf1);
for i = 1:length(Pf1)                                % ����
    dB_binary= power(10,dB/10);                      %�����dbת��ʮ����
    Lambda(i)= (n1+sqrt(2*n1)*qfuncinv(Pf1(i)));     %Ԥ���龯��������_TO Pd_theory
    Nd=0;
       for n=1:time              
            t=(0:n1-1)/100; 
            Apt=sqrt(dB_binary);
            x=Apt*cos(2*300*pi*t).*cos(2*pi*1000*t); %�����ź�     
            noise=randn(1,n1);                       %��˹�����������ֵ������1
            y=x+noise;                               %�����ź�=�ź�+����
            Y=fft(y,n1);       
            signalPower = (sum(x).^2)/length(x);
            T_simu(i)=sum(abs(Y).^2)/n1;
                if  T_simu(i) > Lambda(i)
                    Nd=Nd+1;
                end
       end
   Pd_simu(i)=Nd/n;
   w(i+1)=w(i)*Pd_simu(i)*sqrt(w(i)*Pd_simu(i));
   y=Pd_simu;
   x=w;
end
%figure(1)
%plot(Pf1,Pd_simu,'dk',Pf1,Pd_theory,'-k');
%xlabel('Pf')
%ylabel('Pd_single')
%grid on