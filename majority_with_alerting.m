%   感谢亲亲使用此代码，此代码解决您的问题了吗~(@^_^@)~
%   没解决的话告诉亲亲一个好消息，我这里可以1分钱帮助代码改错，还提供1分钱成品代码(′▽`〃)哦~
%   登录淘宝店铺“大成软件工作室”便可领取
%   是的，亲亲真的没有看错，挠破头皮的问题真的1分钱就可以解决了\(^o^)/YES!
%   小的这就把传送门给您，记得要收藏好哦(づ￣3￣)づ╭～
%   传送门：https://item.taobao.com/item.htm?spm=a1z10.1-c.w4004-15151018122.5.uwGoq5&id=538759553146
%   如果传送门失效，亲亲可以来店铺讨要，客服MM等亲亲来骚扰哦~(*/ω╲*)
function majority_with_alerting(system_para);

% 在本方案中，我们忽略对错警概率的要求，只确保检测概率大于某门限即可。


n = 0;

QM_num = 0; % 记录漏警的次数，即PU发出alerting的次数

while(n < system_para.num_of_test)
    n = n + 1;
    k = floor ( system_para.M / 2 ) - 1;
    
    % step1 将alpha由高到低排序
    system_para.alpha = system_para.alpha_d * system_para.Qd - system_para.alpha_f * system_para.Qf;
    [system_para.alpha, IX] = sort( system_para.alpha, 'descend');

    % step2 generate u(i), i = 1,..,M
    PU_existance = randsrc(1,1,[-1,1;system_para.P0,system_para.P1]);
    
    for i = 1:1:system_para.M
        if(PU_existance == 1)
            system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pd(i),1 - system_para.Pd(i)]);
        else
            system_para.u(i) = randsrc(1,1,[1,-1;system_para.Pf(i),1 - system_para.Pf(i)]);
        end
    end
    system_para.u = system_para.u(IX); % 将u按关键字alpha排序
    
    % step3

     % sub_step1 judgement and update PD, PF
     U = 0;
    
     for i = 1:1:system_para.M
         U = U + system_para.u(i);
     end
    
     threshold = 2 * k - system_para.M;
     if(U >= threshold) % 0比threshold效果好
         u_global = 1;
     else
         u_global = -1;
     end
     
     if(u_global == 1 && PU_existance == 1)
         system_para.H1_num = system_para.H1_num + 1;
         system_para.u1_num = system_para.u1_num + 1;
         system_para.N1 = system_para.N1 + 1;
         
     elseif(u_global == -1 && PU_existance == 1) % 漏警
         system_para.H1_num = system_para.H1_num + 1;
         QM_num = QM_num + 1;
         system_para.N0 = system_para.N0 + 1;    
                
     elseif(u_global == 1 && PU_existance == -1) % 虚警
         system_para.H0_num = system_para.H0_num + 1;
         system_para.u0_num = system_para.u0_num + 1;
         system_para.N1 = system_para.N1 + 1;  
         
     elseif(u_global == -1 && PU_existance == -1)
         system_para.H0_num = system_para.H0_num + 1;
         system_para.N0 = system_para.N0 + 1;                
     end
      
     system_para.PD = system_para.u1_num / system_para.H1_num;
     system_para.PF = system_para.u0_num / system_para.H0_num; 
     
     % sub_step2: update Qd and Qf
     for i = 1:1:system_para.M
         if( system_para.u(i) == 1 && u_global == 1 )
             system_para.Qd(i) = ( system_para.Qd(i) * (system_para.N1 - 1) + 1 ) / system_para.N1;            
             
         elseif( system_para.u(i) == -1 && u_global == 1 )
             system_para.Qd(i) = ( system_para.Qd(i) * (system_para.N1 - 1) )/ system_para.N1;
             
         elseif( system_para.u(i) == 1 && u_global == -1 )
             system_para.Qf(i) = ( system_para.Qf(i) * (system_para.N0 - 1) + 1 )/ system_para.N0;       
             
         elseif( system_para.u(i) == -1 && u_global == -1 )
             system_para.Qf(i) = ( system_para.Qf(i) * (system_para.N0 - 1) )/ system_para.N0;    
         end
     end 

     % sub_step4: update QD and QF    
     system_para.QD = system_para.u1_num / system_para.N1;
     system_para.QF = system_para.u0_num / system_para.N0;

     
     display( 'program is running....just finished a test' );
     log.n=n;
     log.M=system_para.M;
     log.u1_num=system_para.u1_num;
     log.QM_num=QM_num;
     log.u0_num=system_para.u0_num;
     log.N1=system_para.N1;
     log.N0=system_para.N0;
     log.H1_num=system_para.H1_num;
     log.H0_num=system_para.H0_num;
     log.QD=system_para.QD;
     log.QF=system_para.QF;
     log.PD=system_para.PD;
     log.PF=system_para.PF;
     log.alpha=system_para.alpha(1:10);
     log.u=system_para.u(1:10);
     log.PU_exist=PU_existance;
     log.u_global=u_global;
     log
     system_para.record(n) = system_para.M;
     
     % sub_step5: update M
     if(system_para.QD > system_para.PD_t + system_para.e1)
         system_para.M = system_para.M - system_para.step;
     elseif(system_para.QD < system_para.PD_t - system_para.e2)
         system_para.M = system_para.M + system_para.step;
     end 
     
     if(system_para.M > system_para.N)
         system_para.M = system_para.N;
     elseif(system_para.M <= 1)
         system_para.M = 1;
     end
     

 
end % end of while

display( 'running totally finished..............' );
plot(system_para.record);
title('M vs n');
 