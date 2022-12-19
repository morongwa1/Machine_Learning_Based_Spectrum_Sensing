

clear all;
close all;
randn( 'state', 0 );

option = 2; 



system_para.num_of_test = 1000;
system_para.record = zeros(1, system_para.num_of_test);

system_para.N = 100;
system_para.M = system_para.N;

system_para.step = 1;
system_para.H1_num = 1; %
system_para.H0_num = 1; %
system_para.u1_num = 1; %
system_para.u0_num = 0; %

system_para.ui1_num = zeros(1,system_para.M); % f
system_para.ui0_num = zeros(1,system_para.M); % f

system_para.N1 = 1; %fusion 
system_para.N0 = 1; %fusion 

system_para.P0 = 0.8;
system_para.P1 = 0.2;

system_para.PD_t = 0.90;
system_para.PF_t = 0.15;

system_para.Pd = randint(1,system_para.M,[500,800]); 
system_para.Pd = system_para.Pd/1000; 
system_para.Pf = randint(1,system_para.M,[200,500]); 
system_para.Pf = system_para.Pf/1000;

system_para.Qd_init = 0.7; % ¶ÔÓÚmajority rule£¬£¬
system_para.Qf_init = 0.3; 

system_para.Qd = ones(1,system_para.M); 
system_para.Qd = system_para.Qd * system_para.Qd_init;
system_para.Qf = ones(1,system_para.M); 
system_para.Qf = system_para.Qf * system_para.Qf_init;

system_para.QD = 0;
system_para.QF = 0;
system_para.PD = 0;
system_para.PF = 0;

system_para.alpha_d = 1;
system_para.alpha_f = 0;

system_para.u = zeros(1,system_para.M);

system_para.e1 = 0.001;
system_para.e2 = 0.03;


if(option == 0)
    CV(system_para);
elseif(option == 1)
    CV_with_alerting(system_para);
elseif(option == 2)
    majority(system_para);
elseif(option == 3)
    majority_with_alerting(system_para);
end

