clc;
clear all;
load WPBC;
% load ionosphere;
% load Sonar;
%  X=x;  %iono needs start it
 L=10;
 %Har
% pars.C=32;
% pars.rho1=64;
% pars.rho2=64;
% pars.rho3=4;
%qita 
pars.C=4;
pars.rho1=4;
pars.rho2=8;
pars.rho3=8;
% Xynew,ionosphere
% pars.C=4;
% pars.rho1=4;
% pars.rho2=4;
% pars.rho3=16;
[m1,n1]=size(X);
fprintf('------------------------------------------------------------------------------\n');
  fprintf('    C    rho1   rho2   rho3    acc    TACC    NSV   iter   time\n');
fprintf('------------------------------------------------------------------------------\n');

numRuns=30;
result=zeros(numRuns,1);
result2=zeros(numRuns,1);
result3=zeros(numRuns,1);
result4=zeros(numRuns,1);
  for i=1:numRuns
      TACC=0;
      % d=0;
     out=MKL01ADMM(X,y,pars);  
 if out.flag==3
            b=out.b;
            pars.b=b;
            d=out.d;
           pars.d=d;
           w=out.w;
           X1=out.X1;
           y1=out.y1;
           Xt=out.Xt;
           yt=out.yt;
           pars.w=w;
           pars.u=out.u;
           pars.z=out.z;
           pars.theta=out.theta;
           pars.alpha=out.alpha;
           pars.lam=out.lam;
            iter=out.iter;
            nsv=out.nsv;
            time=out.time;
            lam=out.lam;
            K=out.K;
            res1=0;
             acc=out.acc;
             TACC=out.TACC;
              
  
     fprintf('|%5.2f| %5.2f| %5.2f|%5.4f |%5.4f|%5.4f |%5.2f| %4d|%5.2fsec |\n',...
           pars.C,pars.rho1,pars.rho2,pars.rho3,acc,TACC,nsv,iter,time)

 end
 result(i,1)=TACC;
 result2(i,1)=nsv;
 result3(i,1)=time;
  result4(i,1)=nnz(d);
  end
  acc_junzhi=mean(result(11:numRuns))
  acc_std=std(result(11:numRuns))
  Nsv_junzhi=mean(result2(11:numRuns))
  Nsv_std=std(result2(11:numRuns))
  time_junzhi=mean(result3(11:numRuns))
  time_std=std(result3(11:numRuns))
 d_junzhi=mean(result4(11:numRuns))
 d_std=std(result4(11:numRuns))


    
    
    
    
    
    
    
    
    
