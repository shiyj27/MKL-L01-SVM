clc;
clear all;
load Synthetic;
% load WPBC;
% load ionosphere;
% X=x;
 L=10;
 %Har
% pars.C=32;
% pars.rho1=64;
% pars.rho2=64;
% pars.rho3=4;
%qita 
% pars.C=4;
% pars.rho1=4;
% pars.rho2=8;
% pars.rho3=8;
% Synthetic,ionosphere
pars.C=4;
pars.rho1=4;
pars.rho2=4;
pars.rho3=16;
[m1,n1]=size(X);
numRuns=30;
resiult=zeros(numRuns,1);
TACC=0;
acc=0;
for i=1:numRuns
     out=MKL01ADMM_decision(X,y,pars);  
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
             result(i,1)=TACC;
 end
 if TACC>0.9
     break;
 end
end
 
 d1 = 0.02; % Step size of the grid
[x1Grid,x2Grid] = meshgrid(min(Xt(:,1)):d1:max(Xt(:,1)),...
    min(Xt(:,2)):d1:max(Xt(:,2)));
xGrid = [x1Grid(:),x2Grid(:)]; 
 norms1=sum(xGrid'.^2);


sig=[0.1 0.2 0.3 0.5 0.7 1 1.2 1.5 1.7 2];
            [m2,~]=size(X1); 
            norms2=sum(X1'.^2);
            [mt,~]=size(xGrid);   
            g_l_test=zeros(mt,L);
            KK=cell(L,1);
            Dy1=diag(y1);

    for   l=1:L
       KK{l}=exp((-norms1'*ones(1,m2)-ones(mt,1)*norms2+2*xGrid*(X1'))/(2*sig(l)^2));
      g_l_test(:,l)=-d(l)*KK{l}*w;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    end 
      res2=0;
            for l=1:L
                res2=res2+d(l)*KK{l}*w;
            end




figure;
  h(1:2) = gscatter(X1(:,1),X1(:,2),y1);
   h(1:2) = gscatter(Xt(:,1),Xt(:,2),yt);
hold on
    
   contour(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)), [ 0 0],'k')
 
   title('Scatter Diagram with the Decision Boundary')
 legend({'-1','1','decision boundary'},'Location','northeast');
 grid on ;
% legend({'-1','1'},'Location','Best');
% hold off  
% 设置坐标轴标签和图例的字体大小
set(gca, 'FontSize', 12); % 12为你想要的字体大小
lgd = legend({'-1','1','decision boundary'},'Location','northeast');
set(lgd, 'FontSize', 12);

% 如果想要设置标题的字体大小
titleFontSize = 16;
title('Scatter Diagram with the Decision Boundary', 'FontSize', titleFontSize)


    
    
    
    
    
    
    
    
    
