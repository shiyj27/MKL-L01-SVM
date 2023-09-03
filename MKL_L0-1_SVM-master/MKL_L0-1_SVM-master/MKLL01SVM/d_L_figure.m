clc;
clear all;
load WPBC;
% load ionosphere;
% X=x;
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
% fprintf('------------------------------------------------------------------------------\n');
%   fprintf('    C    rho1   rho2   rho3    acc    TACC    NSV   iter   time\n');
% fprintf('------------------------------------------------------------------------------\n');
%  for i=2:1:8
%    pars.C=1;
% %     pars.rho1=6;
% %     pars.rho2=1;
% %     pars.rho3=6;
% pars.C=2^i;
% for g=2:1:8
%     pars.rho1=2^g;
%     for h=2:1:8
%     pars.rho2=2^h;
%     for j=2:1:8
%     pars.rho3=2^j;
numRuns=2;
resiult=zeros(numRuns,1);
TACC=0;
% acc=0;
for i=1:numRuns
    pars.numRuns=i;
     out=MKL01ADMM_d(X,y,pars);  
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
%             Xt=X;
%   X1=X;
%   yt=y;
%   y1=y;
%             Xt=X1;
%             yt=y1;
%             norms1=sum(Xt'.^2)';
%             [m2,~]=size(X1);
%             norms2=sum(X1'.^2);
%             [mt,~]=size(Xt);   
%             g_l_test=zeros(mt,L);
%             KK=cell(L,1);
%             Dy1=diag(y1);
%   
%             for l=1:L
%             KK{l}=exp((-norms1*ones(1,m2)-ones(mt,1)*norms2+2*Xt*(X1'))/(2*sig(l)^2));
%             end
%             for l=1:L
%                 res1=res1+d(l)*KK{l}*w;
%             end
%             
%          
%              TACC=1-sum(abs(sign((res1+b))-yt))/(2*length(yt));
%  end
% % %             
% % %             
     % fprintf('|%5.2f| %5.2f| %5.2f|%5.4f |%5.4f|%5.4f |%5.2f| %4d|%5.2fsec |\n',...
     %       pars.C,pars.rho1,pars.rho2,pars.rho3,acc,TACC,nsv,iter,time)
% %        fprintf('|%5.2f| %5.2f| %5.2f|%5.4f |%5.4f|%5.2f| %4d|%5.2fsec |\n',...
% %            pars.C,pars.rho1,pars.rho2,pars.rho3,TACC,nsv,iter,time)
 end
% end

%     end  
%     end
%     end
% end
% if nnz(result)==numRuns
%     Meanacc=mean(result)
% end
% %     
% %     
% %    
% if TACC>0.90
%     break;
% end
%     end
%     if TACC>0.7
%         break;
%     end
% %     end
%     if TACC>0.7
%         break;
%     end
% % end
% if TACC>0.7
%         break;
% end
% end
% 
% % Xt=X1;
% % yt=y1;
% %     if out.flag==3
% %             b=out.b;
% %            d=out.d;
% %             iter=out.iter;
% %             nsv=out.nsv;
% %             time=out.time;
% %             lam=out.lam;
% % %            acc=out.acc;       
% %          norms1=sum(Xt'.^2)';
% %             [m2,~]=size(X1);
% %             norms2=sum(X1'.^2);
% %             [mt,~]=size(Xt);   
% %             g_l_test=zeros(mt,L);
% %             KK=cell(L,1);
% %             Dy1=diag(y1);
% %        for l=1:L
% %            if d(l)<1e-4
% %                  d(l)=0;
% %            end  
% % %           d= [0.0021  0.4547 0  0.5290 0.0113 0 0 0 0 ]
% % %          d= [0  0 0.4547 0  0.5290 0 0 0 0 0 ]
% %         KK{l}=exp((-norms1*ones(1,m2)-ones(mt,1)*norms2+2*Xt*(X1'))/(2*sig(l)^2));
% %         g_l_test(:,l)=-d(l)*KK{l}*(y1.*lam);
% %         end      
% %                TACC=1-sum(abs(sign(sum(g_l_test,2)+b)-yt))/(2*length(yt));
% %                fprintf('|%5.2f| %5.2f| %5.2f | %5.2f|%5.4f |%5.4f | %4d| %4d|%5.2fsec |\n',...
% %            pars.C,pars.rho1,pars.rho2,pars.rho3,acc,TACC,nsv,iter,time)
% %        fprintf('|%5.2f| %5.2f| %5.2f|%5.4f |%5.2f | %4d|%5.2fsec |\n',...
% %            pars.C,pars.rho1,pars.rho2,pars.rho3,nsv,iter,time)
% % %       
%      %end
%  d1 = 0.02; % Step size of the grid
% [x1Grid,x2Grid] = meshgrid(min(Xt(:,1)):d1:max(Xt(:,1)),...
%     min(Xt(:,2)):d1:max(Xt(:,2)));
% xGrid = [x1Grid(:),x2Grid(:)]; 
%  norms1=sum(xGrid'.^2);
% 
% % %  [m2,~]=size(xGrid); 
% % %             norms2=sum(xGrid'.^2);
% % %             [mt,~]=size(xGrid);   
% % %             g_l_test=zeros(mt,L);
% % %             KK=cell(L,1);
% % %             Dy1=diag(y1);
% % %  for   l=1:L
% % %        KK{l}=exp((-norms1'*ones(1,m2)-ones(mt,1)*norms2+2*xGrid*(xGrid'))/(2*sig(l)^2));
% % %       g_l_test(:,l)=-d(l)*KK{l}*(y1.*lam);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
% % %     end
% % 
% sig=[0.1 0.2 0.3 0.5 0.7 1 1.2 1.5 1.7 2];
%             [m2,~]=size(X1); 
%             norms2=sum(X1'.^2);
%             [mt,~]=size(xGrid);   
%             g_l_test=zeros(mt,L);
%             KK=cell(L,1);
%             Dy1=diag(y1);
% 
%     for   l=1:L
%        KK{l}=exp((-norms1'*ones(1,m2)-ones(mt,1)*norms2+2*xGrid*(X1'))/(2*sig(l)^2));
%       g_l_test(:,l)=-d(l)*KK{l}*w;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%     end 
%       res2=0;
%             for l=1:L
%                 res2=res2+d(l)*KK{l}*w;
%             end
% 
% 
% 
% 
% figure;
%   h(1:2) = gscatter(X1(:,1),X1(:,2),y1);
%    h(1:2) = gscatter(Xt(:,1),Xt(:,2),yt);
% hold on
%     % plot(X1(out.T,1),X1(out.T,2),'ko','MarkerSize',10) ;
% %    contour(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)),[ 1  1],'b'); 
% % contour(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)),[ 0  0],'k'); 
%    % contour(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)), [-1 0 1],'r'); 
%    contour(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)), [ 0 0],'k')
%  %    % mesh(x1Grid,x2Grid,reshape(sum(g_l_test,2)+b,size(x1Grid)));
%    title('Scatter Diagram with the Decision Boundary')
%  legend({'-1','1','decision boundary'},'Location','northeast');
%  grid on ;
% % legend({'-1','1'},'Location','Best');
% % hold off  
% % 设置坐标轴标签和图例的字体大小
% set(gca, 'FontSize', 12); % 12为你想要的字体大小
% lgd = legend({'-1','1','decision boundary'},'Location','northeast');
% set(lgd, 'FontSize', 12);
% 
% % 如果想要设置标题的字体大小
% titleFontSize = 16;
% title('Scatter Diagram with the Decision Boundary', 'FontSize', titleFontSize)

% % % %     
% % %     
% % %     
% % %     
% % %     
% %     
    
    
    
    
    
    
    
    
    
    
