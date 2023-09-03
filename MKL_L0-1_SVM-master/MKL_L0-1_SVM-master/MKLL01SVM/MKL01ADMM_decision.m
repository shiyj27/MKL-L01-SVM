function  out= MKL01ADMM_decision(x,y,pars)


verbose=1;
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
kernelt={'gaussian' 'gaussian'};
kerneloptionvect={[0.1 0.2 0.3 0.5 0.7 1 1.2 1.5 1.7 2] [0.1 0.2 0.3 0.5 0.7 1 1.2 1.5 1.7 2]};
variablevec={'all' 'single'};
nbiter=1;
ratio=0.7;
classcode=[1 -1];
[nbdata,dim]=size(x);
nbtrain=floor(nbdata*ratio);
% rand('state',0);
  
    [xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(x, y, nbtrain,classcode);%划分训练集和测试集，并记录训练和测试的索引
    [xapp,xtest]=normalizemeanstd(xapp,xtest);%对输入数据进行标准化处理，均值设置为0，标准差设置为1.
    [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
    kernel={'gussian'};
     [nbk,Infokernel]=Unit1(kernel,kerneloptionvec,variableveccell);
     nnbk=nbk-1;
   K1=MMKLkernel(xapp,Infokernel,nnbk);
K=cell(1,10);
    for i=1:10
        K{i}=K1(:,:,i);
    end
L=10;
x=xapp;  
y=yapp;
TACC=0;
[m,n]=size(xapp);
if nargin<3;               pars  = [];                             end
if isfield(pars,'maxit');  maxit = pars.maxit; else; maxit = 2000;  end
if isfield(pars,'sig');  sig = pars.sig; else; sig = 1;    end
if isfield(pars,'rho1');  rho1 = pars.rho1; else; rho1 = 1;    end
if isfield(pars,'rho2');  rho2 = pars.rho2; else; rho2 = 1;    end
if isfield(pars,'rho3');  rho3 = pars.rho3; else; rho3 = 1;    end
if isfield(pars,'tol');    tol   = pars.tol;   else; tol   = 1e-3; end
if isfield(pars,'C');      C     = pars.C;     else; C     = 100;    end
if isfield(pars,'w');      w     = pars.w;     else; w =zeros(m,1);   end
if isfield(pars,' theta'); theta = pars.theta; else; theta=zeros(L,1);  end
if isfield(pars,'alpha');  alpha = pars.alpha; else; alpha=0;    end
if isfield(pars,'lam');     lam  = pars.lam;   else; lam=zeros(m,1);   end
if isfield(pars,'b');         b  = pars.b;     else; b=1;       end
if isfield(pars,'z');          z = pars.z;     else; z=zeros(L,1);  end
if isfield(pars,'d');          d= pars.d;     else; d=1/L*ones(L,1);  end
Dy=diag(yapp);%求对角阵
y=yapp;
X=xapp;
beta=zeros(1,8);
lamrho1=lam/rho1;
flag=0;%定义一个标签变量，后续flag变化时候进入此（与u的工作集T的更新有关）
flag1=0;
to=tic; %程序计时器
T      = GetT(y,m,n);%数据集的划分，选取正负标签各一半作为初始的工作集
TS=(1:L)';  %z的初始工作集
TSb=[];
ss=(d+theta/rho2); %与z的更新有关
u=zeros(m,1);   % u的初始化
t=1:L;   
% fig = figure;
% FirstDvalues=zeros(1,200);
% SecondDvalues=zeros(1,200);
% ThirdDvalues=zeros(1,200);
% maxit=2000;
for iter=1:1:maxit    %开始迭代
    t=1:L;
%     update T   （工作集的更新，初始工作集为给定的，当满足条件时进入此次循环）  
    if flag==1 
       con = sqrt(2*C/rho1);
       T   = find( s>0 & s <= con);   %找出满足0<s<cons 的下标T,对T更新
    end
    nT=nnz(T);%当结束循环时候（收敛的时候），对应的支持向量的个数
%%
%update u
sum2=0;
for l=1:L 
    sum2=sum2+d(l)*K{l};
end
s=ones(m,1)-Dy*sum2*w-b*y-lam/rho1 ;%prox的变量（与u的更新有关）
temp1=u;
u = Prox(s,C,rho1,T);   %计算temp1-u来判断是否收敛
%update w
temp2=w;
w=solutionw(L,rho1,u,d,K,Dy,b,y,lam,m);
w;
% %update b  
temp3=b;
b=solutionB(y,u,rho1,m,lam,L,Dy,d,K,w);  
% %update TS
    if flag1==1 
       TS   = find(ss>0);
       TSb=find(ss<=0);
    end
% %update z
temp4=z;
ss=(d+theta/rho2);
z=mapp(TS,ss);

% % %  update d
temp5=d;
% a1=zeros(L);
% for i=1:L
%     for j=1:L
% %         a1(i,j)=rho1*w'*K(:,:,i)*K(:,:,j)*w;
%         a1(i,j)=rho1*w'*K{i}*K{j}*w;      
%     end 
% end
% b1=rho2*eye(L,L);
% c1=rho3*ones(L,L);
% aaa=(a1+b1+c1);
% myvector=zeros(L,1);
% for l=1:L
%     myvector(l)=-0.5*w'*K{l}*w-lam'*Dy*K{l}*w-rho1*w'*K{l}*Dy*(u+b*y-ones(m,1));
% %     myvector(l)=-0.5*w'*K(:,:,l)*w-lam'*Dy*K(:,:,l)*w-rho1*w'*K(:,:,l)*Dy*(u+b*y-ones(m,1));
% end
% bbb=myvector-theta+rho2*z+(rho3-alpha)*ones(L,1);
% d=aaa\bbb;
% d;
for l=1:L
res=0;
    for t1=1:L-1
      k={K{t~=l}};
        d1=d(t~=l);
        res=res+d1(t1)*k{t1};  
    end
d(l)=-(0.5*w'*K{l}*w+lam'*Dy*K{l}*w+...
rho1*(Dy*K{l}*w)'*(u+b*y-ones(m,1))+rho1*w'*K{l}*(res)*w...
 +theta(l)-rho2*z(l)+alpha+rho3*(sum(d(t~=l))-1))/(rho1*w'*K{l}*K{l}*w+rho2+rho3);
end

 temp6=theta;
theta=solutionTheta(d,z,theta,rho2,TS);

% %update alpha
temp7=alpha;
alpha=solutionAlpha(alpha,rho3,L,d);
dd=SimplexProj(d');
d=dd';
d(TSb)=0;
dmaxthree=sort(d,'descend');
FirstDvalues(iter)=dmaxthree(1);
SecondDvalues(iter)=dmaxthree(2);
ThirdDvalues(iter)=dmaxthree(3);
% % %update lamda
temp8=lam;
lam=solutionLamda(lam,u,Dy,rho1,T,m,b,y,L,d,K,w);

% % %stopping crterion
 sumK=0;
for l=1:L
    sumK=sumK+d(l)*K{l};
end
s=ones(m,1)-Dy*sumK*w-b*y-lam/rho1 ;%prox的变量（与u的更新有关）
 ss=(d+theta/rho2);
 beta(1)=norm(u-temp1);
 beta(2)=norm(w-temp2);
 beta(3)=norm(b-temp3);
 beta(4)=norm(z-temp4);
 beta(5)=norm(d-temp5);
 beta(6)=norm(theta-temp6);
 beta(7)=norm(alpha-temp7);
 beta(8)=norm(lam-temp8);
 error=max(beta);
 CPU=toc(to);
  ACC=0;
  TACC=ACC;
 if error<tol 
     flag=3;
     res2=0;
 for l=1:L
    res2=res2+d(l)*K{l}*w;
 end       
         ACC=mean(sign(res2+b)==yapp);
     X1=xapp;
     y1=yapp;
        Xt=xtest;
            yt=ytest;
 sig=[0.1 0.2 0.3 0.5 0.7 1 1.2 1.5 1.7 2];
            Kt=MMKLkernel(xapp,Infokernel,nnbk,xtest,d);
             ypred=Kt'*w+b;
      TACC=mean(sign(ypred)==ytest);
      break;
 end

if iter>5 

         flag=1;
         flag1=1;
     end
end





out.iter = iter;
out.time = CPU;
out.w  = w;
out.u    = u;
out.z=z;
out.b=b;
out.theta=theta;
out.d=d;
out.alpha=alpha;
out.lam  = lam;
out.T=T;
out.nsv  = nT;
out.acc  = ACC;
out.error= error;
out.flag = flag;
out.K=K;
out.TACC=TACC;
out.X1=xapp;
out.y1=yapp;
out.Xt=xtest;
out.yt=ytest;

end

function theta=solutionTheta(d,z,theta,rho2,TS)
  % theta=theta+rho2*(d-z);
theta(TS)=theta(TS)+rho2*(d(TS)-z(TS));
end
function lam=solutionLamda(lam,u,Dy,rho1,T,m,b,y,L,d,K,w)
lamT   = lam(T);
lam    = zeros(m,1);
sumA1d=0;
 for l=1:L
     sumA1d=sumA1d+Dy*(d(l)*K{l});
% sumA1d=sumA1d+Dy*(d(l)*K(:,:,l));
 end
 r=u+sumA1d*w+b*y-ones(m,1);
lamT   = lamT + rho1*r(T); 
lam(T) = lamT;
end
function w=solutionw(L,rho1,u,d,K,Dy,b,y,lam,m)
sum1=0;
for l=1:L
  sum1=sum1+d(l)*K{l};  
% sum1=sum1+d(l)*K(:,:,l); 
end
w=-(eye(m,m)+rho1*sum1)\(Dy*(lam+rho1*(u+b*y-ones(m,1)))); 
end
 function b=solutionB(y,u,rho1,m,lam,L,Dy,d,K,w)
 sumAd=0;
 for l=1:L
     sumAd=sumAd+Dy*(d(l)*K{l});
% sumAd=sumAd+Dy*(d(l)*K(:,:,l));
 end
b=-(y'*(lam+rho1*(u+sumAd*w-ones(m,1))))/(m*rho1);
 end
function alpha=solutionAlpha(alpha,rho3,L,d)
alpha=alpha+rho3*(ones(L,1)'*d-1);
end
function u= Prox(s,C,rho1,T)
if isempty(T)
    cons=sqrt(2*C/rho1);
    T=find(s>0 &s<=cons);
    u=s;
    u(T)=0;
else 
    u=s;
    u(T)=0;
end
end 
function z=mapp(TS,ss)
if isempty(TS)
    TS=(find (ss>0));
    z=zeros(size(ss));
    z(TS)=ss(TS);
   
else
    z=zeros(size(ss));
    z(TS)=ss(TS); 
end
end
function T = GetT(y,m,n)
if  m>n   
s0      = ceil(n*(log(m/n))^2); 
T1      = find(y==1);  nT1= nnz(T1);
T2      = find(y==-1); nT2= nnz(T2);
if  nT1 < s0
    T  = [T1; T2(1:(s0-nT1))];   
elseif nT2 < s0
    T  = [T1(1:(s0-nT2)); T2]; 
else 
    T  = [T1(1:ceil(s0/2)); T2(1:(s0-ceil(s0/2)))];    
end
T      = sort(T(1:s0));
else
    T  = 1:length(y);
end
end
