function K=MMKLkernel(xapp,InfoKernel,nnbk,xtest,beta)

if nargin<4
beta=[];
xtest=xapp;

 for k=1:nnbk

        Kr=createkernel(xapp(:,InfoKernel(k).variable),InfoKernel(k).kerneloption, xtest(:,InfoKernel(k).variable));
        % Kr=Kr*Weight(k);
%         if options.efficientkernel
%             Kr=build_efficientK(Kr);
%         end;

        K(:,:,k)=Kr;


 end
else
    ind=find(beta);
    K=zeros(size(xapp,1),size(xtest,1));
    for i=1:length(ind)
        k=ind(i); 
        Kr=createkernel(xapp(:,InfoKernel(k).variable),InfoKernel(k).kerneloption, xtest(:,InfoKernel(k).variable));
        % Kr=Kr*Weight(k);
        K=K+ Kr*beta(k);
    end

end




