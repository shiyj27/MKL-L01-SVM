



function K=createkernel(xapp,sigma,xtest)

metric = diag(1./sigma.^2);
    ps = xapp*metric*xtest'; 
    [nps,pps]=size(ps);
    normx = sum(xapp.^2*metric,2);
    normxsup = sum(xtest.^2*metric,2);
    ps = -2*ps + repmat(normx,1,pps) + repmat(normxsup',nps,1) ; 
    
    
    K = exp(-ps/2);

 

end