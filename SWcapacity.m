function  [SW]=SWcapacity(h,theta_s,theta_r,alpha,p,q)
global elem 

for i=1:size(elem,1)

    if h(i)<0
        SW(i)=(theta_s-theta_r)*p*(1-1/p)*alpha*(-alpha*h(i))^(p-1)*(1+(-alpha*h(i))^p)^(1/p -2);

    else
       SW(i)=0;
    end

end