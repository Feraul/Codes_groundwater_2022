

function  theta=thetafunction(h,theta_s,theta_r,alpha,pp,q)
global elem 


for i=1:size(elem,1)

    if h(i,1)<0

        theta(i,1)= theta_r +((theta_s -theta_r)/(1+abs(alpha*h(i,1))^pp)^q);

    else
        theta(i,1)=theta_s;
    end
end