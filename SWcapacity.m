function  [Capw]=SWcapacity(h,theta_s,theta_r,alpha,pp,q)
global elem numcase  
% calculo da capicidade da agua
for i=1:size(elem,1)
if numcase==431
    if h(i,1)<0 || h(i,1)==0
        Capw(i,1)=-q*pp*(alpha^pp)*((theta_s-theta_r)/((1+(alpha*abs(h(i)))^pp)^(q +1)))*(abs(h(i))^pp)*(h(i)^-1);
    else
       Capw(i,1)=0;
    end
else
    if h(i,1)<0 || h(i,1)==0
        Capw(i,1)=(theta_s-theta_r)*pp*q*alpha*((-alpha*h(i))^(pp-1))*(1+...
            (-alpha*h(i))^pp)^((1/pp) -2);
    else
       Capw(i,1)=0;
    end
end
end