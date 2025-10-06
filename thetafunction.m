

function  theta=thetafunction (h,theta_s,theta_r,alpha,p,q,iterinicial)
global elem numcase


for i=1:size(elem,1)
    if numcase==431 || numcase==435
        if h(i,1)<0

            theta(i,1)= theta_r +((theta_s -theta_r)/(1+abs(alpha*h(i,1))^p)^q);

        else
            theta(i,1)=theta_s;
        end
    elseif numcase==432 || numcase==434
        if h(i,1)<0 || h(i,1)==0
            theta(i,1)= theta_r +(theta_s -theta_r)*(1/(1+(-alpha*h(i,1))^p))^q;
        else
            theta(i,1)=theta_s;
        end

    elseif numcase==433 
        if h(i,1)<0

            theta(i,1)= theta_r +((theta_s -theta_r)/(1+abs(alpha*h(i,1))^p)^q);

        else
            theta(i,1)=theta_s;
        end

    end
end