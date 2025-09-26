function [theta_s,theta_r,alpha,pp,q,h_init,dt]=preRE
global elem numcase centelem

switch numcase
    case 431
        theta_s=0.43;
        theta_r=0.045;
        alpha=1.45;
        pp=2.68;
        q=1-(1/pp);
        h_init=-100*ones(size(elem,1),1);
        dt=0.05;
    case 432
        theta_s=0.396;
        theta_r=0.131;
        alpha=0.423;
        pp=2.06;
        q=1-(1/pp);
        for i=1:size(elem,1)
            h_init(i,1) =1-centelem(i,2);
        end
        %dt=3.3333;
        dt=0.021;
end
end