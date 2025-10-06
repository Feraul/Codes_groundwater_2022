function [theta_s,theta_r,alpha,pp,q,h_init,dt]=preRE
global elem numcase centelem

switch numcase
    case 431
        
        %------------------------------------------------------------------
        theta_s=0.363;
        theta_r=0.186;
        alpha=0.01;
        pp=1.53;
        % K=0.0001
        % h=0  CC. Dirichlet upper
        q=1-(1/pp);
        h_init=-800*ones(size(elem,1),1);
        dt=195;

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
    case 433
         theta_s=0.43;
         theta_r=0.078;
         alpha=0.036;
         pp=1.56;
         q=1-(1/pp);
         % K= 24.96
         %flux 20  CC. Neumann
         h_init=-51.3949*ones(size(elem,1),1);
         dt=0.05;
    case 434
         theta_s=0.3;
         theta_r=0.01;
         alpha=0.033;
         pp=4.1;
         q=1-(1/pp);
        for i=1:size(elem,1)
            h_init(i,1) =65-centelem(i,2);
        end
         dt=5;
    case 435
        %------------------------------------------------------------------
        theta_s=0.43;
        theta_r=0.078;
        alpha=0.036;
        pp=1.56;
        % K=0.0001
        % h=20  CC. Dirichlet
        q=1-(1/pp);
        for i=1:size(elem,1)
            a =10-centelem(i,2);
            if a>0
                h_init(i,1)=10;
            else
                h_init(i,1)=-90;
            end
        end
        dt=0.03;


end
end