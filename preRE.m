function [theta_s,theta_r,alpha,pp,q,h_init,dt]=preRE
global elem
theta_s=0.43;
theta_r=0.045;
alpha=1.45;
pp=2.68;
q=1-(1/pp);

h_init=-100*ones(size(elem,1),1);
dt=0.05;

end