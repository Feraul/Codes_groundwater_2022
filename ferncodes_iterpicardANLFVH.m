function [p,flowrate,flowresult]=ferncodes_iterpicardANLFVH(M_old,RHS_old,...
    parameter,p_old,nflagface,wells,mobility,weightDMP,w,s,nflagno,contnorm)
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);


%% Acelerador de Anderson
%p_old=M_old\RHS_old;
 [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));

%[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
[p_old,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
%[p_new]=ferncodes_andersonacc(p_old,parameter,w,s,nflagno,mobility,wells,...
%    weightDMP,nflagface);
[p,erro,iter]=ferncodes_andersonacc2(p_old,1e-6,parameter,w,s,...
    nflagface,weightDMP,wells,mobility,R0,contnorm);
% se o campo de ressão é negativo ele coloca zero
%Message to user:
fprintf('\n Iteration number, iterations = %d \n',iter)
fprintf('\n Residual error, error = %d \n',erro)
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpHP(p,nflagface,parameter,weightDMP,mobility);
[flowrate,flowresult]=ferncodes_flowrateNLFVH(p, pinterp, parameter,mobility);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
end
