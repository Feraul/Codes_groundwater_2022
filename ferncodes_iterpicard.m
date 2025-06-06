function [p,flowrate,flowresult,flowratedif]=ferncodes_iterpicard(M_old,RHS_old,...
    parameter,w,s,p_old,nflagno,wells,viscosity,Con,nflagc,wightc,sc,...
    dparameter,contnorm,SS,dt,h,MM,gravrate,source)
global nltol maxiter
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para itera��o Picard
step=0;
er=1;
pinterp_new=0;
while (nltol<er || nltol==er) && (step<maxiter)
    % atualiza itera��es
    step=step+1;
    
    %--------------------------------------------------------------------------
    %Solve global algebric system
    
    % calculo das press�es
    p_new = solver(M_old,RHS_old);
    
    %% plotagem no visit
    %S=ones(size(p_new,1),1);
    %ferncodes_postprocessor(p_new,S,step)
    [pinterp_new,]=ferncodes_pressureinterpNLFVPP(p_new,nflagno,w,s,Con,...
        nflagc,wightc,sc);
    %% Calculo da matriz global
    
    [M,I]=ferncodes_assemblematrixNLFVPP(pinterp_new,parameter,viscosity,...
        contnorm,SS,dt,h,MM,gravrate,nflagno);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    
    %Often it may change the global matrix "M"
    [M_new,RHS_new] = addsource(sparse(M),I,wells);
    % Often with source term
    [RHS_new]=sourceterm(RHS_new,source);
    %% Calculo do residuo
    
    R = norm(M_new*p_new - RHS_new);
    
    if (R0 ~= 0.0)
        er = abs(R/R0);
    else
        er = 0.0; %exact
    end
    errorelativo(step)=er;
    
    % atualizar
    M_old=M_new;
    RHS_old=RHS_new;
    
end
%--------------------------------------------------------------------------
%Solve global algebric system

% calculo das press�es
p = p_new;
%Message to user:
fprintf('\n Iteration number, iterations = %d \n',step)
fprintf('\n Residual error, error = %d \n',er)
%Message to user:
disp('>> The Pressure field was calculated with success!');
[pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflagno,w,s,Con,...
                                                         nflagc,wightc,sc);
%Get the flow rate (Diamond)
[flowrate,flowresult,flowratedif]=ferncodes_flowrateNLFVPP(p, pinterp,...
      parameter,viscosity,Con,nflagc,wightc,sc,dparameter,cinterp,gravrate);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end