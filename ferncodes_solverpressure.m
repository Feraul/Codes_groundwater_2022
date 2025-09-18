%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in poorus media
%Modified: Fernando Contreras, 2021
function [p,flowrate,flowresult,flowratedif] = ferncodes_solverpressure( viscosity,...
    wells,Hesq,Kde,Kn,Kt,Ded,nflag,nflagface,weight,s,Con,Kdec,Knc,Ktc,Dedc,...
    nflagc,wc,sc,SS,dt,h,MM,gravrate,P,kmap,tempo,N,h_kickoff,source,...
    theta_s,theta_r,alpha,pp,q)
global numcase acel elem interptype

%--------------------------------------------------------------------------
%Solve global algebric system: pressure or hydraulic head
if numcase==331 || numcase==431
    
    %----------------------------------------------------------------------
    % Montagem da matriz global
    [M_old,I_old,elembedge] = ferncodes_globalmatrix(weight,s,Kde,Ded,Kn,...
        Kt,Hesq,viscosity,nflag,nflagface,SS,dt,h,MM,gravrate,theta_s,...
        theta_r,alpha,pp,q,h_kickoff);
    %----------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"

    %Often it may change the global matrix "M" with wells
    [M_old,I_old] = addsource(sparse(M_old),I_old,wells);

    % Often with source term
    [I_old]=sourceterm(I_old,source);
    parameter=0; dparameter=0; contnorm=0; weightc=0;
    %----------------------------------------------------------------------
    if strcmp(acel,'FPI')
        [p,flowrate,flowresult,]=ferncodes_iterpicard(M_old,I_old,...
            parameter,weight,s,h_kickoff,nflag,wells,viscosity,Con,...
            nflagc,weightc,sc,dparameter,contnorm,SS,dt,h,MM,gravrate,...
            source,kmap,nflagface,N,theta_s,theta_r,alpha,pp,q,Kde,Ded,Kn,Kt,Hesq);
    elseif strcmp(acel,'AA')
        %% Picard-Anderson Acceleration
        [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
            parameter,weight,s,h_kickoff,nflag,wells,viscosity,0,contnorm,Con,nflagc,...
            weightc,sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate,source);
    end
else
    % Montagem da matriz global
    [M,I,elembedge] = ferncodes_globalmatrix(weight,s,Kde,Ded,Kn,Kt,Hesq,...
        viscosity,nflag,nflagface,SS,dt,h,MM,gravrate,kmap);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    %Often it may change the global matrix "M" with wells
    [M,I] = addsource(sparse(M),I,wells);

    % Often with source term
    [I]=sourceterm(I,source);

    p = solver(M,I);
    wc=0;Kdec=0;Knc=0;Ktc=0;Dedc=0;
    % auxiliary variables interpolation
    [pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflagno,weight,s,Con,nflagc,wc,sc);
    %Get the flow rate (Diamond)
    [flowrate,flowresult,flowratedif] = ferncodes_flowrate(p,pinterp,cinterp,Kde,...
        Ded,Kn,Kt,Hesq,viscosity,nflagno,Con,Kdec,Knc,Ktc,Dedc,nflagc,gravrate);
end





