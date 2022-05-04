%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 18/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras,
function [p,flowrate,flowresult,flowratedif]=ferncodes_solverpressureNLFVPP(nflag,...
                                  parameter,kmap,wells,viscosity,V,Sw,N,...
                                  p_old,contnorm,wight,s,Con,nflagc,wightc,sc,dparameter)
%Define global parameters
global acel;

% interpolação nos nós ou faces
[pinterp,]=ferncodes_pressureinterpNLFVPP(p_old,nflag,wight,s,Con,nflagc,wightc,sc);
[M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,contnorm);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);
%% full Picard iteration
if strcmp(acel,'FPI')
    [p,flowrate,flowresult,flowratedif]=ferncodes_iterpicard(M_old,RHS_old,...
        parameter,wight,s,p_old,nflag,wells,viscosity,Con,nflagc,wightc,sc,dparameter);
elseif strcmp(acel,'AA')
    %% Picard-Anderson Acceleration
    [p,flowrate,flowresult,flowratedif]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
        parameter,wight,s,p_old,nflag,wells,viscosity,0,contnorm,Con,nflagc,wightc,sc,dparameter);
end
end