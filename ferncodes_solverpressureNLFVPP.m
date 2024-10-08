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
%Programer: M�rcio Souza
%Modified: Fernando Contreras,
function [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_solverpressureNLFVPP(nflag,...
                                  parameter,kmap,wells,viscosity,V,N,...
                                  p_old,contnorm,weight,s,Con,nflagc,...
                                  weightc,sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate)
                            
%Define global parameters
global acel;

% interpola��o das pressoes e concentracoes nos vetrtices ou faces
[pinterp,]=ferncodes_pressureinterpNLFVPP(p_old,nflag,weight,s,Con,nflagc,weightc,sc);

% montagem das matriz global 
[M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,contnorm,SS,dt,h,MM,gravrate,p_old);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);
%% full Picard iteration
if strcmp(acel,'FPI')
    [p,flowrate,flowresult,flowratedif]=ferncodes_iterpicard(M_old,RHS_old,...
        parameter,weight,s,p_old,nflag,wells,viscosity,Con,nflagc,weightc,sc,...
        dparameter,contnorm,SS,dt,h,MM,gravrate);
elseif strcmp(acel,'AA')
    %% Picard-Anderson Acceleration
    [p,flowrate,flowresult,flowratedif,flowresultc]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
        parameter,weight,s,p_old,nflag,wells,viscosity,0,contnorm,Con,nflagc,...
        weightc,sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate);
end
end