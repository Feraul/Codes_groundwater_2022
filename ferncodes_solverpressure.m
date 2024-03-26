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

%Modified: Fernando Contreras, 2021
function [p,flowrate,flowresult,flowratedif] = ferncodes_solverpressure(viscosity,...
    wells,Hesq,Kde,Kn,Kt,Ded,nflag,nflagface,weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,...
    wc,sc,SS,dt,h,MM,gravrate,P,kmap)


% Montagem da matriz global
[M,I,elembedge] = ferncodes_globalmatrix(weight,s,Kde,Ded,Kn,Kt,Hesq,viscosity,nflag,nflagface,...
    SS,dt,h,MM,gravrate,kmap);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells,P,elembedge);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p = solver(M,I);
%Message to user:
disp('>> The Pressure field was calculated with success!');

% auxiliary variables interpolation 
[pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflag,weight,s,Con,nflagc,wc,sc);
%Get the flow rate (Diamond)
[flowrate,flowresult,flowratedif] = ferncodes_flowrate(p,pinterp,cinterp,Kde,...
    Ded,Kn,Kt,Hesq,viscosity,nflag,Con,Kdec,Knc,Ktc,Dedc,nflagc,gravrate);


%Message to user:
disp('>> The Flow Rate field was calculated with success!');
