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
    wells,Hesq,Kde,Kn,Kt,Ded,nflag,weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,wc,sc)

% Montagem da matriz global
[M,I] = ferncodes_globalmatrix(weight,s,Kde,Ded,Kn,Kt,Hesq,viscosity,nflag);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p = solver(M,I);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%Get the flow rate (Diamond)
[pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflag,weight,s,Con,nflagc,wc,sc);
[flowrate,flowresult,flowratedif] = ferncodes_flowrate(p,pinterp,cinterp,Kde,...
    Ded,Kn,Kt,Hesq,viscosity,nflag,Con,Kdec,Knc,Ktc,Dedc,nflagc);


%Message to user:
disp('>> The Flow Rate field was calculated with success!');
