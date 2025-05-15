%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in poorus media
%Modified: Fernando Contreras, 2021
function [p,flowrate,flowresult,flowratedif] = ferncodes_solverpressure(viscosity,...
    wells,Hesq,Kde,Kn,Kt,Ded,nflag,nflagface,weight,s,Con,Kdec,Knc,Ktc,Dedc,...
    nflagc,wc,sc,SS,dt,h,MM,gravrate,P,kmap,tempo,source)
% Montagem da matriz global
[M,I,elembedge] = ferncodes_globalmatrix(weight,s,Kde,Ded,Kn,Kt,Hesq,...
    viscosity,nflag,nflagface,SS,dt,h,MM,gravrate,kmap);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M" with wells
[M,I] = addsource(sparse(M),I,wells);

% Often with source term
[I]=sourceterm(I,source);

%--------------------------------------------------------------------------
%Solve global algebric system: pressure or hydraulic head 
p = solver(M,I);
%Message to user:
disp('>> The Pressure or hydraulic head field was calculated with success!');

% auxiliary variables interpolation 
[pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflag,weight,s,Con,nflagc,wc,sc);
%Get the flow rate (Diamond)
[flowrate,flowresult,flowratedif] = ferncodes_flowrate(p,pinterp,cinterp,Kde,...
    Ded,Kn,Kt,Hesq,viscosity,nflag,Con,Kdec,Knc,Ktc,Dedc,nflagc,gravrate);


%Message to user:
disp('>> The Flow Rate field was calculated with success!');
