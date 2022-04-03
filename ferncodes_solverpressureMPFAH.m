function [p,flowrate,flowresult]=ferncodes_solverpressureMPFAH(nflagface,parameter,weightDMP,wells)


[M,I]=ferncodes_assemblematrixMPFAH(parameter,nflagface,weightDMP);
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
[pinterp]=ferncodes_pressureinterpHP(p,nflagface,parameter,weightDMP);
%Get the flow rate 
[flowrate,flowresult]=ferncodes_flowratelfvHP(parameter,weightDMP,pinterp,p);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
end