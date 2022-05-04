
%Modified: Fernando Contreras, 2021

function [p,flowrate,flowresult]=ferncodes_solverpressureMPFAQL(nflagno,parameter,kmap,weightDMP,wells,mobility,V,Sw,N)
%Define global parameters
global interptype;

%It switches according to "interptype"
switch char(interptype)
    %LPEW 1
    case 'lpew1'
        % calculo dos pesos que correspondem ao LPEW1
        [wight,s] = ferncodes_Pre_LPEW_1(kmap,mobility,V,Sw,N);
    %LPEW 2
    case 'lpew2'
        % calculo dos pesos que correspondem ao LPEW2
        [wight,s] = ferncodes_Pre_LPEW_2(kmap,N);
end  %End of SWITCH

[M,I]=ferncodes_assemblematrixMPFAQL(parameter,wight,s,nflagno,weightDMP,mobility);
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
[pinterp]=ferncodes_pressureinterpMPFAQL(p,nflagno,wight,s);
%Get the flow rate (Diamond)
[flowrate,flowresult]=ferncodes_flowratelfvMPFAQL(parameter,weightDMP,mobility,pinterp,p);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end