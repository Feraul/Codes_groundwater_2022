%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media
%Type of file: FUNCTION
%Criate date: 16/04/2015 (My wife is a PHD since yesterday)
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure,flowrate,flowresult,flowratedif] = ferncodes_solvePressure_TPFA(Kde, Kn, nflag, Hesq,wells,viscosity, Kdec, Knc,nflagc,Con )
%Define global parameters:

global inedge bedge elem coord centelem bcflag numcase

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
mvector=zeros(size(elem,1),1);
bedgesize = size(bedge,1);
% Loop de faces de contorno

for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(ifacont,1);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        
        c1=nflag(bedge(ifacont,1),2);
        
        
        %Preenchimento
        
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- visonface*A*(norm(v0)^2);
        
        mvector(bedge(ifacont,3))=mvector(bedge(ifacont,3))-visonface*c1*A*(norm(v0)^2);
        
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        mvector(bedge(ifacont,3))=mvector(bedge(ifacont,3))- norm(v0)*bcflag(r,2);
        
        
    end
end

for iface=1:size(inedge,1),
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(bedgesize + iface,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-visonface*Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+visonface*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-visonface*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+visonface*Kde(iface,1);
end




%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M,mvector] = addsource(M,mvector,wells);

%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver".
%This function returns the pressure field with value put in each colocation
%point.
[pressure] = solver(M,mvector);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%--------------------------------------------------------------------------
%Once the pressure was calculated, the "flowrate" field is also calculated

%Calculate flow rate through edge. "satkey" equal to "1" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[flowrate, flowresult,flowratedif]=ferncodes_flowrateTPFA(pressure,Kde,Kn,Hesq,nflag,viscosity,Kdec, Knc,nflagc,Con);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
