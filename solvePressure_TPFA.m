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
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure,flowrate,flowresult,flowratedif] = solvePressure_TPFA(transmvecleft,...
    knownvecleft,viscosity,wells,Fg,bodyterm,Con,transmvecleftc,aa,SS,dt,h,MM)
%Define global parameters:
global elem bedge inedge phasekey numcase elemarea;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize the global matrix which have order equal to parameter
%"size(elem)".
M = zeros(size(elem,1));
%Initialize "mvector" which is the independent vector of algebric system.
mvector = zeros(size(elem,1),1);
%Swept "bedge"
for ibedg = 1:bedgesize
    %Get "leftelem"
    leftelem = bedge(ibedg,3);
    %viscosidade
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(ibedg,1);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    
    %Fill the global matrix "M" and known vector "mvector"
    M(leftelem,leftelem) = M(leftelem,leftelem) + ...
        visonface*transmvecleft(ibedg);
    %Update "transmvecleft"
    transmvecleft(ibedg) = visonface*transmvecleft(ibedg);
     if numcase==330
        % Letf
        transmvecleft(ibedg) = transmvecleft(ibedg);
    end
    %Fill "mvector"
    mvector(leftelem) = ...
        mvector(leftelem) + visonface*knownvecleft(ibedg);
end  %End of FOR ("bedge")

%Swept "inedge"
for iinedg = 1:inedgesize
    %Get "leftelem" and "right"
    leftelem = inedge(iinedg,3);
    rightelem = inedge(iinedg,4);
    %viscosidade
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(bedgesize + iinedg,1);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    
    %Fill the global matrix "M" and known vector "mvector"
    %Contribution from the element on the LEFT to:
    %Left
    M(leftelem,leftelem) = M(leftelem,leftelem) + ...
        visonface*transmvecleft(bedgesize + iinedg);
    %Right
    M(leftelem,rightelem) = M(leftelem,rightelem) - ...
        visonface*transmvecleft(bedgesize + iinedg);
    %Contribution from the element on the RIGHT to:
    %Right
    M(rightelem,rightelem) = M(rightelem,rightelem) + ...
        visonface*transmvecleft(bedgesize + iinedg);
    %Right
    M(rightelem,leftelem) = M(rightelem,leftelem) - ...
        visonface*transmvecleft(bedgesize + iinedg);
    if numcase==330
        % Letf
        M(leftelem,leftelem) = M(leftelem,leftelem) + ...
            coeficiente;
        %Right
        M(rightelem,rightelem) = M(rightelem,rightelem) + ...
            coeficiente;
    end
    %Update "transmvecleft"
    transmvecleft(bedgesize + iinedg) = ...
        visonface*transmvecleft(bedgesize + iinedg);
end  %End of FOR ("inedge")
%==========================================================================
% para calcular a carga hidraulica
% para calcular a carga hidraulica
if numcase>300
    if numcase~=336 && numcase~=334 && numcase~=340
        if numcase==333 || numcase==331
            coeficiente=dt^-1*SS.*elemarea(:);
        else
            coeficiente=dt^-1*MM*SS.*elemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            M=M+coeficiente.*eye(size(elem,1));
            transmvecleft=transmvecleft+coeficiente.*eye(size(elem,1))*h;
        else
            % Crank-Nicolson method
            transmvecleft=transmvecleft+coeficiente.*eye(size(elem,1))*h-0.5*M*h;
            M=0.5*M+coeficiente.*eye(size(elem,1));
            
        end
        
    end
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
[flowrate,flowresult,flowratedif] = calcflowrateTPFA(transmvecleft,pressure,Con,transmvecleftc,aa);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
