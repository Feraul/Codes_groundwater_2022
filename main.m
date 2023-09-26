%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: MAIN
%Criate date: 29/02/2012
%Modify data:   /  /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals: Do the manegement of simulator. This is a MAIN program. 

%--------------------------------------------------------------------------
%In this numerical routine the flow may be simulated with one or two phase 
% or contaminants or groundwater. The functions below are organized in order
%give flexibility to software resourcers.
%For example: the saturation and pressure fields are calculated by IMPES,
%but it may be calculated also by a fully implicit scheme. This change is
%done in the present rountine just call each function.

%|--------------------------------------|
%|  leia com muita atencao as indicacoes|
%|--------------------------------------|

%Clear the screem 
clc;
%Clear all the memory of the matlab
clear all;
%Define the format of out data
format short;
%It begins the time counter and "profile".
tic
% profile on

%--------------------------------------------------------------------------
%Define the global variables:
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath resfolder ...
    numcase pmethod smethod phasekey order timeorder auxcvfactor ...
    interptype multdopt goefreeopt lsneightype lsexp ...
    recovtype keygravity g benchkey rowposit nltol maxiter acel bcflagc...
    methodhydro SS MM P;

%--------------------------------------------------------------------------
%Call the "preprocessor" function

%This function obtains from gmsh file all data structure.
[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,pmethod,...
    smethod,xyrz,r0,symaxe,keymsfv,coarseratio,auxcvfactor,interptype,...
   multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,filepath,resfolder,benchkey,...
    kmap,wells,klb,limiterflag,rowposit,nltol,maxiter,acel] = preprocessor;
% As concicoes abaixo utilize para malhas Tipo1Malha4.msh, Tipo1Malha5.msh, Tipo1Malha6.msh
% x=bedge(:,1);
% y=bedge(:,2);
% bedge(:,1)=y;
% bedge(:,2)=x;
% necesitamos ajeitar para o campo de pressao quando tem "restart"
% this numcase is used to simulate concentration solute in aquifers
if 200<numcase && numcase <300
    % Flags adequation for contamination and groundwater problem
    [bedge,bcflagc,wellsc,auxpar,velmedio,dmap,Dmedio,gamma]=...
        preconcentration(bedge,wells);
    % calculate permeability tensor
    if numcase==247 || numcase==249 || numcase==250
        % permeability field obtained by: Nicolaides, Cueto-Filgueroso,
        % Juanes 2015.
        % reorganização da matriz de permeabilidade
        load('Perm_Var0p1.mat')
        %==================================================================
        perm=flipud(perm);
        auxperm2=perm(1:125,1:125);
        m=1;
        for i=1:125
            for j=1:125
                auxperm3(m,1)=auxperm2(j,i);
                m=m+1;
            end
        end
        kmap=auxperm3;
        %adeSPE; % para um campo de permeabilidade da SPE active descomente.
        %=================================================================
    elseif numcase==251
        kmap=kmap;
        %adeSPE; % para um campo de permeabilidade da SPE active descomente.
    end
    SS=0; h_old=0; MM=0;wells=0; dt=0;P=0; 
    % This numcase is used to simulate head hydraulic in aquifers
elseif numcase>300
    
    auxpar=0;
    dmap=0;
    Dmedio=0;
    gamma=0;
    velmedio=0;
    wellsc=0;
   % Flags adequation for hydrological head problem
   [SS,h_old,MM,wells,dt,P]=prehydraulic;
   % Choose Backward method or Crank-Nicolson
   % use oara caso aquifero confinado
   methodhydro='backward';
   % use para o caso aquifero nao confinado
   %methodhydro='cranknicolson'; 
   % This numcase is used for two-phase (water-oil) flow problems  
else
    auxpar=0;
    dmap=0;
    Dmedio=0;
    gamma=0;
    velmedio=0;
    wellsc=0;
    SS=0; h_old=0; MM=0; dt=0; P=0; 
end
%adeSPE; % para um campo de permeabilidade da SPE active descomente.
%--------------------------------------------------------------------------
%Call "setmethod"
%elem(:,5)=1;
%Initialize general parameters:
elemsize = size(elem,1);
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%It preprocess the schemes and set a One-phase or Two-phase simulation.
setmethod(kmap,wells,'i',8,limiterflag,klb,elemsize,bedgesize,...
    inedgesize,auxpar, wellsc,velmedio,dmap,Dmedio,gamma,SS,h_old,MM,dt,P);

