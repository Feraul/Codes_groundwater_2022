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
%--------------------------------------------------------------------------

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
    recovtype keygravity g benchkey rowposit nltol maxiter acel bcflagc;

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
% utilize Tipo1Malha4.msh, Tipo1Malha5.msh, Tipo1Malha6.msh

% x=bedge(:,1);
% y=bedge(:,2);
% bedge(:,1)=y;
% bedge(:,2)=x;

% necesitamos ajeitar para o campo de pressao quando tem "restart"
if numcase >200
    % Flags adequation for contamination and groundwater problem
     [bedge,bcflagc,wellsc,auxpar,velmedio,dmap,Dmedio,gamma,kmap]=...
         preconcentration(bedge,wells,kmap);
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
    inedgesize,auxpar, wellsc,velmedio,dmap,Dmedio,gamma);

