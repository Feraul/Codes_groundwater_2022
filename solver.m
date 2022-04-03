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
%Goals: %This function solves the global algebric system.
%"globalM" is a global matrix

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure] = solver(globalmatrix,vector)
%Solver the algebric system using matlab's routines (black box!)
pressure = globalmatrix\vector;
%[L,U] = ilu(sparse(globalmatrix),struct('type','ilutp','droptol',1e-6));
       
  %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
 % [pressure,fl1,rr1,it1,rv1]=gmres(globalmatrix,vector,10,1e-9,10000,L,U);