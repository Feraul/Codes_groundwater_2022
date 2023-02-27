%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 10/01/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or
%heterogen domain such as isotropic and anisotropic media for each time
%step or in the steady state when will be important.

%--------------------------------------------------------------------------
%This routine receives geometry and physical data.

%--------------------------------------------------------------------------

function hydraulic(wells,overedgecoord,V,N,Hesq,Kde,Kn,Kt,Ded,kmap,nflag,...
               parameter,h_old,contnorm,SS,MM,weight,s,dt)
%Define global parameters:
global timew  totaltime  pmethod filepath ;

%--------------------------------------------------------------------------
%Initialize parameters:

%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
time = 0;
stopcriteria = 0;

%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
timew = 0;

satonvertices=0;
producelem=0;
h=h_old;

Con=0;
Kdec=0;
Knc=0;
nflagc=0;
viscosity=1;
contiterplot=0;
auxkmap=0;
mobility=1;
Ktc=0;Dedc=0;wightc=0;sc=0;dparameter=0;
tic
while stopcriteria < 100
    if strcmp(pmethod,'tpfa')
    % calculo das vaz�es
     %Get "hydraulic charge" and "flowrate"
    [h_new,flowrate,] = ferncodes_solvePressure_TPFA(Kde, Kn,...
        nflag, Hesq,wells,viscosity, Kdec, Knc,nflagc,Con,SS,dt,h,MM);
    elseif strcmp(pmethod,'mpfad')
    %Calculate "pressure", "flowrate" and "flowresult"
        [h_new,flowrate,] = ferncodes_solverpressure(...
            mobility,wells,Hesq,Kde,Kn,Kt,Ded,nflag,...
            weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,wightc,sc,SS,dt,h,MM);
    elseif strcmp(pmethod,'nlfvpp')
         [h_new,flowrate,]=...
            ferncodes_solverpressureNLFVPP(nflag,parameter,kmap,wells,...
            mobility,V,Con,N,p_old,contnorm,weight,s,Con,nflagc,wightc,...
            sc,dparameter,SS,dt,h,MM);
    end
        %% caculo do passo de tempo
    
    time = time + dt;
    
    concluded = time*100/finaltime;
    stopcriteria = concluded;
    contiterplot=contiterplot+1
    h=h_new;
    postprocessor(h,flowrate,Con,contiterplot,overedgecoord,'i',1,auxkmap);
end

toc
%Write data file ("ProdutionReport.dat" and others)

plotandwrite(producelem,Con,h,satonvertices,0,0,0,0);

%--------------------------------------------------------------------------
%Write data file ("ProdutionReport.dat" and others)

% plotandwrite(producelem,Sw,pressure,overedgecoord(:,1),injecelem);

% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Hydraulic head extrema values [hmax hmin]:');
max_satval = max(h)
min_satval = min(h)
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
