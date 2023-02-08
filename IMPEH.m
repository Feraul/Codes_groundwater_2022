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
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or
%heterogen domain such as isotropic and anisotropic media for each time
%step or in the steady state when will be important.

%--------------------------------------------------------------------------
%This routine receives geometry and physical data.

%--------------------------------------------------------------------------

function IMPEH(Sw,injecelem,producelem,satinbound,wells,klb,satonvertices,...
    satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,lsw,...
    transmvecleft,transmvecright,knownvecleft,knownvecright,mapinv,...
    maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,Fg,overedgecoord,...
    bodyterm,normk,limiterflag,massweigmap,othervertexmap,V,N,Hesq,Kde,Kn,...
    Kt,Ded,kmap,nflag,swsequence,ntriang,areatriang,lastimelevel,...
    lastimeval,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound,elemsize,bedgesize,inedgesize,parameter,...
    weightDMP,nflagface,h_old,contnorm,SS,B)
%Define global parameters:
global timew elemarea totaltime timelevel pormap numcase pmethod smethod ...
    filepath benchkey resfolder order inedge bedge;

%--------------------------------------------------------------------------
%Initialize parameters:
c = 0;
%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
time = lastimeval;
stopcriteria = 0;
%"timelevel" is a parameter used to plot each time step result. Image
%1,2,...,n. This parameter is incremented and sent to "postprocessor"
timelevel = lastimelevel + 1;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
%"earlyswonedge" is a key used for define the strategy to calculate the
%mobility. See the "getmobility" function to more detail.
earlysw = 0;
timew = 0;
%Initialize "flagtoplot".  Initialy it is "0" and when is "1", plot the
%vtk. It avoids too much vtk files.
flagtoplot = 0;
%"contiterplot" is the number of the vtk created.
contiterplot = 1;

%Parameters to plot
%It is an auxiliary counter. When it get "10", the production parameters
%and saturation field are stored in a file *.dat
countstore = 0;
oilaccum = 0;
S_extrema_old = zeros(1,2);

Sleft = 0;
Sright = 0;

h=h_old;
tic
while stopcriteria < 100
    
    % calculo das vazões
     %Get "hydraulic charge" and "flowrate"
    [h,flowrate,flowresult] = solvePressure_TPFA(Kde, Kn,nflag, Hesq,wells);
    %% caculo do passo de tempo
    %dt=steptime(flowrate,CFL,auxflag,order);
    dt=0.0001;
    
    %% aproximação upwind
    RHS=zeros(size(elem,1),1);
    for iface = 1 : size(inedge,1)
        
        lef = inedge(iface,3); % elemento a esquerda
        rel = inedge(iface,4); % elemento a direita
        ve_mais  = (flowrate(iface+size(bedge,1)) + abs(flowrate(iface+size(bedge,1))))/2;
        ve_menos = (flowrate(iface+size(bedge,1)) - abs(flowrate(iface+size(bedge,1))))/2;
        
        RHS(rel)  = RHS(rel)  + ve_mais  + ve_menos;
        RHS(lef)  = RHS(lef)  - ve_mais  - ve_menos;
        
    end
    %% contribuição do contorno
    for ifacont = 1:size(bedge,1)
        
        RHS(bedge(ifacont,3)) = RHS(bedge(ifacont,3)) - flowrate(ifacont);
        
    end
        t_old=t_old+dt
    %% calculo da concentração
    for i = 1:size(elem,1)
        
        h_new(i,1) = h(i,1) + dt*BB*((RHS(i))/(elemarea(i)*SS));
    end
    cont=cont+1
    h=h_new;
    postprocessor(h,solanal,cont);
end

toc
%Write data file ("ProdutionReport.dat" and others)

plotandwrite(producelem,Sw,pressure,overedgecoord(:,1),injecelem);

%--------------------------------------------------------------------------
%Write data file ("ProdutionReport.dat" and others)

% plotandwrite(producelem,Sw,pressure,overedgecoord(:,1),injecelem);

% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Saturation extrema values [Smax Smin]:');
max_satval = max(maxminsatval(:,1))
min_satval = min(maxminsatval(:,2))
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
