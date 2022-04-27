
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 10/01/2012
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

function IMPEC(Con,injecelem,producelem,satinbound,wells,klb,satonvertices,...
    satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,lsw,...
    transmvecleft,transmvecright,knownvecleft,knownvecright,mapinv,...
    maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,Fg,overedgecoord,...
    bodyterm,normk,limiterflag,massweigmap,othervertexmap,V,N,Hesq,Kde,Kn,...
    Kt,Ded,kmap,nflag,swsequence,ntriang,areatriang,lastimelevel,...
    lastimeval,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound,elemsize,bedgesize,inedgesize,parameter,...
    weightDMP,nflagface,p_old,contnorm,dmap,dparameter,nflagc,gamma,...
    Dmedio,Kdec,Knc,Ktc,Dedc,wightc,sc,nflagfacec,weightDMPc,wellsc)
%Define global parameters:
global timew elemarea totaltime timelevel pormap numcase pmethod smethod ...
    filepath benchkey resfolder order bcflagc bedge inedge elem interptype ;

%--------------------------------------------------------------------------
%Initialize parameters:
c = 0;
Conaux=Con;
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
timew = 0;
%Initialize "flagtoplot".  Initialy it is "0" and when is "1", plot the
%vtk. It avoids too much vtk files.
flagtoplot = 0;
%"contiterplot" is the number of the vtk created.
contiterplot = 1;
earlysw = 0;
oilaccum = 0;

Sleft = 0;
Sright = 0;
%Parameters to plot
%It is an auxiliary counter. When it get "10", the production parameters
%and saturation field are stored in a file *.dat

C_extrema_old = zeros(1,2);
%--------------------------------------------------------------------------
%Verify if there exists a restart

%Verify if the "lastimelevel" is bigger than zero
%In this case, the production parameters must be cought

%Call the postprocessor:
if numcase==247
    postprocessor(ones(elemsize,1),0,Con,contiterplot - 1,...
        order*ones(elemsize,1),'i',1,log(kmap(:,2)));
else
    postprocessor(ones(elemsize,1),0,Con,contiterplot - 1,...
        order*ones(elemsize,1),'i',1,normk);
end

if numcase~=246 && numcase~=246 && numcase~=247 && numcase~=248 && numcase~=249 && numcase~=250
    if strcmp(pmethod,'tpfa') && numcase~=31.1
        
        %Get "pressure" and "flowrate"
        [pressure,flowrateadvec,flowresult] = solvePressure_TPFA(Kde, Kn, nflag, Hesq,wells);
      
        %MPFA-D (Gao and Wu, 2010)
    elseif strcmp(pmethod,'mpfad') &&  numcase~=31.1
        
        %Calculate "pressure", "flowrate" and "flowresult"
        [pressure,flowrateadvec,flowresult] = ferncodes_solverpressure(kmap,...
            1,wells,Con,V,N,Hesq,Kde,Kn,Kt,Ded,nflag);
       
    elseif strcmp(pmethod,'mpfaql')
        [pressure,flowrateadvec,flowresult]=ferncodes_solverpressureMPFAQL(nflag,...
            parameter,kmap,weightDMP,wells,1,V,Con,N);
       
    elseif strcmp(pmethod,'mpfah')
        
        [pressure,flowrateadvec,flowresult]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells);
       
    elseif strcmp(pmethod,'nlfvpp')
        
        [pressure,flowrateadvec,flowresult]=ferncodes_solverpressureNLFVPP(nflag,...
            parameter,kmap,wells,1,V,Con,N,p_old,contnorm);
        
        %Any other type of scheme to solve the Pressure Equation
    elseif strcmp(pmethod,'nlfvh')&& numcase~=31.1 % revisar com cuidado
        
        [pressure,flowrate,flowresult]=ferncodes_solverpressureNLFVH(nflagface,...
            parameter,wells,1,weightDMP,p_old,0,0,0,contnorm);
    elseif strcmp(pmethod,'nlfvdmp')&& numcase~=31.1
        
        [pressure,flowrate,flowresult]=ferncodes_solverpressureDMP(nflagface,...
            parameter,wells,1,weightDMP,p_old,0,0,0,contnorm);
    elseif numcase~=31.1
        
        %Calculate the PRESSURE field (Two-Phase context):
        [pressure,flowrate,flowresult] = solvePressure(transmvecleft,...
            transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
            Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
            bodyterm);
    end  %End of IF (type of pressure solver)
end
%It switches according to "interptype"
switch char(interptype)
    %LPEW 1
    case 'lpew1'
        % calculo dos pesos que correspondem ao LPEW1
        [weight,s] = ferncodes_Pre_LPEW_1(kmap,1,V,Sw,N);
        %LPEW 2
    case 'lpew2'
        % calculo dos pesos que correspondem ao LPEW2
        [weight,s] = ferncodes_Pre_LPEW_2(kmap,N);
        
end  %End of SWITCH

if numcase==233
    velmedio=0.5; % quando a pessão no contorno é 5
elseif numcase==235
    velmedio=2;    % quando a pressão no contorno é 15
elseif numcase==231 || numcase==232 || numcase==243 || ...
        numcase==236 || numcase==237 || numcase==238 ||...
        numcase==239 || numcase==241 || numcase==242 ||...
        numcase==244  || numcase==234|| numcase==233 || numcase==248
    velmedio=1;
end

q=zeros(size(elem,1),1);
while stopcriteria < 100
    
    % while time< finaltime
    %User message:
    %Jump a row (in prompt)
    disp(' ');
    disp('---------------------------------------------------');
    disp('>> Show timelevel:')
    timelevel
    
    if numcase==246 || numcase==245 || numcase==247 || numcase==248 || numcase==249 || numcase==250
        [viscosity] = ferncodes_getviscosity(satinbound,injecelem,Con,earlysw,smethod,...
            timelevel,numcase,Sleft,Sright,c,overedgecoord,nflagc,nflagfacec);
        % calculo da pressão e fluxo
    if strcmp(pmethod,'nlfvpp')    
        [pressure,flowrateadvec,flowresult,flowratedif]=ferncodes_solverpressureNLFVPP(nflag,...
            parameter,kmap,wells,viscosity,V,Con,N,p_old,contnorm,weight,s,Con,nflagc,wightc,sc,dparameter);
    elseif strcmp(pmethod,'mpfad')
       [pressure,flowrateadvec,flowresult,flowratedif] = ferncodes_solverpressure(...
            viscosity,wells,Hesq,Kde,Kn,Kt,Ded,nflag,weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,wightc,sc); 
    else
       [pressure,flowrateadvec,flowresult] = solvePressure_TPFA(transmvecleft,...
            knownvecleft,viscosity,wells,Fg,bodyterm,Con);  
    end
       
    end
    
    %----------------------------------------------------------------------
    %Calculate "dt" using the function "calctimestep"
    
    %This function obtains the time step using the "Courant" number.
    %The necessity of calculate the time step is ensure the stability of
    %explicit saturation formulation.
    
    dt = calctimestep(flowrateadvec,satinbound,gamma,Dmedio)
    
    %----------------------------------------------------------------------
    
    %Calculate the CONCENTRATION field (choose saturation method):
    [newC,cflowrate,Sleft,Sright] = ...
        solveSaturation(Con,flowrateadvec,flowratedif,dt,injecelem,producelem,satinbound,...
        Fg,flagknownvert,satonvertices,flagknownedge,satonedges,flowresult,...
        wvector,wmap,constraint,lsw,limiterflag,1,massweigmap,...
        othervertexmap,swsequence,ntriang,areatriang,prodwellbedg,...
        prodwellinedg,mwmaprodelem,vtxmaprodelem,coordmaprodelem,...
        amountofneigvec,rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
        isonbound,elemsize,bedgesize,inedgesize,gamma,time);
    
    %Update the saturation field
    Con = newC;
    
    %----------------------------------------------------------------------
    %Define PVI or DIMENTIONAL TIME
    
    %Dimentional (s, h, day, etc)
    
    time = time + dt
    concluded = time*100/finaltime;
    %It is used for restart activation
    percentdt = dt*100/finaltime;
    stopcriteria = concluded;
    %Define a flag to plot the vtk file
    flagtoplot = flagtoplot + dt*100/finaltime;
    concluded = num2str(concluded);
    status = [concluded '% concluded']
    
    
    if any(producelem)
                
        %Write table (Time (VPI), Oil Flow rate, Accumulated Oil and Water Cut)
        %Create the file name
        prfilename = [resfolder '_' 'ProdutionReport.dat'];
        
        %Select the "letter" for "fopen"
        if timelevel == 1
            letter = 'w';
        else
            letter = 'a';
        end  %End of IF
        
        %Open the file
        writeproductionreport = fopen([filepath '\' prfilename],letter);
        
        %Write "productionreport" according to "producelem" size
        %There is one producer well
        if length(producelem) == 1
            fprintf(writeproductionreport,'%26.16E %26.16E\r\n',[time cflowrate]);
            %There are more than one producer wells (but just write two)
        else
            fprintf(writeproductionreport,'%26.16E %26.16E\r\n',[time cflowrate]);
        end  %End of IF
        
        %Close the file "writeproductionreport.dat"
        fclose(writeproductionreport);
    end  %End of IF (when there exists producer well)
    
    %----------------------------------------------------------------------
    %Call the "postprocessor" (plot results in each time step)
    
    %Just create the vtk file if "flagtoplot" reaches 0.1.
    if flagtoplot >= 1
        %This function create the "*.vtk" file used in VISIT to posprocessing
        %the results
        if numcase==247
            postprocessor(pressure,flowrateadvec,Con,contiterplot,overedgecoord,'i',1,log(kmap(:,2)));
        else
            postprocessor(pressure,flowrateadvec,Con,contiterplot,overedgecoord,'i',1,normk);
            %Update "flagtoplot"
        end
        flagtoplot = 0;
        %Update "contiterplot"
        contiterplot = contiterplot + 1;
    end  %End of IF
    
    %User mesage
    disp('>> Concentration field calculated with success!');
    disp('>> Concentration extrema values [Con_max con_min]:');
    %Show extrema values
    C_extrema = [max(Con); min(Con)]
    
    %Store maximum and minimum saturation values:
    maxminconval = [max(C_extrema(1),C_extrema_old(1)) min(C_extrema(2),C_extrema_old(2))];
    %Update "S_extrema"
    C_extrema_old = C_extrema;
    
    %Increment the parameters "timelevel" and "countstore"
    timelevel = timelevel + 1;
    
    c=c+1;
    %It gives the time spent per "timelevel"
    if time>1000 && (numcase==242 || numcase==243 || numcase==249 || numcase==250)
        %if time>0.5 && (numcase==242 || numcase==243 || numcase==245)
        if numcase==245 || numcase==247
            wellsc(1:2,4)=0;
            Con(wells(1,1),1)=wells(1,4);
            Con(wells(2,1),1)=wells(2,4);
            %Initialize and preprocess the parameters:
            [nflagc,nflagfacec] = ferncodes_calflag_con(0);
            %Define elements associated to INJECTOR and PRODUCER wells.
            [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                Conaux,klb);
            
            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,1);
            
        else
            bcflagc(2,2)=0; % valor dirich
            %Initialize and preprocess the parameters:
            [nflagc,nflagfacec] = ferncodes_calflag_con(0);
            %Define elements associated to INJECTOR and PRODUCER wells.
            [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                Conaux,klb);
            
            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,0);
            
        end
    end
    
    if numcase==248
        %Initialize and preprocess the parameters:
        nflag = ferncodes_calflag(time);
        [nflagc,nflagfacec] = ferncodes_calflag_con(time);
        %Define flags and known saturation on the vertices and edges.
        [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
            getsatandflag(satinbound,injecelem,Con,nflagc,nflagfacec);
    end
    %[analsol]=anasolaux(velmedio,Dmedio,time);
end  %End of While
toc
%Write data file ("ProdutionReport.dat" and others)

plotandwrite(producelem,Con,pressure,satonvertices,Dmedio,velmedio,gamma);


%--------------------------------------------------------------------------

% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Concentration extrema values [Con_max Con_min]:');
max_conval = max(maxminconval(:,1))
min_conval = min(maxminconval(:,2))
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system

system(command);
