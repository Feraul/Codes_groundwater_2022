%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Programer: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or
%heterogen domain such as isotropic and anisotropic media for each time
%step or in the steady state when will be important.

%--------------------------------------------------------------------------
%This routine receives geometry and physical data.

%--------------------------------------------------------------------------

function hydraulic(wells,overedgecoord,V,N,Hesq,Kde,Kn,Kt,Ded,kmap,nflag,...
    parameter,h_old,contnorm,SS,MM,weight,s,dt,gravrate,...
    nflagface,weightDMP,P)
%Define global parameters:
global timew  totaltime coord pmethod filepath elem numcase inedge ...
        bedge interptype centelem ;

%---------------------------------------------------------------------------
% h: representa a carga hidraulica 
% h_new: carga hidraulica atualizado a cada passo de tempo
% h_old: carga hidraulica inicial
%--------------------------------------------------------------------------
%Initialize parameters:
hanalitica=0;
hauxiliar=0;
%"time" is a parameter which add "dt" in each looping (dimentional or adm.)

time = 0;

stopcriteria = 0;
orderintimestep = ones(size(elem,1),1)*0;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
timew = 0;
% inicialization variables
satonvertices=0;producelem=0;h=h_old;Con=0;Kdec=0;Knc=0;nflagc=0;viscosity=1;
contiterplot=1;auxkmap=0;mobility=1;Ktc=0;Dedc=0;wightc=0;sc=0;dparameter=0;
% armazena os vtk s no tempo 0
postprocessor(h_old,zeros(size(inedge,1)+size(bedge,1),1),Con,1-Con,...
    contiterplot,overedgecoord,orderintimestep,'i',1,auxkmap);
if numcase==342
        %----------------------------------------------------------------------
        vx=63; % malha quadrilateral ortogonal e distorcida
        %vx=212;% malha triangular nao-estruturada
        hauxiliar(1,1)=time;
        hauxiliar(1,2)=h_old(vx);
        hanalitica(1,1)=time;
        hanalitica(1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*time/(3.28*10^-3))));
end
dtaux=dt;
tic
%Get "hydraulic head" and "flowrate"
while stopcriteria < 100
    % utiliza o metodo TPFA para aproximar a carga hidraulica
    if strcmp(pmethod,'tpfa')
        [h_new,flowrate,] = ferncodes_solvePressure_TPFA(Kde, Kn,...
            nflagface, Hesq,wells,viscosity, Kdec, Knc,nflagc,Con,SS,dt,h,MM,P,time);
        % utiliza o metodo MPFA-D para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfad')
        % Calculate hydraulic head and flowrate using the MPFA with diamond pacth
        [h_new,flowrate,] = ferncodes_solverpressure(...
            mobility,wells,Hesq,Kde,Kn,Kt,Ded,nflag,nflagface,...
            weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,wightc,sc,SS,dt,h,MM,...
            gravrate,P,kmap,time);
        % utiliza o metodo MPFA-H para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfah')
        % Calculate hydraulic head and flowrate using the MPFA with harmonic
        % points
        [h_new,flowrate,]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells,SS,dtaux,h,MM,gravrate,viscosity,P,time);
        % utiliza o metodo NL-TPFA para aproximar a carga hidraulica
    elseif strcmp(pmethod,'nlfvpp')
        [h_new,flowrate,]=...
            ferncodes_solverpressureNLFVPP(nflag,parameter,kmap,wells,...
            mobility,V,Con,N,p_old,contnorm,weight,s,Con,nflagc,wightc,...
            sc,dparameter,SS,dt,h,MM);
    end
    %% time step calculation
    disp('>> Time evolution:');
    time = time + dt
    % calcula a evoluicao da simulacao em porcentagem
    concluded = time*100/finaltime;
    stopcriteria = concluded;
    concluded = num2str(concluded);
    status = [concluded '% concluded']
    % contador
    contiterplot=contiterplot+1
    % atualiza a carga hidraulica 
    h=h_new;
    %----------------------------------------------------------------------
    % case unconfined aquifer
    % case 1 and 4 of the article Qian, et al 2023
    if numcase==333 || numcase==331
        if strcmp(pmethod,'mpfah')
            [facelement]=ferncodes_elementfacempfaH;
            [~,kmap] = calcnormk(kmap,MM,h);
            [pointarmonic]=ferncodes_harmonicopoint(kmap);
            [parameter,]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
            [weightDMP]=ferncodes_weightnlfvDMP(kmap,elem);
        elseif strcmp(pmethod,'mpfad')
            [~,kmap] = calcnormk(kmap,MM,h);
            %Get preprocessed terms:
            [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem);
            %It switches according to "interptype"
            switch char(interptype)
                %LPEW 1
                case 'lpew1'
                    % calculo dos pesos que correspondem ao LPEW1
                    [weight,s] = ferncodes_Pre_LPEW_1(kmap,N);
                    %LPEW 2
                case 'lpew2'
                    % calculo dos pesos que correspondem ao LPEW2
                    [weight,s] = ferncodes_Pre_LPEW_2(kmap,N,zeros(size(elem,1),1));
            end  %End of SWITCH
        else % TPFA
           [~,kmap] = calcnormk(kmap,MM,h);
            %Get preprocessed terms:
            [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem); 
        end
    elseif numcase==342 && strcmp(pmethod,'mpfah')
        %----------------------------------------------------------------------
        vx=63; % malha quadrilateral ortogonal distorcida
        %vx=212; % malha triangular nao estruturada
           hauxiliar(contiterplot+1,1)=log10(time);
           hauxiliar(contiterplot+1,2)=h_new(vx);
           hanalitica(contiterplot+1,1)=log10(time);
           hanalitica(contiterplot+1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*(time)/(3.28*10^-3)))); 
        % a cada passo de tempo atualiza a condicao de contorno do outro da
        % face a direita
        vv=find(bedge(:,5)==102);
        gg=(coord(bedge(vv,1),:) +coord(bedge(vv,2),:))*0.5;
        %hbound=h(gg);
        hbound=3*erfc(gg(:,1)/(2*sqrt(30.5*(time)/(3.28*10^-3))));
        nflagface(vv,2)=hbound;
        
    end
    dtaux=dt;
    % armanzena os vtks e calcula alguns erros  
    postprocessor(h,flowrate,Con,1-Con,contiterplot,overedgecoord,orderintimestep,'i',1,auxkmap);
    
end

toc
% plotagem dos graficos em determinados regioes do dominio
plotandwrite(producelem,Con,h,satonvertices,0,0,0,0,overedgecoord,hanalitica,hauxiliar);

%--------------------------------------------------------------------------
% activate if you want to visualize the hydraulic field in the final time

%h

%--------------------------------------------------------------------------
% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Hydraulic head extrema values [hmax hmin]:');
max_hyval = max(h)
min_hyval = min(h)
if numcase==342
erro=norm(hauxiliar(:,2)-hanalitica(:,2))
end
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
