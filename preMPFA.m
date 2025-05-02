
%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in porous media
%Type of file: FUNCTION
%Programer: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the

%--------------------------------------------------------------------------

function [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
    Bleft,Bright,overedgecoord,normk,Fg,mapinv,maptransm,mapknownvec,...
    pointedge,bodyterm,Hesq,Kde,Kn,Kt,Ded,V,N,kmap,nflag,parameter,weightDMP,...
    nflagface,p_old,contnorm,weight,s,transmvecleftcon,transmvecrightcon,...
    knownvecleftcon,knownvecrightcon,storeinvcon,Bleftcon,Brightcon,Fgcon,...
    mapinvcon,maptransmcon,mapknownveccon,pointedgecon,bodytermcon,Kdec,...
    Knc,Ktc,Dedc,wightc,sc,weightDMPc,dparameter,nflagnoc,nflagfacec,Con,...
    lastimelevel,lastimeval,gravrate,source] = preMPFA(kmap,klb,dmap,MM,h)
%Define global parameters:
global pmethod elem interptype phasekey keygravity numcase

%Obtain the coordinate of both CENTER and AUXILARY nodes of elements which
%constitute the mash. The AREA of each element is also calculated.
%Message to user:
disp(' ');
disp('---------------------------------------------------');
disp('>> Preprocessing Pressure or Hydraulic head Equation...');
disp(' ');

%Initialize all parameters of output
transmvecleft = 0; transmvecright = 0; knownvecleft = 0; knownvecright = 0;
storeinv = 0; Bleft = 0; Bright = 0; Fg = 0; mapinv = 0; maptransm = 0;
mapknownvec = 0; pointedge = 0; bodyterm = 0; Hesq = 0; Kde = 0; Kn = 0;
Kt = 0; Ded = 0; V = 0; N = 0; parameter=0; weightDMP=0;nflagface=0;p_old=0;
contnorm=0;weight=0;s=0; transmvecleftcon=0; transmvecrightcon=0;
knownvecleftcon=0; knownvecrightcon=0; storeinvcon=0; Bleftcon=0;
Brightcon=0; Fgcon=0; mapinvcon=0; maptransmcon=0; mapknownveccon=0;
pointedgecon=0; bodytermcon=0; Kdec=0;Knc=0;Ktc=0;Dedc=0;wightc=0;sc=0;
weightDMPc=0;nflag=0;dparameter=0; nflagnoc=0;nflagfacec=0;Con=0; gravrate=0;
lastimeval=0;lastimelevel=0;source=0;
%Define parametric variables:
%Parameter Used in Full Pressure Support (FPS)
%"p" quadrature point to flux in the auxilary sub interaction region
p = 1;
%Parameter Used in Full Pressure Support (FPS) and Triangle Pressure
%Support (TPS)
%"q" quadrature point to flux in the sub interaction region
q = 1;
%==========================================================================
% When creating a new case or benchmark, add the physical properties in 
% functions (1), (3), (4) and (5)
%==========================================================================
%Fill the matrix "overedgecoord"
overedgecoord = overedge;
%% (1) Define the norm of permeability or conductivity hidraulic tensor ("normk")
[normk,kmap] = calcnormk(kmap,MM,h);
%Get the length of the edge with non-null Neumann Boundary Condition.
knownboundlength = getknownboundlength(klb);
% source term
%Catch "source" came from "PLUG_sourcefunction"
if (numcase >= 10 && numcase <= 30) || numcase == 1.6 ||...
        numcase==336 || numcase==333 || numcase==335 || numcase==337 ||...
        numcase==338 || numcase==341 || numcase==342 || numcase==347 || ...
        numcase==341.1
    if numcase==341

        [P]=ferncodes_calcfonte;
    elseif numcase==341.1
        [P]=ferncodes_calcfonte_1D;
    end
    %% (2)Catch "source" came from "PLUG_sourcefunction"
    source = PLUG_sourcefunction(P,0,wells);
end

% calculate the weight
if strcmp(pmethod,'mpfad')|| strcmp(pmethod,'nlfvpp')|| strcmp(pmethod,'mpfaql')
    %% (3) adequacao dos flags de contorno
    [nflag,nflagface] = ferncodes_calflag(0);
    %Call another parameters that I don't know.
    [V,N,] = ferncodes_elementface(nflag);
    %It switches according to "interptype"
    switch char(interptype)
        %LPEW 1
        case 'lpew1'
            % calculo dos pesos que correspondem ao LPEW1
            [weight,s] = ferncodes_Pre_LPEW_1(kmap,N);
            %LPEW 2
        case 'lpew2'
            % calculo dos pesos que correspondem ao LPEW2
            [weight,s] = ferncodes_Pre_LPEW_2(kmap,N,zeros(size(elem,1),1),nflagface,nflag);
    end  %End of SWITCH
end
% Physical and geometric parameters in relation to methods based on harmonic points
if strcmp(pmethod,'mpfah')|| strcmp(pmethod,'nlfvh')|| strcmp(pmethod,'nlfvdmp')

    p_old=1e1*ones(size(elem,1),1); % inicializando a pressao para metodo nao linear
    % utilize tolerancia menor que 10^-12 para testes de convergencia
    % Para teste de monotonicidade ou problemas bifásicos utilize 10^-8.
    % # iterações de Picard
    % faces alrededor de um elemento
    [facelement]=ferncodes_elementfacempfaH;
    % calculo dos pontos harmonicos
    [pointarmonic]=ferncodes_harmonicopoint(kmap);
    % calculo dos parametros ou constantes (ksi)
    % temos usado este parametro durante muito tempo em muitos testes
    [parameter,auxface]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
    %% (4) adequacao dos flag no ponto medio do contorno
    nflagface= ferncodes_contflagface;
    %% (5)adequacao dos flags no vertice do contorno
    nflag = ferncodes_calflag(0);
    % calculo dos pesos da interpolacao dos pontos harmonicos
    [weightDMP]=ferncodes_weightnlfvDMP(kmap, elem);

end

% calculation of the gravitational flux
if strcmp(keygravity,'y')
    [vec_gravelem,vec_gravface,]=PLUG_Gfunction;
    if 200<numcase && numcase<300
        [gravrate,]=gravitation(kmap,vec_gravelem,vec_gravface);
    elseif numcase<200
        [gravrate,]=gravitation(kmap,vec_gravelem,vec_gravface);
    end
end
%--------------------------------------------------------------------------
%Calculate the TRANSMISSIBILITY parameters:

%Chose the type of MPFA according "pmethod"
switch char(pmethod)
    %Calculate the transmissibilities from TPFA
    case 'tpfa'
        %[transmvecleft,knownvecleft,Fg,bodyterm] = transmTPFA(kmap,0);
        %nflag = ferncodes_calflag(0);
        nflagface= ferncodes_contflagface;
        %Get preprocessed terms:
        [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem);
        %Calculate the little matrices for MPFA-TPS (Aavatsmark et al., 1998)
    case 'mpfao'
        %Get the initial condition
        [Con,lastimelevel,lastimeval] = applyinicialcond;
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmTPS(kmap,q,knownboundlength);
        %-------------------------------------------------------------------
        if phasekey==3
            %Get the length of the edge with non-null Neumann Boundary Condition.
            knownboundlengthcon = getknownboundlengthcon(klb);
            [transmvecleftcon,transmvecrightcon,knownvecleftcon,knownvecrightcon,storeinvcon,...
                Bleftcon,Brightcon,Fgcon,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
                bodytermcon] = transmTPScon(dmap,q,knownboundlengthcon);
        end
        %Calculate the little matrices for MPFA-FPS (Edwards and Zheng, 2008)
    case 'fps'
        %Get the initial condition
        [Con,lastimelevel,lastimeval] = applyinicialcond;
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmFPS(kmap,q,p,knownboundlength);
        %-------------------------------------------------------------------
        if phasekey==3
            %Get the length of the edge with non-null Neumann Boundary Condition.
            knownboundlengthcon = getknownboundlengthcon(klb);
            [transmvecleftcon,transmvecrightcon,knownvecleftcon,knownvecrightcon,storeinvcon,...
                Bleftcon,Brightcon,Fgcon,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
                bodytermcon] = transmFPScon(dmap,q,p,knownboundlengthcon);
        end
        %Calculate the little matrices for MPFA-Enriched (Chen et al., 2008)
    case 'empfa'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmEnriched(kmap,knownboundlength);

        %Calculate geometrical and physical terms to be used in MPFA-Diamond
    case 'mpfad' %(Gao and Wu, 2010)
        %Get preprocessed terms:
        [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap, elem);
        % for the concentration transport with pressure
        if 200<numcase && numcase<300
            %Get the initial condition
            [Con,lastimelevel,lastimeval] = applyinicialcond;
            % calculate the auxiliary parameters
            [Hesq,Kdec,Knc,Ktc,Dedc,wightc,sc,weightDMPc,dparameter ]=...
                parametersauxiliary(dmap,N);
            % flags boundary conditions
            [nflagnoc,nflagfacec] = ferncodes_calflag_con(lastimeval);
        elseif 350<numcase && numcase<400
            % for the concentration transport with hydraulic head
            %Get the initial condition
            [Con,lastimelevel,lastimeval] = applyinicialcond;
            % calculate the auxiliary parameters
            [Hesq,Kdec,Knc,Ktc,Dedc,wightc,sc,weightDMPc,dparameter ]=...
                parametersauxiliary(dmap,N);
            % flags boundary conditions
            [nflagnoc,nflagfacec] = ferncodes_calflag_con(lastimeval);
        end
        
    case 'mpfaql' % Contreras et al, 2019
        % calculo dos parametros ou constantes (ksi)
        [parameter]=ferncodes_coefficient(kmap);
        % calculo dos pesos DMP
        [weightDMP]=ferncodes_weightnlfvDMP(kmap);
    case 'nlfvpp' % Contreras et al., 2021
        if 200<numcase && numcase<300
            %Get the initial condition
            [Con,lastimelevel,lastimeval] = applyinicialcond;
            % calculate the auxiliary parameters
            [Hesq,Kdec,Knc,Ktc,Dedc,wightc,sc,weightDMPc,dparameter ]=...
                parametersauxiliary(dmap,N);
            % flags boundary conditions
            [nflagnoc,nflagfacec] = ferncodes_calflag_con(lastimeval);
        end
        % contreras et al, 2016
    case 'nlfvh'
        % contreras et al, 2016
        if 200<numcase && numcase<300
            %==================================================================
            % variables relativos a concentracao
            %Get the initial condition
            [Con,lastimelevel,lastimeval] = applyinicialcond;
            % calculoa dos pontos armonicos
            [conpointarmonic]=ferncodes_harmonicopoint(dmap);
            %temos usado para muitos estes o seguinte rutina
            [dparameter,]=ferncodes_coefficientmpfaH(facelement, conpointarmonic,dmap);
            % calculo dos pesos DMP
            [weightDMPc]=ferncodes_weightnlfvDMP(dmap);
            % flags boundary conditions
            [nflagnoc,nflagfacec] = ferncodes_calflag_con(lastimeval);
        end
end  %End of SWITCH

%Message to user:
disp('>> "preMPFA" was finished with success!');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
% calculate the auxiliary parameters
function       [Hesq,Kdec,Knc,Ktc,Dedc,weightc,sc,weightDMPc,dparameter ]=...
    parametersauxiliary(dmap,N)
global interptype elem

% dmap: Valor medio do tensor de difusao molecular, usualmente esse tersor
% consiera-se como homogeneo em todo o dominio
if max(max(dmap))~=0
    elem(:,5)=1;
end
%Get preprocessed terms:
[Hesq,Kdec,Knc,Ktc,Dedc] = ferncodes_Kde_Ded_Kt_Kn(dmap, elem);
%temos usado para muitos estes o seguinte rutina
[dparameter,]=ferncodes_coefficient(dmap,elem);
% calculo dos pesos DMP
[weightDMPc]=ferncodes_weightnlfvDMP(dmap,elem);

%It switches according to "interptype"
switch char(interptype)
    %LPEW 1
    case 'lpew1'
        % calculo dos pesos que correspondem ao LPEW1
        [weightc,sc] = ferncodes_Pre_LPEW_1(dmap,N);
        %LPEW 2
    case 'lpew2'
        % calculo dos pesos que correspondem ao LPEW2
        [weightc,sc] = ferncodes_Pre_LPEW_2(dmap,N,1);

end  %End of SWITCH
end
%--------------------------------------------------------------------------
%FUNCTION "overedgesection"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordsection] = overedgesection(edgematrix)
%Define global parameters:
global coord;

%Initialize the matrix used in this function
coordsection = zeros(size(edgematrix,1),3);
%This loop swept among the limits stabelished above
%"firstlim" and "lastlim" are parameters which define where the loop begins
%and where its finish

iover = 1:size(edgematrix,1);
coordsection(iover,1:3) = 0.5*(coord(edgematrix(iover,1),:) + ...
    coord(edgematrix(iover,2),:));
end
%--------------------------------------------------------------------------
%Function "overedge"
%--------------------------------------------------------------------------

function [overedgecoord] = overedge
%Define global parameters:
global bedge inedge;

%Initialize the matrix "overedgecoord"
overedgecoord = zeros((size(bedge,1) + size(inedge,1)),3);

%Fill "overedgecoord" (just edge over boundary)
overedgecoord(1:size(bedge,1),:) = overedgesection(bedge);
%Fill the "overedgecoord" rest (just edge inside domain)
%"continedge" is an internal edge's counter
overedgecoord(size(bedge,1) + 1:size(overedgecoord,1),:) = ...
    overedgesection(inedge);
%--------------------------------------------------------------------------
end

