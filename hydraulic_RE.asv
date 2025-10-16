%--------------------------------------------------------------------------
%Subject: obtain the hydraulic head from the Richards' equation
%Type of file: FUNCTION
%Programer: Fernando Contreras,
%--------------------------------------------------------------------------
%This routine receives geometry and physical data.
function hydraulic_RE(wells,overedgecoord,V,N,Hesq,Kde,Kn,Kt,Ded,kmap,nflag,...
    parameter,h_init,SS,MM,weight,s,dt,gravrate,nflagface,weightDMP,P,...
    p_old,source,theta_s,theta_r,alpha,pp,q,gravresult,flowrateZ,flowresultZ)
%Define global parameters:
global timew  totaltime coord pmethod filepath elem numcase inedge ...
    bedge interptype centelem bcflag ;

%---------------------------------------------------------------------------
% h: representa a carga hidraulica
% h_new: carga hidraulica atualizado a cada passo de tempo
% h_old: carga hidraulica inicial (a partir da condicao inicial)
%--------------------------------------------------------------------------
%Initialize parameters:
hanalit=0;
haux=0;
%"time" is a parameter which add "dt" in each looping (dimentional or adm.)

time = 0;

stopcriteria = 0;
orderintimestep = ones(size(elem,1),1)*0;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
timew = 0;
% inicialization paramenters 
satonvertices=0;producelem=0;h=h_init;Con=0;Kdec=0;Knc=0;nflagc=0;viscosity=1;
count=1;auxkmap=0;mobility=1;Ktc=0;Dedc=0;wightc=0;sc=0;dparameter=0; count_aux=1;
% storage the file vtk s in the time 0
postprocessor(h_init,zeros(size(inedge,1)+size(bedge,1),1),Con,1-Con,...
    count,overedgecoord,orderintimestep,'i',1,kmap(:,2),0);
if numcase==432 
    auxelem22=bedge(57:80,3);
    plot(h_init(auxelem22), kmap(auxelem22,2))
    xlabel('h(m)')
    ylabel('hidraulic conductivity')
end
if numcase==342
    %----------------------------------------------------------------------
    vx=63; % malha quadrilateral ortogonal e distorcida
    %vx=212;% malha triangular nao-estruturada
    haux(1,1)=time;
    haux(1,2)=h_init(vx);
    hanalit(1,1)=time;
    hanalit(1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*time/(3.28*10^-3))));
end
dtaux=dt;
iterinicial=1;
theta_init=thetafunction(h_init,theta_s,theta_r,alpha,pp,q);
sum2=0;
sum1=0;
tic
%Get "hydraulic head" and "flowrate"
while stopcriteria < 100
    % utiliza o metodo TPFA para aproximar a carga hidraulica
    if strcmp(pmethod,'tpfa')
        [h_new,flowrate,] = ferncodes_solvePressure_TPFA(Kde, Kn,...
            nflagface, Hesq,wells,viscosity, Kdec, Knc,nflagc,Con,SS,dt,h,...
            MM,P,time);
        % utiliza o metodo MPFA-D para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfad')
        % Calculate hydraulic head and flowrate using the MPFA with diamond pacth
        [h_new,flowrate,flowresult,flowratedif,dt_opcional,kmap_h,faceaux] = ...
            ferncodes_solverpressure(mobility,wells,Hesq,Kde,Kn,Kt,Ded,...
            nflag,nflagface,weight,s,Con,Kdec,Knc,Ktc,Dedc,nflagc,wightc,...
            sc,SS,dt,h,MM,gravrate,P,kmap,time,N,p_old,source,theta_s,...
            theta_r,alpha,pp,q,iterinicial,gravresult,flowrateZ,flowresultZ);
        dt=dt_opcional;
        if count==1
            sum1=sum(flowresult);
        end
        
        % utiliza o metodo MPFA-H para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfah')
        % Calculate hydraulic head and flowrate using the MPFA with harmonic
        % points
        [h_new,flowrate,]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells,SS,dtaux,h,MM,gravrate,viscosity,P,time);
        % utiliza o metodo NL-TPFA para aproximar a carga hidraulica
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
    count=count+1
    % update the hydraulic head
    h=h_new;
    % delta t auxiliar
    dtaux=dt;
    %----------------------------------------------------------------------
    % case unconfined aquifer
    % case 1 and 4 of the article Qian, et al 2023
    if numcase==333 || numcase==331 || numcase==347
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
        haux(count+1,1)=log10(time);
        haux(count+1,2)=h_new(vx);
        hanalit(count+1,1)=log10(time);
        hanalit(count+1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*(time)/(3.28*10^-3))));
        % a cada passo de tempo atualiza a condicao de contorno do outro da
        % face a direita
        vv=find(bedge(:,5)==102);
        gg=(coord(bedge(vv,1),:) +coord(bedge(vv,2),:))*0.5;
        %hbound=h(gg);
        hbound=3*erfc(gg(:,1)/(2*sqrt(30.5*(time)/(3.28*10^-3))));
        nflagface(vv,2)=hbound;
    elseif numcase==343
        % calculate the element with coordinate (0.5, 0.5) 
         b2=find((0.4<centelem(:,1)& centelem(:,1)<0.56) & ...
            (0.4<centelem(:,2)& centelem(:,2)<0.56));
         h_time(count,1)=time;
         h_time(count,2)=h(b2);
    elseif numcase==432
         %if time<1 || time==1 % Beit Netofa clay
         %bcflag(3,2)=-2+(2.2*(time/1));
         if time<0.0625 || time==0.0625 %silt loam
         bcflag(3,2)=-2+(2.2*(time/0.0625));
         else
            bcflag(3,2)=0.2;
         end
         % modifica a condicao de Dirichlet
           [nflag,nflagface] = ferncodes_calflag(0);
         
    end
    if numcase==431
        p_oldaux1=logical(h<=0 );
        p_oldaux2=logical(h>0);
        
        p_old=10*p_oldaux2-10*p_oldaux1;
        % balanco de massa
    
    theta_n=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    sum2=sum2+ (flowrate(1,1)-flowrate(402,1));
    MBE(count_aux)=100*((sum(theta_n-theta_init)*0.25- sum(flowresult)*time)/(-sum1*time));
    DT(count_aux)=time;
    count_aux=count_aux+1;
    elseif numcase==433
        p_oldaux1=logical(h<=0 );
        p_oldaux2=logical(h>0);
        p_old=25*p_oldaux2-25*p_oldaux1;
    elseif numcase==432 
       p_oldaux1=logical(h<=0 );
        p_oldaux2=logical(h>0);
        p_old=2*p_oldaux2-2*p_oldaux1; 
    elseif numcase==434
        p_oldaux1=logical(h<=0 );
        p_oldaux2=logical(h>0);
        p_old=50*p_oldaux2-50*p_oldaux1; 
    elseif  numcase==435
        p_oldaux1=logical(h<=0 );
        p_oldaux2=logical(h>0);
        
        p_old=5*p_oldaux2-5*p_oldaux1;
    end
    iterinicial=iterinicial+1;
    % storage the vtks and calculate errors
    postprocessor(h,flowrate,Con,1-Con,count,overedgecoord,orderintimestep,...
        'i',1,kmap_h(:,2),time);
    if numcase==432 && count==2
    auxelem22=bedge(49:80,3);
    plot(h(auxelem22), kmap(auxelem22,2))
    xlabel('h(m)')
    ylabel('hidraulic conductivity')
    end
    if stopcriteria <100
        faceaux(:,1)=0;
    end
end
flowrate (faceaux(:,1),1)
if numcase==435
    theta=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    theta_init=thetafunction(h_init,theta_s,theta_r,alpha,pp,q);
    figure(1)
    plot(centelem(:,2),theta)
    xlabel('Elevation (cm)')
    ylabel('\theta(h)')
    
    hold on
    plot(centelem(:,2),theta_init )
    grid
    figure(2)
    plot(centelem(:,2),h)
    grid
    xlabel('Elevantion (cm)')
    ylabel('h(cm)')
    hold on
elseif numcase==431
    
    figure(1)
    theta=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    plot(h, theta)
    xlabel('h(cm)')
    ylabel('\theta(h)')
    grid 

    figure(2)
    plot(h,centelem(:,2))
    hold on

% Warrick solution 11700
% BBB=[-343.291	73.1783
% -170.127	74.5736
% -78.9873	77.5194];

% Warrick solution 23400
 % BBB=[-345.291	60.7752
 % -170.127	62.3256
 % -78.9873	65.8915];
% Warrick solution 46800
BBB=[-343.291	39.8450
-170.127	41.7054
-78.9873	45.8915];
plot(BBB(:,1),BBB(:,2), 's')
xlabel('h (cm)')
ylabel('Elevation(cm)')
hold on
grid
figure(3)
plot(DT, MBE)
xlabel('Time(s)')
ylabel('Massa Balance Error (%)')
grid
end
toc
% plotagem dos graficos em determinados regioes do dominio
plotandwrite(producelem,Con,h,satonvertices,0,0,0,time,overedgecoord,...
    hanalit,haux);
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
% if numcase==342
%     erro=norm(haux(:,2)-hanalit(:,2))
% elseif numcase==343
%    plot(h_time(:,1),h_time(:,2))
%    hold on
%    grid
% end
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
