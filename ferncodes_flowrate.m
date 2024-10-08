%It is called by "ferncodes_solverpressure.m"

%fun�ao que calcula os fluxos nas arestas internas
%equacoes 28 e 29 (heterogeneo) ou 15 e 16 (homogeneo)

function [flowrate,flowresult,flowratedif] = ferncodes_flowrate(p,pinterp,cinterp,Kde,Ded,Kn,Kt,...
    Hesq,viscosity,nflag,Con,Kdec,Knc,Ktc,Dedc,nflagc,gravrate)



%Define global variables:
global coord  bedge inedge centelem bcflag phasekey ...
    smethod numcase bcflagc keygravity dens;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
flowratedif = zeros(bedgesize + inedgesize,1);
%Swept "bedge"
for ifacont = 1:bedgesize
     if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 || ...
                numcase==248 || numcase==249 ||numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(ifacont,:);
            %It is a Two-phase flow
        else
            visonface = 1;
        end  %End of IF
    elseif 30<numcase && numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(ifacont,:));
    else
        visonface=1;
    end
    
    lef = bedge(ifacont,3);
    
    O=centelem(lef,:); % baricentro do elemento a esuqerda
    B1=bedge(ifacont,1);
    B2=bedge(ifacont,2);
    nor = norm(coord(bedge(ifacont,1),:) - coord(bedge(ifacont,2),:));
    % Termo gravitacional
    if strcmp(keygravity,'y')
        if numcase<200
            % escoamento bifasico oleo-agua
            averagedensity=(viscosity(ifacont,:)*dens')/visonface;
            m=averagedensity*gravrate(ifacont);
        else
            % concentracao soluto-solvente
            m=dens(1,1)*gravrate(ifacont)/visonface;
        end
        
    else
        m=0;
    end
    if bedge(ifacont,5) < 200 % se os n�s esteverem na fronteira de DIRICHLET
        c1 = nflag(bedge(ifacont,1),2);
        c2 = nflag(bedge(ifacont,2),2);
        A=(Kn(ifacont)/(Hesq(ifacont)*nor));
        flowrate(ifacont) =-A*(((O-coord(B2,:)))*(coord(B1,:)-coord(B2,:))'*c1+...
            (O-coord(B1,:))*(coord(B2,:)-coord(B1,:))'*c2-(nor^2)*p(lef))-(c2-c1)*Kt(ifacont);
        
        flowrate(ifacont) = visonface*flowrate(ifacont)-visonface*m;
        
    else
        x = logical(bcflag(:,1) == bedge(ifacont,5));
        flowrate(ifacont)= -nor*bcflag(x,2);
    end
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    %% ===================================================================
    % calculo do fluxo dispersivo
    if (200<numcase && numcase<300) || (379<numcase && numcase<400)
        if bedge(ifacont,7)>200
            x=bcflagc(:,1)==bedge(ifacont,7);
            r=find(x==1);
            %flowrate(ifacont,1)= normcont*bcflag(r,2);% testes feitos em todos os
            %problemas monof�sico
            flowratedif(ifacont,1)= nor*bcflagc(r,2);% problema de buckley leverett Bastian
        else
            Ac=(Knc(ifacont)/(Hesq(ifacont)*nor));
            c1aux = nflagc(bedge(ifacont,1),2);
            c2aux = nflagc(bedge(ifacont,2),2);
            flowratedif(ifacont,1)=-Ac*(((O-coord(B2,:)))*(coord(B1,:)-coord(B2,:))'*c1aux+...
                (O-coord(B1,:))*(coord(B2,:)-coord(B1,:))'*c2aux-(nor^2)*Con(lef))-...
                (c2aux-c1aux)*Ktc(ifacont)-visonface*m;
        end
    end
    
end  %End of FOR ("bedge")

%Swept "inedge"
for iface = 1:inedgesize
    if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 ||...
                numcase==248 || numcase==249 || numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(bedgesize + iface,:);
            %It is a Two-phase flow
        else
            visonface = 1;
        end  %End of IF
    elseif 30<numcase && numcase<200 % monofasico:numcase menor que 30
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(bedgesize + iface,:));
    else
        visonface=1;
    end
    
     % contribuicao do termo gravitacional
    if strcmp(keygravity,'y')
        if numcase<200
            % escoamento bifasico oleo-agua
            averagedensity=(viscosity(bedgesize + iface,:)*dens')/visonface;
            m=averagedensity*gravrate(bedgesize + iface,1);
        else
            % concentracao soluto-solvente
            m=dens(1,1)*gravrate(bedgesize + iface,1)/visonface;
        end
        
    else
        m=0;
       
    end
    
    lef = inedge(iface,3); %indice do elemento a direita da aresta i
    rel = inedge(iface,4); %indice do elemento a esquerda da aresta i
    
    p1=pinterp(inedge(iface,1),1);
    p2=pinterp(inedge(iface,2),1);
    
    %calculo das vaz�es
    
    flowrate(bedgesize + iface) =visonface*Kde(iface)*(p(rel)-p(lef)-Ded(iface)*(p2 - p1))-visonface*m;
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    %% ====================================================================
    if (200<numcase && numcase<300) || (379<numcase && numcase<400)
        % calculo do fluxo dispersivo
        
        conno1=cinterp(inedge(iface,1),1);
        conno2=cinterp(inedge(iface,2),1);
        flowratedif(bedgesize + iface) =Kdec(iface)*(Con(rel) - Con(lef) -...
            Dedc(iface)*(conno2 - conno1))-visonface*m;
    end
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%When some multiD schemes are chosen, it is necessary attribute flow rate
%for each half-edge.

%Verify if the scheme is MultiD and which type of that one
if phasekey == 2 && (strcmp(smethod,'mwec') || strcmp(smethod,'mwic') || ...
        strcmp(smethod,'rtmd'))
    %Initialize "auxflowrate"
    auxflowrate = zeros(2*length(flowrate),1);
    %Initialize auxiliary counter
    c = 0;
    %Distribute the flowrate calculated to whole edge in the half-edges.
    for i = 1:length(flowrate)
        auxflowrate(c + 1:c + 2) = 0.5*flowrate(i);
        %Update "c"
        c = c + 2;
    end  %End of FOR
    
    %Finaly, it update "flowrate"
    flowrate = auxflowrate;
    %Clear "auxflowrate"
    clear auxflowrate;
end  %End of IF





