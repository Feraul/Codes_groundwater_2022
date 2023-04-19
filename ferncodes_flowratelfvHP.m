function [flowrate,flowresult]=ferncodes_flowratelfvHP(parameter,weightDMP,pinterp,p)
global inedge coord bedge bcflag centelem numcase
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
    
for ifacont=1:bedgesize
        
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        %flowrate(ifacont,1)=normcont*bcflag(r,2);% testes feitos em todos os
        %problemas monofásico
        flowrate(ifacont,1)=-normcont*bcflag(r,2);% problema de buckley leverett Bastian
    else
        facelef1=parameter(1,3,ifacont);
        facelef2=parameter(1,4,ifacont);
        % average hydralic head 
        % unconfined aquifer
        if numcase==301
            h_avg=p(lef);
        else
           h_avg=1; 
        end
        flowrate(ifacont,1)= h_avg*normcont*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*p(lef)-...
            parameter(1,1,ifacont)*pinterp(facelef1)-parameter(1,2,ifacont)*pinterp(facelef2));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    
end
% Montagem da matriz global

for iface=1:inedgesize
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
     % average hydraulic head 
     % unconfined aquifer
    if numcase==301
           h_avg=(p(lef)+p(rel))/2;
        else
           h_avg=1; 
    end
    
    % orientation vector 
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
   
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxfacelef1=parameter(1,3,ifactual);
    auxfacelef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxfacerel1=parameter(2,3,ifactual);
    auxfacerel2=parameter(2,4,ifactual);
    % calculo dos fluxo parcial a esquerda
    fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*p(lef)-...
                   parameter(1,1,ifactual)*pinterp(auxfacelef1)-parameter(1,2,ifactual)*pinterp(auxfacelef2));
    % calculo dos fluxo parcial a direita        
    fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*p(rel)-...
                       parameter(2,1,ifactual)*pinterp(auxfacerel1)-parameter(2,2,ifactual)*pinterp(auxfacerel2));
    % calculo do fluxo unico na face
    flowrate(iface+size(bedge,1),1)=h_avg*(murel*fluxesq-mulef*fluxdireit);
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
end
end