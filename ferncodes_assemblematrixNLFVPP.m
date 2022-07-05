function [M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,contnorm)
global inedge coord bedge bcflag elem numcase
%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I = zeros(size(elem,1),1);

for ifacont=1:bedgesize
    
    if numcase == 246 || numcase == 245 || numcase==247 || ...
            numcase==248 || numcase==249 || numcase==251
        % vicosity on the boundary edge
        visonface = viscosity(ifacont,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        %I(lef)=I(lef)- normcont*bcflag(r,2); % testes feitos em todos os
        %problemas monofásico
        I(lef)=I(lef)+ normcont*bcflag(r,2); % problema de buckley leverett Bastian
    else
        %% calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        
        alef= visonface*normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
            parameter(1,2,ifacont)*pinterp(parameter(1,4,ifacont)));
        
        Alef=visonface*normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        
        %% implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ Alef;
        I(lef,1)=I(lef,1)+alef;
    end
end

%% Montagem da matriz global

for iface=1:inedgesize
    if numcase == 246 || numcase == 245 || numcase==247 ||...
            numcase==248 || numcase==249 || numcase==251
        % vicosity on the boundary edge
        visonface = viscosity(bedgesize + iface,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
        
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma= sqrt(vd1(1,1)^2+vd1(1,2)^2);
    ifactual=iface+size(bedge,1);
    
    % calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    % tirei o norma*  parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
    %    parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual))
    alef=parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual));
    
    % direita
    % tirei norma*parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
    %    parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual));
    arel= parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual));
    
     mulef=(abs(arel)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
     
     murel=(abs(alef)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
%      mulef=(abs(arel)+contnorm^2)/(abs(alef)+abs(arel)+2*contnorm^2);
     
%      murel=(abs(alef)+contnorm^2)/(abs(alef)+abs(arel)+2*contnorm^2);

    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    
    ARR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ visonface*ALL;
    M(lef,rel)=M(lef,rel)- visonface*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ visonface*ARR;
    M(rel,lef)=M(rel,lef)- visonface*ALL;
end
end