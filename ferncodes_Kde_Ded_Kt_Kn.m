function [Hesq,Kde,Kn,Kt,Ded,flowrateZ,flowresultZ] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem,...
    theta_r,theta_s,pp,alpha,q)

global bedge inedge coord centelem numcase bcflag kmapaux 
% Retorna alguns parâmetros da expressão dos fluxos na face interna e no contorno.
% podemos verificar na pag. 5 e na pag. 6 do paper chines

%Prealocação das matrizes.%

Hesq = zeros(size(bedge,1),1);
Kde = zeros(size(inedge,1),1);
Ded = zeros(size(inedge,1),1);
Kn = zeros(size(bedge,1),1);
Kt = zeros(size(bedge,1),1);
flowresultZ=zeros(size(elem,1),1);
% determinar o flag do nó interior e fronteira de Neumann
K1 = zeros(3,3);
K2 = zeros(3,3);
K = zeros(3,3);

%Loop de arestas de contorno. ("bedge")
for ifacont = 1:size(bedge,1)
    %Determinação do baricentro a elemento à esquerda
    C1 = centelem(bedge(ifacont,3),:);
    lef=bedge(ifacont,3);
    %Determinação das alturas dos elementos à esquerda.
    ve1 = coord(bedge(ifacont,2),:) - coord(bedge(ifacont,1),:); % face
    vm=0.5*(coord(bedge(ifacont,2),:)+coord(bedge(ifacont,1),:));
    ve2 = coord(bedge(ifacont,2),:) - C1; %Do centro esquerdo ao fim da face.

    ce = cross(ve1,ve2); % produto vetorial
    Hesq(ifacont) = norm(ce)/norm(ve1); % altura a relativo as faces do contorno

    if (numcase==432 || numcase==431) && bedge(ifacont,5)<200
        % the fourth column of bedge is vertex flag
        x = logical(bcflag(:,1) == bedge(ifacont,4));
        % label of the vertex
        vertex = bedge(ifacont,1);
        if bcflag(x,1)==101
            if numcase==432
                hbedge = PLUG_bcfunction(vertex,x,1);
            else
                hbedge = bcflag(x,2);
            end
        else
            hbedge = bcflag(x,2);
        end
        if hbedge<0 || hbedge==0
            Se=1/(1+(-alpha*hbedge)^pp)^q;
        else
            Se=1;
        end
        if numcase==431
            coefi=kmapaux(1,2)*(Se^(0.5))*(1-(1-Se^(1/q))^q)^2;
        else
            coefi=kmapaux(1,2)*(Se^0.5)*(1-(1-Se^(pp/(pp-1)))^(q))^2;
        end
        K(1,1) = coefi;
        K(1,2) = 0;
        K(2,1) = 0;
        K(2,2) = coefi;
        %------------------------------------------------------------------
        %Essa é UMA maneira de construir os tensores
        K1(1,1) = kmap(elem(bedge(ifacont,3),5),2);
        K1(1,2) = kmap(elem(bedge(ifacont,3),5),3);
        K1(2,1) = kmap(elem(bedge(ifacont,3),5),4);
        K1(2,2) = kmap(elem(bedge(ifacont,3),5),5);
        Kn1 = (RotH(ve1)'*K1*RotH(ve1))/norm(ve1)^2;
        A=-Kn1/(Hesq(ifacont));
    else
        %Essa é UMA maneira de construir os tensores
        K(1,1) = kmap(elem(bedge(ifacont,3),5),2);
        K(1,2) = kmap(elem(bedge(ifacont,3),5),3);
        K(2,1) = kmap(elem(bedge(ifacont,3),5),4);
        K(2,2) = kmap(elem(bedge(ifacont,3),5),5);
        Kn2 = (RotH(ve1)'*K*RotH(ve1))/norm(ve1)^2;
        A=-Kn2/(Hesq(ifacont));
    end
    %Cálculo das constantes tangenciais e normais
        Kn(ifacont) = (RotH(ve1)'*K*RotH(ve1))/norm(ve1)^2;
        Kt(ifacont) = (RotH(ve1)'*K*(ve1)')/norm(ve1)^2;
    
    %----------------------------------------------------------------------
    if  430<numcase && numcase <450
        if bedge(ifacont,5)==201
            if numcase==433
                flowrateZ(ifacont,1)=A*(norm(ve1))*(vm(1,2)-C1(1,2));
            elseif numcase==435
                flowrateZ(ifacont,1)=0;
            elseif numcase==432 || numcase==431
                flowrateZ(ifacont,1)=0;
            else
                flowrateZ(ifacont,1)=A*(norm(ve1))*(vm(1,2)-C1(1,2));
            end
        elseif bedge(ifacont,5)==202 || bedge(ifacont,5)==203
            flowrateZ(ifacont,1)=0;
        else
            if numcase==435 && bedge(ifacont,5)==101
                flowrateZ(ifacont,1)=-A*(norm(ve1))*(vm(1,2)-C1(1,2));
            else
                if numcase==431 && bedge(ifacont,5)==101
                    flowrateZ(ifacont,1)=A*(norm(ve1))*(vm(1,2)-C1(1,2));
                elseif numcase==432 && bedge(ifacont,5)==102
                    flowrateZ(ifacont,1)=-A*(norm(ve1))*(vm(1,2)-C1(1,2));
                else
                    flowrateZ(ifacont,1)=A*(norm(ve1))*(vm(1,2)-C1(1,2));
                end

            end
        end
        flowresultZ(lef,1)=flowresultZ(lef,1)-flowrateZ(ifacont,1);
    end
    %Adequação dos flags dos nós de Dirichlet
end  %End of FOR ("bedge")

%Loop de arestas internas. ("inedges")
for iface = 1:size(inedge,1),
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    C1 = centelem(inedge(iface,3),:); % baricentro do elemento a esquerda
    C2 = centelem(inedge(iface,4),:); % baricentro do elemento direito

    lef=inedge(inedge(iface,3),3);
    rel=inedge(inedge(iface,3),4);

    vcen = C2 - C1;
    vd1 = coord(inedge(iface,2),:) - coord(inedge(iface,1),:);

    %Determinação das alturas dos centróides dos elementos à direita e à%
    %esquerda.                                                          %

    vd2 = C2 - coord(inedge(iface,1),:);     %Do início da aresta até o
    %centro da célula da direita.
    cd = cross(vd1,vd2);
    H2 = norm(cd)/norm(vd1); % altura a direita

    ve2 = C1 - coord(inedge(iface,1),:);

    ce = cross(vd1,ve2);
    H1 = norm(ce)/norm(vd1); % altura a esquerda

    %Cálculo das constantes.%
    %A segunda entrada será tal que: 1=dir, 2=esq.

    %Essa é UMA maneira de construir os tensores.
    %Permeability on the Left

    K1(1,1) = kmap(elem(inedge(iface,3),5),2);
    K1(1,2) = kmap(elem(inedge(iface,3),5),3);
    K1(2,1) = kmap(elem(inedge(iface,3),5),4);
    K1(2,2) = kmap(elem(inedge(iface,3),5),5);

    %Permeability on the Right

    K2(1,1) = kmap(elem(inedge(iface,4),5),2);
    K2(1,2) = kmap(elem(inedge(iface,4),5),3);
    K2(2,1) = kmap(elem(inedge(iface,4),5),4);
    K2(2,2) = kmap(elem(inedge(iface,4),5),5);

    % calculo das constantes tangenciais e normais em cada face interna
    Kn1 = (RotH(vd1)'*K1*RotH(vd1))/norm(vd1)^2;
    Kt1 = (RotH(vd1)'*K1*(vd1)')/norm(vd1)^2;

    Kn2 = (RotH(vd1)'*K2*RotH(vd1))/norm(vd1)^2;
    Kt2 = (RotH(vd1)'*K2*(vd1)')/norm(vd1)^2;
    % calculo das constantes nas faces internas
    Kde(iface) = -norm(vd1)*((Kn1*Kn2))/(Kn1*H2 + Kn2*H1);
    % Ded: é uma constante que tem constantes geometricas + contantes
    % tangeciais
    Ded(iface) = (dot(vd1,vcen)/norm(vd1)^2) - ...
        (1/norm(vd1))*((Kt2/Kn2)*H2 + (Kt1/Kn1)*H1);
    %----------------------------------------------------------------------
    if  430<numcase && numcase <450
    flowrateZ(iface+size(bedge,1),1)=Kde(iface)*(centelem(rel,2)-centelem(lef,2));

    flowresultZ(lef,1)=flowresultZ(lef,1)-flowrateZ(iface+size(bedge,1),1);
    flowresultZ(rel,1)=flowresultZ(rel,1)+flowrateZ(iface+size(bedge,1),1);
    end
end  %End of FOR ("inedge")



