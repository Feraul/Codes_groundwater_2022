%It is called by "ferncodes_solvepressure.m"

function [M,I] = ferncodes_globalmatrix(w,s,Kde,Ded,Kn,Kt,Hesq,viscosity,...
    nflag,nflagface,SS,dt,h,MM,gravrate,kmap)
%Define global variables:
global coord elem esurn1 esurn2 bedge inedge centelem bcflag ...
    numcase methodhydro keygravity dens elemarea normals modflowcompared...
    Nmod;

%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I = zeros(size(elem,1),1);
m=0;
jj=1;
% viscosidade ou mobilidade
visonface = 1;
%Evaluate "bedge"
for ifacont = 1:bedgesize
    
    if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 || ...
                numcase==248 || numcase==249 ||numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(ifacont,:);
            %It is a Two-phase flow
        end  %End of IF
    elseif 30<numcase && numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(ifacont,:));
    end
    
    %Get element on the left:
    lef = bedge(ifacont,3);
    %Get another parameters:
    v0 = coord(bedge(ifacont,2),:) - coord(bedge(ifacont,1),:); %fase.
    v1 = centelem(bedge(ifacont,3),:) - coord(bedge(ifacont,1),:);
    v2 = centelem(bedge(ifacont,3),:) - coord(bedge(ifacont,2),:);
    
    %Calculate the nom of the edge
    nor = norm(coord(bedge(ifacont,1),:) - coord(bedge(ifacont,2),:));
    
    % Tratamento do nó nos vértices 2 e 4%
    
    %Dirichlet Boundary
    if bedge(ifacont,5) < 200
        % valor da pressao no contorno
        c1 = nflag(bedge(ifacont,1),2);
        c2 = nflag(bedge(ifacont,2),2);
        if strcmp(modflowcompared,'y')
            elembedge(jj,1)=bedge(ifacont,3);
            elembedge(jj,2)=nflag(bedge(ifacont,2),2);
            jj=jj+1;
        end
        A = -Kn(ifacont)/(Hesq(ifacont)*norm(v0));
        % contribuicao do termo gravitacional
        if strcmp(keygravity,'y')
            if numcase<200
                % escoamento bifasico oleo-agua
                averagedensity=(viscosity(ifacont,:)*dens')/visonface;
                m=averagedensity*gravrate(ifacont);
            else
                % concentracao soluto-solvente
                m=dens(1,1)*gravrate(ifacont)/visonface;
            end
        end
        % average hydraulic head
        % unconfined aquifer
        
        %------------------------------------------------------------------
        % ambos os nos pertenecem ao contorno de Dirichlet
        if nflag(bedge(ifacont,2),1)<200 && nflag(bedge(ifacont,1),1)<200
            
            %montagem da matriz global
            M(lef,lef)=M(lef,lef)-visonface*A*(norm(v0)^2);
            % termo de fonte
            I(lef)=I(lef)-visonface*A*(dot(v2,-v0)*c1+dot(v1,v0)*c2)+visonface*(c2-c1)*Kt(ifacont)+visonface*m;
            
        else
            % quando um dos nos da quina da malha computacional
            % pertence ao contorno de Neumann
            if nflag(bedge(ifacont,1),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-visonface*A*(norm(v0)^2)+visonface*Kt(ifacont)+visonface*A*dot(v2,-v0);
                % termo de fonte
                I(lef)=I(lef)-visonface*A*(dot(v1,v0)*c2)+visonface*(c2)*Kt(ifacont)+visonface*m;
            elseif nflag(bedge(ifacont,2),1)>200
                %montagem da matriz global
                M(lef,lef)=M(lef,lef)-visonface*A*(norm(v0)^2)-visonface*Kt(ifacont)+visonface*A*dot(v1,v0);
                % termo de fonte
                I(lef)=I(lef)-visonface*A*(dot(v2,-v0)*c1)+visonface*(-c1)*Kt(ifacont)+visonface*m;
            end
        end
        %------------------------------------------------------------------
        %Preenchimento
        %
        %         M(bedge(ifacont,3),bedge(ifacont,3)) = M(bedge(ifacont,3),...
        %             bedge(ifacont,3)) - visonface*A*(norm(v0)^2);
        %
        %         %!!!!!!!!!!!!!!!!!!!!Vericar no monofásico o sinal
        %         I(bedge(ifacont,3)) = I(bedge(ifacont,3)) - ...
        %             visonface*(dot(v2,-v0)*c1 + dot(v1,v0)*c2)*A + ...
        %             visonface*(c2 - c1)*Kt(ifacont);
        
        %Neumann boundary
    else
        if numcase==341
            %-----------------------------------------------------------
            
            [phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
            %[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
            phi = phiExpNmod10000(1:Nmod);
            C(:,1) = wavenumberExp0Nmod10000(1:Nmod);
            C(:,2) = wavenumberExp1Nmod10000(1:Nmod);
            KMean = 15;
            aaa=0.5*(coord(bedge(ifacont,1),:) + coord(bedge(ifacont,2),:));
            auxkmap = K(aaa(1,1),aaa(1,2),KMean,C(:,1),C(:,2),phi);
            %----------------------------------------------------------
            
            
            %auxkmap=kmap(lef, 2);
            I(lef) = I(lef)+ normals(ifacont,2)*auxkmap*nflagface(ifacont,2);
            
        else
            x = logical(bcflag(:,1) == bedge(ifacont,5));
            I(lef) = I(lef) + nor*bcflag(x,2);
        end
    end  %End of IF
end  %End of FOR

% end  %End of IF

% contribuição nas faces internas
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
    elseif 30<numcase && numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(bedgesize + iface,:));
    else
        visonface=1;
    end
    
    % pressão prescrita no elemento do poço injetor
    %Contabiliza as contribuições do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
    
    M(inedge(iface,3),inedge(iface,3)) = ...
        M(inedge(iface,3),inedge(iface,3)) - visonface*Kde(iface);
    M(inedge(iface,3),inedge(iface,4)) = ...
        M(inedge(iface,3),inedge(iface,4)) + visonface*Kde(iface);
    M(inedge(iface,4),inedge(iface,4)) = ...
        M(inedge(iface,4), inedge(iface,4)) - visonface*Kde(iface);
    M(inedge(iface,4),inedge(iface,3)) = ...
        M(inedge(iface,4), inedge(iface,3)) + visonface*Kde(iface);
    
    if nflag(inedge(iface,1),1) < 200
        I(inedge(iface,3)) = I(inedge(iface,3)) - visonface*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
        I(inedge(iface,4)) = I(inedge(iface,4)) + visonface*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
    end
    if nflag(inedge(iface,2),1) < 200
        I(inedge(iface,3)) = I(inedge(iface,3)) + visonface*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
        I(inedge(iface,4)) = I(inedge(iface,4)) - visonface*Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
    end
    % quando o nó pertece ao contorno de Neumann
    if nflag(inedge(iface,1),1) == 202 || nflag(inedge(iface,1),1) == 201
        
        I(inedge(iface,3)) = I(inedge(iface,3)) - visonface*Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
        
        I(inedge(iface,4)) = I(inedge(iface,4)) + visonface*Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
    end
    if nflag(inedge(iface,2),1) == 202 || nflag(inedge(iface,2),1) == 201
        
        I(inedge(iface,3)) = I(inedge(iface,3)) + visonface*Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
        I(inedge(iface,4)) = I(inedge(iface,4)) - visonface*Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
    end
    
    %Contabilização das contribuições dos nós que não estão na
    %fronteiras de Dirichlet.
    
    if nflag(inedge(iface,1),1) > 200
        for j = 1:(esurn2(inedge(iface,1) + 1) - esurn2(inedge(iface,1)))
            
            post_cont = esurn2(inedge(iface,1)) + j;
            
            M(inedge(iface,3),esurn1(post_cont)) = M(inedge(iface,3),...
                esurn1(post_cont)) + visonface*Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4),esurn1(post_cont)) = M(inedge(iface,4),...
                esurn1(post_cont)) - visonface*Kde(iface)*Ded(iface)*w(post_cont);
            
        end
    end
    if nflag(inedge(iface,2),1) > 200
        for j = 1:(esurn2(inedge(iface,2) + 1) - esurn2(inedge(iface,2)))
            
            post_cont = esurn2(inedge(iface,2)) + j;
            
            M(inedge(iface,3), esurn1(post_cont)) = M(inedge(iface,3),...
                esurn1(post_cont)) - visonface*Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont)) = M(inedge(iface,4),...
                esurn1(post_cont)) + visonface*Kde(iface)*Ded(iface)*w(post_cont);
        end
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
        I(inedge(iface,3))=I(inedge(iface,3))+visonface*m;
        I(inedge(iface,4))=I(inedge(iface,4))-visonface*m;
    end
end  %End of FOR ("inedge")
if strcmp(modflowcompared,'y')
    for iw = 1:size(elembedge,1)
        M(elembedge(iw,1),:)=0*M(elembedge(iw,1),:);
        M(elembedge(iw,1),elembedge(iw,1))=1;
        I(elembedge(iw,1))=elembedge(iw,2);
    end
end
%==========================================================================
% calcula um problema transiente
if numcase>300
    %
    if numcase~=336 && numcase~=334 && numcase~=335 &&...
            numcase~=337 && numcase~=338 && numcase~=339 &&...
            numcase~=340 && numcase~=341
        if numcase==333 || numcase==331
            coeficiente=dt^-1*SS.*elemarea(:);
        else
            coeficiente=dt^-1*MM*SS.*elemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            M=M+coeficiente.*eye(size(elem,1));
            I=I+coeficiente.*eye(size(elem,1))*h;
        else
            % Crank-Nicolson method
            I=I+coeficiente.*eye(size(elem,1))*h-0.5*M*h;
            M=0.5*M+coeficiente.*eye(size(elem,1));
            
        end
        
    end
end

end
function sol = K(x,y,KMean,C1,C2,phi)
global Nmod varK
coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end
