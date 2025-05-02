function [M,I]=ferncodes_implicitandcranknicolson(M,I,SS,dt,MM,h)
global numcase elemarea methodhydro 

if (300<numcase && numcase<379) && numcase~=347
    %
    if numcase~=336 && numcase~=334 && numcase~=335 &&...
            numcase~=337 && numcase~=338 && numcase~=339 &&...
            numcase~=340 && numcase~=341 && numcase~=380 && numcase~=341.1 
        if numcase==333 || numcase==331 %|| numcase==347
            %para aquifero nao confinado
            coeficiente=dt^-1*SS.*elemarea(:);
        else
            % para quifero confinado
            coeficiente=dt^-1*MM*SS.*elemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            % equacao 30 Qian et al 2023
            M=coeficiente.*eye(size(elem,1))+M;
            I=I+coeficiente.*eye(size(elem,1))*h;
        else
            % Crank-Nicolson method
            % equacao 33 Qian et al 2023
            I=I+(coeficiente.*eye(size(elem,1))-0.5*M)*h;
            M=  (coeficiente.*eye(size(elem,1))+0.5*M);
        end
    end
end
end