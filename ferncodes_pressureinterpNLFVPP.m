function [pressurinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,nflagno,...
    w,s,Con,nflagc,wightc,sc)
global coord esurn1 esurn2 numcase
% pressões interpolados:pressurinterp
% concentrações interpolados:cinterp

cinterp=0;
%% interpolação das pressões nos vertices
for no=1:size(coord,1)
    nec1=esurn2(no+1)-esurn2(no);
    c1=0;
    auxflag=202; % quando a vazão é diferente de 0
    if nflagno(no,1) >200
        if nflagno(no,1)==auxflag
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                c1=c1+w(esurn2(no)+j)*p(element1);
            end
            % s(no,1) é contribuição do termo fluxo distinto de zero no nó
            % "no", quando o fluxo é zero o s(no,1) é zero.
            c1=c1+s(no,1);
        else
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                c1=c1+w(esurn2(no)+j)*p(element1);
            end
        end
        
    else
        c1=nflagno(no,2);
    end
    
    pressurinterp(no,1)=c1;
    %% interpolacao concentracao nos vertices
    
    if 200<numcase && numcase<300
        c1aux=0;
        auxflagc=202; % quando a vazão é diferente de 0
        if nflagc(no,1) >200
            if nflagc(no,1)==auxflagc
                for j=1:nec1
                    element1=esurn1(esurn2(no)+j);
                    c1aux=c1aux+wightc(esurn2(no)+j)*Con(element1);
                end
                % s(no,1) é contribuição do termo fluxo distinto de zero no nó
                % "no", quando o fluxo é zero o s(no,1) é zero.
                c1aux=c1aux+sc(no,1);
            else
                for j=1:nec1
                    element1=esurn1(esurn2(no)+j);
                    c1aux=c1aux+wightc(esurn2(no)+j)*Con(element1);
                end
            end
            
        else
            c1aux=nflagc(no,2);
        end
        
        cinterp(no,1)=c1aux;
    end
end

end