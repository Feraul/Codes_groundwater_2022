%It is called by "preMPFA.m"

function [nflag,nflagface] = ferncodes_calflag(a)
% determinar o flag do nó interior e fronteira de Neumann
global coord bedge bcflag numcase;
nflag = 5000*ones(size(coord,1),2);
nflagface=0;
for ifacont = 1:size(bedge,1)
    if numcase==341 || numcase==341.1
       
       % the fourth column of bedge is vertex flag
        x = logical(bcflag(:,1) == bedge(ifacont,4));
        % label of the vertex
        vertex = bedge(ifacont,1);
        %Second column receives the boundary condition value.
        nflag(vertex,2) = PLUG_bcfunction(vertex,x,a);
        %First column receives the boundary condition flag.
        nflag(vertex,1) = bcflag(x,1); 
        %=======================================================
        % the fifth column of bedge is face flag
        x1 = logical(bcflag(:,1) == bedge(ifacont,5));
        
        % label of the edge
        vertex1 = bedge(ifacont,1:2);
        nflagface(ifacont,2)=PLUG_bcfunction(vertex1,x1,a);
        nflagface(ifacont,1) = bcflag(x1,1); 
    else
        % the fourth column of bedge is vertex flag
        x = logical(bcflag(:,1) == bedge(ifacont,4));
        % label of the vertex
        vertex = bedge(ifacont,1);
        %Second column receives the boundary condition value.
        nflag(vertex,2) = PLUG_bcfunction(vertex,x,a);
        %First column receives the boundary condition flag.
        nflag(vertex,1) = bcflag(x,1);
    end
end  %End of FOR