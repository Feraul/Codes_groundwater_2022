function [normk,kmap] = calcnormk(kmap,MM,h)
%Define global parameters:
global elem centelem numcase;

%Initialize "normk" (it is a vector)
normk = zeros(size(centelem,1),1);
%Define the norm of permeability tensor
%Obtain "kmap" for each case
kmap = PLUG_kfunction(kmap,h);
%Swept all elements
for ik = 1:length(normk)
    %Define the material pointer in "elem"
    pointer = elem(ik,5);
    %It catches only the permeability components
    permcompon = [kmap(pointer,2) kmap(pointer,3); ...
        kmap(pointer,4) kmap(pointer,5)];
    %Calculate the norm of tensor
    normk(ik) = norm(permcompon);
end
if 300<numcase && numcase~=336
    % see equation (2), of the article -- A local grid-refined
    % numerical groundwater model based on the vertex centered finite
    % volume method
    if numcase==330 || numcase==332 
        kmap(:,2:5)=MM*kmap(:,2:5);
    elseif numcase==333
        kmap(:,2:5)=kmap(:,2:5);
    elseif numcase==331
        kmap=kmap;
    end
end
end%End of FOR