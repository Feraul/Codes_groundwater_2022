
% Esta funcao determina os parametros iniciais, adequacao dos pocos
% injetores e produtores, e condicao inicial hydrologica
function [SS,h_old,MM,wells]=prehydrological

global numcase elem centelem

switch numcase
    case 330

        % initially hydraulic charge
        h_old=100*ones(size(elem,1),1);
        % coeficiente de armazenamento
        SS=0.001;
        % topo
        MM=1;
        % find the well element
        b=find(abs(centelem(:,1)-500)<1e-9 & abs(centelem(:,2)-500)<1e-9);
        wells(1,1)=b;
        wells(1,2)=2;
        wells(1,3)=1;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-10000;
end

end