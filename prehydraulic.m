
% Esta funcao determina os parametros iniciais como adequacao dos pocos
% injetores e produtores, condicao inicial hydrologica alguns outros
% parametro fisicos necessarios para rodar os problemas
function [SS,h_old,MM,wells,dt,P]=prehydraulic

global numcase elem centelem
SS=0;
h_old=0;
MM=0;
wells=0;
dt=0;
P=0;
switch numcase
    % The cases 330-333 were obtained of the article -- A local grid-refined
    % numerical groundwater model based on the vertex centered finite
    % volume method
    case 330
        % Case 1: single pumping well in a confined aquifer (Y. Qian et al)
        
        % initially hydraulic charge
        h_old=100*ones(size(elem,1),1);
        % the espeficied storage
        SS=0.001;
        % aquifer thickness
        MM=3;
        % step time
        dt=0.01;
        % find the well element
        %b=find(abs(centelem(:,1)-500)<1e-9 & abs(centelem(:,2)-500)<1e-9);
        b=find((abs(centelem(:,1)-500)/500)<1e-2 & (abs(centelem(:,2)-500)/500)<1e-2);
        wells(1,1)=b;
        wells(1,2)=2;
        wells(1,3)=1;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-10000;
    case 331
       % Case 2: single pumping well in a unconfined aquifer (Y. Qian et al)
        
        % initially hydraulic charge
        h_old=90*ones(size(elem,1),1);
        % the specified yield 
        SS=0.1;
        % aquifer thickness
        MM=100;
        % step time
        dt=0.01;
        % find the well element
        b=find(abs(centelem(:,1)-500)<1e-9 & abs(centelem(:,2)-500)<1e-9);
        wells(1,1)=b;
        wells(1,2)=2;
        wells(1,3)=1;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-40000;
    case 332
        % Case 3: multiple pumping wells in a confined aquifer (Y. Qian et al)
        
        % initially hydraulic charge
        h_old=100*ones(size(elem,1),1);
        % coeficiente de armazenamento especifico
        SS=0.001;
        % espesura do aquifero
        MM=3;
        % step time
        dt=0.01;
        % number of divisions
        const=25; 
        % find the well element
        
        % the flow rate value in the well is divided by number wells
        
        %b1=find(abs(centelem(:,1)-250)<1e-9 & abs(centelem(:,2)-250)<1e-9);
        b1=find((const*9<centelem(:,1)& centelem(:,1)<const*10) & (const*9<centelem(:,2)& centelem(:,2)<const*10));
        %b1=find((abs(centelem(:,1)-250)/250)<1e-2 & (abs(centelem(:,2)-250)/250)<1e-2);
        wells(1,1)=b1;
        wells(1,2)=1;
        wells(1,3)=0;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-10000;
        %-----------------------------------------------------------------
        b2=find((const*30<centelem(:,1)& centelem(:,1)<const*31) & (const*9<centelem(:,2)& centelem(:,2)<const*10));
        %b2=find((abs(centelem(:,1)-750)/750)<1e-2 & (abs(centelem(:,2)-250)/250)<1e-2);
        %b2=find(abs(centelem(:,1)-750)<1e-9 & abs(centelem(:,2)-250)<1e-9);
        wells(2,1)=b2;
        wells(2,2)=1;
        wells(2,3)=0;
        wells(2,4)=0;
        wells(2,5)=0;
        wells(2,6)=-10000;
        %------------------------------------------------------------------
        b3=find((const*30<centelem(:,1)& centelem(:,1)<const*31) & (const*30<centelem(:,2)& centelem(:,2)<const*31));
        %b3=find((abs(centelem(:,1)-750)/750)<1e-2 & (abs(centelem(:,2)-750)/750)<1e-2);
        
        wells(3,1)=b3;
        wells(3,2)=1;
        wells(3,3)=0;
        wells(3,4)=0;
        wells(3,5)=0;
        wells(3,6)=-10000; 
        %------------------------------------------------------------------
        b4=find((const*9<centelem(:,1)& centelem(:,1)<const*10) & (const*30<centelem(:,2)& centelem(:,2)<const*31));
        %b4=find((abs(centelem(:,1)-250)/250)<1e-2 & (abs(centelem(:,2)-750)/750)<1e-2);
        
        wells(4,1)=b4;
        wells(4,2)=1;
        wells(4,3)=0;
        wells(4,4)=0;
        wells(4,5)=0;
        wells(4,6)=-10000;
    case 333
        % Case 4:Two parallel canals (Qian et al)
        elem(:,5)=1:size(elem,1);
        % initially hydraulic charge
        h_old=2*ones(size(elem,1),1);
        % coeficiente de armazenamento especifico
        SS=0.1;
        % espesura do aquifero
        MM=5;
        % step time no considered 
        dt=0.5;
        % precipitation infiltration 
        P=0.002;
end

end