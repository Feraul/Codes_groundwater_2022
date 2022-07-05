function viscosidade=calc_viscosity(Con,nflagfacec)
global bedge inedge visc elemarea 
for i=1:size(bedge,1)+size(inedge,1)
if i<= size(bedge,1)
    if bedge(i,7)<200
        
        concen=nflagfacec(i,2);
        expo=((1-concen)*visc(2)^-0.25+(visc(1)^(-0.25))*concen)^4;
    else
        concen=0;
        expo=((1-concen)*visc(2)^-0.25+(visc(1)^(-0.25))*concen)^4;
    end
    viscosidade(i,1)=expo;   
else
    lef=inedge(i-size(bedge,1),3);
    rel=inedge(i-size(bedge,1),4);
   concenlef=Con(lef);
   concenrel=Con(rel);
   
    expolef=((1-concenlef)*visc(2)^-0.25+(visc(1)^(-0.25))*concenlef)^4;
    exporel=((1-concenrel)*visc(2)^-0.25+(visc(1)^(-0.25))*concenrel)^4;
    
    
    viscosidade(i,1)=(elemarea(lef)*expolef+elemarea(rel)*exporel)/(elemarea(lef)+elemarea(rel));    
    
end


    
end

end