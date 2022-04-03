function[analsolution,analpressure]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,x,gamma)
global totaltime numcase elem centelem
switch numcase
    
    case 231
        t=totaltime(2);
        analsolution=(max(satonvertices)./2).*(erfc((x-velmedio.*t)./(2.*sqrt(Dmedio.*t)))+exp((velmedio.*x)./Dmedio).*erfc((x+velmedio.*t)./(2.*sqrt(Dmedio.*t))));
    case 232
        t=totaltime(2);
       analsolution=(max(satonvertices)/2)*(erfc((x - velmedio*t)/(2*sqrt(Dmedio*t)))); 
    case 233 
        t=totaltime(2);
        analsolution=(max(satonvertices)./2).*(erfc((x-velmedio.*t)./(2.*sqrt(Dmedio.*t)))+exp((velmedio.*x)./Dmedio).*erfc((x+velmedio.*t)./(2.*sqrt(Dmedio.*t))));
    case 234
        t=totaltime(2);
        analsolution=(max(satonvertices)./2).*(erfc((x-velmedio.*t)./(2.*sqrt(Dmedio.*t)))+exp((velmedio.*x)./Dmedio).*erfc((x+velmedio.*t)./(2.*sqrt(Dmedio.*t))));
    case 235
        t=totaltime(2);
        analsolution=(max(satonvertices)./2).*(erfc((x-velmedio.*t)./(2.*sqrt(Dmedio.*t)))+exp((velmedio.*x)./Dmedio).*erfc((x+velmedio.*t)./(2.*sqrt(Dmedio.*t))));
    case 236
        t=totaltime(2);
        f1=(x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio);
        f2=(velmedio.*x)./(2*Dmedio);
       ud2=exp(f1).*erfc((x+t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));

       ud1=exp((-x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio)).*erfc((x-t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));
 
      analsolution=(max(satonvertices)./2).*exp(f2).*(ud1-ud2);
    case 237
     t=totaltime(2);
        f1=(x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio);
        f2=(velmedio.*x)./(2*Dmedio);
       ud2=exp(f1).*erfc((x+t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));

       ud1=exp((-x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio)).*erfc((x-t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));
 
      analsolution=(max(satonvertices)./2).*exp(f2).*(ud1-ud2);
    case 238
        t=totaltime(2);
        f1=(x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio);
        f2=(velmedio.*x)./(2*Dmedio);
       ud2=exp(f1).*erfc((x+t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));

       ud1=exp((-x./(2.*Dmedio)).*sqrt((velmedio.^2)+4.*gamma.*Dmedio)).*erfc((x-t.*sqrt((velmedio.^2)+4.*gamma.*Dmedio))./(2.*sqrt(Dmedio.*t)));
 
      analsolution=(max(satonvertices)./2).*exp(f2).*(ud1-ud2);
    case 242
     t=totaltime(2);
     A1= 0.5*(erfc((x-velmedio*t)/(2*sqrt(Dmedio*t)))+        exp((velmedio.*x)./(Dmedio)).*erfc((x+velmedio*t)/(2*sqrt(Dmedio*t))));
     A2= 0.5*(erfc((x-velmedio*(t-1))/(2*sqrt(Dmedio*(t-1))))+exp((velmedio.*x)./(Dmedio)).*erfc((x+velmedio*(t-1))/(2*sqrt(Dmedio*(t-1)))));
     analsolution= A1-A2;
    case 248
        t=totaltime(2);
        for ii=1:size(elem,1)
            
            x=centelem(ii,1);
            y=centelem(ii,2);
            analsolution(ii,1)=sin(0.25*pi*(x+y+2*t));
            analpressure(ii,1)=(0.8/pi)*cos(0.25*pi*(x+y+2*t))+0.5*(x+y);
        end
     
end


end