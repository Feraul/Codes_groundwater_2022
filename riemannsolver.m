function [numflux, earlysw]=riemannsolver(Sright,Sleft,method,bedgesize, inedg,dotvn,dotvg,charvel_rh,dotdif)

switch method
   
    case 'upwind'
        
%         ve_mais = max([charvel_rh; 0],[],1);
%         ve_menos= min([charvel_rh; 0],[],1);
%     
%         numflux= ve_mais*Sleft + ve_menos*Sright + dotdif;
    
           if charvel_rh > 0 || charvel_rh==0
                 %Calculate the numerical flux through interface
                 numflux = Sleft*dotvn + dotdif ;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sleft;
                
                %It uses the saturation on the right
            else
                
                 %Calculate the numerical flux through interface
                 numflux = Sright*dotvn + dotdif;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sright;
                
            end  %End of IF (Upwind flux)
   
end
end