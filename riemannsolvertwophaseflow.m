function [numflux, earlysw]=riemannsolvertwophaseflow(signder_left,signder_right,sign2der_left,...
    sign2der_right,Sright,Sleft,method,bedgesize, inedg,fw,dotvn,dotvg,...
    gama,charvel_rh,dfwdS,dfwdS_rh,dgamadS,mLLF)

switch method
    case 'serna'
        %It use Upwind (Roe)
        if signder_left*signder_right >= 0 && sign2der_left*sign2der_right >= 0 && mLLF==1
            %         %Verify the sign of the characteristic velocity:
            %         %It uses the saturation on the left
            if charvel_rh > 0 || charvel_rh==0
                %Calculate the numerical flux through interface
                numflux = fw(1)*dotvn + gama(1)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sleft;
                
                %It uses the saturation on the right
            else
                
                %Calculate the numerical flux through interface
                numflux = fw(3)*dotvn + gama(3)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sright;
                
            end  %End of IF (Upwind flux)
            
            %It uses the LLF to define the saturation through edge.
        else
            
            alfamax = max(abs([dfwdS*dotvn + dgamadS*dotvg;dfwdS_rh*dotvn] ));
            %alfamax = max(abs(dfwdS*dotvn + dgamadS*dotvg ));
            %Denine the numerical flux
            Fleft = fw(1)*dotvn + gama(1)*dotvg;
            Fright = fw(3)*dotvn + gama(3)*dotvg;
            %Define Local Lax-Friedrichs Flux
            
            LLFlux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
            
            %Calculate the numerical flux through interface using LLF.
            numflux = LLFlux;
            
            earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
        end  %End of IF (type of flux)
    case 'entropyfix'
        %It use Upwind (Roe)
        if signder_left*signder_right >= 0 
            %         %Verify the sign of the characteristic velocity:
            %         %It uses the saturation on the left
            if charvel_rh > 0 || charvel_rh==0
                %Calculate the numerical flux through interface
                numflux = fw(1)*dotvn + gama(1)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sleft;
                
                %It uses the saturation on the right
            else
                
                %Calculate the numerical flux through interface
                numflux = fw(3)*dotvn + gama(3)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sright;
                
            end  %End of IF (Upwind flux)
            
            %It uses the LLF to define the saturation through edge.
        else
            
            alfamax = max(abs([dfwdS*dotvn + dgamadS*dotvg;dfwdS_rh*dotvn] ));
            %alfamax = max(abs(dfwdS*dotvn + dgamadS*dotvg ));
            %Denine the numerical flux
            Fleft = fw(1)*dotvn + gama(1)*dotvg;
            Fright = fw(3)*dotvn + gama(3)*dotvg;
            %Define Local Lax-Friedrichs Flux
            
            LLFlux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
            
            %Calculate the numerical flux through interface using LLF.
            numflux = LLFlux;
            
            earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
        end  %End of IF (type of flux)
    case 'LLF'
        
        alfamax = max(abs(dfwdS*dotvn + dgamadS*dotvg ));
        %Denine the numerical flux
        Fleft = fw(1)*dotvn + gama(1)*dotvg;
        Fright = fw(3)*dotvn + gama(3)*dotvg;
        %Define Local Lax-Friedrichs Flux
        
        LLFlux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
        
        %Calculate the numerical flux through interface using LLF.
        numflux = LLFlux;
        
        earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
    case 'LLFM'
        alfamax = max(abs([dfwdS*dotvn + dgamadS*dotvg;dfwdS_rh*dotvn] ));
        
        %Denine the numerical flux
        Fleft = fw(1)*dotvn + gama(1)*dotvg;
        Fright = fw(3)*dotvn + gama(3)*dotvg;
        %Define Local Lax-Friedrichs Flux
        
        LLFlux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
        
        %Calculate the numerical flux through interface using LLF.
        numflux = LLFlux;
        
        earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
    case 'upwind'
        %         %Verify the sign of the characteristic velocity:
        %         %It uses the saturation on the left
        if charvel_rh > 0 || charvel_rh==0
            %Calculate the numerical flux through interface
            numflux = fw(1)*dotvn + gama(1)*dotvg;
            %Fill "earlysw"
            earlysw(bedgesize + inedg) = Sleft;
            
            %It uses the saturation on the right
        else
            
            %Calculate the numerical flux through interface
            numflux = fw(3)*dotvn + gama(3)*dotvg;
            %Fill "earlysw"
            earlysw(bedgesize + inedg) = Sright;
            
        end  %End of IF (Upwind flux)
    case 'HLL'
        
        if signder_left*signder_right >= 0 && sign2der_left*sign2der_right >= 0
            %Verify the sign of the characteristic velocity:
            %It uses the saturation on the left
            if charvel_rh > 0 || charvel_rh==0
                %Calculate the numerical flux through interface
                numflux = fw(1)*dotvn + gama(1)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sleft;
                
                %It uses the saturation on the right
            else
                
                %Calculate the numerical flux through interface
                numflux = fw(3)*dotvn + gama(3)*dotvg;
                %Fill "earlysw"
                earlysw(bedgesize + inedg) = Sright;
                
            end  %End of IF (Upwind flux)
            
            
        else
            SL=min([dfwdS(1)*dotvn;dfwdS_rh*dotvn;0]);
            
            SR=max([dfwdS(2)*dotvn;dfwdS_rh*dotvn;0]);
            %Denine the numerical flux
            Fleft = fw(1)*dotvn + gama(1)*dotvg;
            Fright = fw(3)*dotvn + gama(3)*dotvg;
            
            HLL = (SR/(SR-SL))*Fleft - (SL/(SR-SL))*Fright + (SR*SL/(SR-SL))*(Sright - Sleft);
            
            %Calculate the numerical flux through interface using LLF.
            numflux = HLL;
        end
        earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
        
end



end