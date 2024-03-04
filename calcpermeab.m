function[permeabi,fonte]=calcpermeab

global centelem
%% 2D exponential case: computation of K fields and source terms f
Nmod=100;
varK=0.1;
[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
%[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
phi = phiExpNmod10000(1:Nmod);
C(:,1) = wavenumberExp0Nmod10000(1:Nmod);
C(:,2) = wavenumberExp1Nmod10000(1:Nmod);
KMean = 15;

fonte=zeros(size(centelem,1),1);
permeabi=zeros(size(centelem,1),1);
for i = 1 : size(centelem,1)
    fonte(i,1) = func(centelem(i,1),centelem(i,2),Nmod,KMean,varK,C(:,1),C(:,2),phi);
    permeabi(i,1) = K(centelem(i,1),centelem(i,2),Nmod,KMean,varK,C(:,1),C(:,2),phi);
end
end
function F = func(x,y,Nmod,KMean,varK,C1,C2,phi)

S1 = sum( (-2*pi)*C1.*sin(phi + (2*pi)*(C1*x + C2*y))) ;

S2 = sum( cos(phi + (2*pi)*(C1*x + C2*y) ) ) ;

S3 = sum( (-2*pi)*C2.*sin(phi + (2*pi)*(C1*x + C2*y))) ;

F = -(2*KMean*exp(-varK/2)*sqrt(varK*2/Nmod)*S1*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y) ...
    - 5*KMean*exp(-varK/2)*exp(sqrt(varK*2/Nmod)*S2)*sin(2*x+y) ...
    + KMean*exp(-varK/2)*sqrt(varK*2/Nmod)*S3*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y));

end
function sol = K(x,y,Nmod,KMean,varK,C1,C2,phi)

coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end

