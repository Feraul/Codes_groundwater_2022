function [mvector]=sourceterm(mvector,elembedge,P,tempo,wells)
global numcase 

if (numcase >= 10 && numcase <= 30) || numcase == 1.6 ||...
        numcase==336 || numcase==333 || numcase==335 || numcase==337 ||...
        numcase==338 || numcase==341 || numcase==342 || numcase==347
    if numcase==341
        [P]=ferncodes_calcfonte;
    end
    %Catch "source" came from "PLUG_sourcefunction"
    source = PLUG_sourcefunction(P,elembedge,tempo,wells);
    %The vector "mvector" is added to vector returned from function below.
    mvector = mvector + source;
end
end 