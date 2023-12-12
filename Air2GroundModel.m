function PL= Air2GroundModel( d, f,h)
    a=0.3
    b=500
    thita=atan(h/d);
    hNLoS=20;
    hLoS=1;
    Propability_LOS=1/(1+a*(exp(-b*(thita-a))));
    Propability_Nlos=1-Propability_LOS;
    if Propability_LOS>Propability_Nlos
    PL = 20*log10(f) + 20*log10(d)+ 32.45+hLoS ;
    else
        PL = 20*log10(f) + 20*log10(d)+hNLoS + 32.45+hLoS 
    end
end 