function kond = fncond( x, P, ii )
% kond = fncond( x, P, ii )
% 
% x=teta:   contenuto d'acqua (ifc=1;3)
% x=h:      potenziale
%           ifc=2 (DURNER);
%           ifc=4 (R&S interacting distributions);
%           ifc=5 (R&S independent distributions [see R&S 1993 eq.21])

if P.ifc(ii)==1
    em              = 1-1/P.en(ii);
    sef             = (x-P.tetar(ii))/(P.tetas(ii)-P.tetar(ii));
    if sef>0.99999
        kond        = P.k0(ii);
    else 
        kond        = P.k0(ii)*(sef)^P.bita(ii)*(1.-(1.-sef^(1/em))^em)^2;
    end 
    
elseif P.ifc(ii)==2
    em              = 1-1/P.en(ii);
    em2             = 1-1/P.en2(ii);
    sef1            = (1+(P.alfvg(ii)*abs(x))^P.en(ii))^-em;
    sef2            = (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^-em2;
    sef             = P.fi(ii)*sef1+(1-P.fi(ii))*sef2;
    if sef>0.99999
        kond        = k0;
    else 
        kr_macr     = P.fi(ii)*P.alfvg(ii)*(1-(1-sef1^(1/em))^em);
        kr_micr     = (1-P.fi(ii))*P.alfvg2(ii)*(1-(1-sef2^(1/em2))^em2);
        kond        = P.k0(ii)*(P.fi(ii)*sef1 + ...
                     (1-P.fi(ii))*sef2)^P.bita(ii)* ( (kr_macr+kr_micr)/...
                     (P.fi(ii)*P.alfvg(ii)+(1-P.fi(ii))*P.alfvg2(ii)) )^2;
    end
elseif P.ifc(ii)==3
    sef             = (x-P.tetar(ii))/(P.tetas(ii)-P.tetar(ii));
    if sef>0.99999
        kond        = P.k0(ii);
    else 
        kond        = P.k0(ii)*(P.fi(ii) * ...
                      exp(-P.bita(ii)*(P.tetas(ii)-x)) + ...
                     (1.-P.fi(ii))*exp(-P.bita2(ii)*(P.tetas(ii)-x)));
    end
    
elseif P.ifc(ii)==4
    em2             = 1-1/P.en2(ii);
    sef1            = (1+P.alfrs(ii)*abs(x))*exp(-P.alfrs(ii)*abs(x));
    sef2            = (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^-em2;
    sef             = P.fi(ii)*sef1 + (1-P.fi(ii))*sef2;
    if sef>0.99999
        kond        = P.k0(ii);
    else
        gmacr       = P.alfrs(ii) * exp(-P.alfrs(ii)*abs(x));        
        gmicr       = P.alfvg2(ii)*P.en2(ii) * betainc(sef2^(1/em2),1,em2);        
        gmacr_0     = P.alfrs(ii);
        gmicr_0     = P.alfvg2(ii)*P.en2(ii);
        kond        = P.k0(ii)*sef^P.bita(ii) * ...
                      (P.fi(ii)*gmacr+(1-P.fi(ii))*gmicr) / ...
                      (P.fi(ii)*gmacr_0+(1-P.fi(ii))*gmicr_0);
    end

elseif P.ifc(ii)==5
    em2             = 1-1/P.en2(ii);
    sef1            = (1+P.alfrs(ii)*abs(x))*exp(-P.alfrs(ii)*abs(x));
    sef2            = (1+(P.alfvg2(ii)*abs(x))^P.en2(ii))^-em2;
    sef             = P.fi(ii)*sef1 + (1-P.fi(ii))*sef2;
    if sef>0.99999
        kond        = P.k0macr(ii);
    else 
        kr_macr     = sef1^P.bita(ii)*exp(-2*P.alfrs(ii)*abs(x));
        kr_micr     = sef2^P.bita(ii)*(1.-(1.-sef2^(1/em2))^em2)^2;
        kond        = P.k0macr(ii)*kr_macr+P.k0(ii)*kr_micr;
    end
end

return