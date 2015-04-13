function kond = multilayer_conductivity_node( x, Psh, ii )
% kond = multilayer_conductivity_node( x, P.sh, ii )
% 
% DESCRIPTION
%   ??
% 
% INPUTs
%   x:      It can be:
%               -teta   contenuto d'acqua, with ifc = { 1 OR 3 }
%               -h      potenziale
%                        *ifc=2 (DURNER);
%                        *ifc=4 (R&S interacting distributions);
%                        *ifc=5 (R&S independent distributions [see R&S
%                               1993 eq.21])
% 
%   Psh:    ??
% 
%   ii:     Current node of soil grid.
% 
% OUTPUTs
%   kond:   Nodal conductivity.

%% main
switch Psh.ifc
    
    case 1% water content
        em          = 1-1/Psh.en(ii);
        % Eq. 2.7 of SWAP 32 manual, page 28:
        sef         = (x-Psh.tetar(ii))/(Psh.tetas(ii)-Psh.tetar(ii));
        if sef>0.99999
            kond    = Psh.k0(ii);
        else
            % Eq. 2.6 of SWAP 32 manual, page 28:
            kond    = Psh.k0(ii)*(sef)^Psh.bita(ii)*(1-(1-sef^(1/em))^em)^2;
        end
    
    case 3% water content
        sef         = (x-Psh.tetar(ii))/(Psh.tetas(ii)-Psh.tetar(ii));
        if sef>0.99999
            kond    = Psh.k0(ii);
        else
            kond    = Psh.k0(ii)*(Psh.fi(ii) * ...
                          exp(-Psh.bita(ii)*(Psh.tetas(ii)-x)) + ...
                         (1.-Psh.fi(ii))*exp(-Psh.bita2(ii)*(Psh.tetas(ii)-x)));
        end
        
    case 2% potential (Durner)
        em          = 1-1/Psh.en(ii);
        em2         = 1-1/Psh.en2(ii);
        sef1        = (1+(Psh.alfvg(ii)*abs(x))^Psh.en(ii))^-em;
        sef2        = (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^-em2;
        sef         = Psh.fi(ii)*sef1+(1-Psh.fi(ii))*sef2;
        if sef>0.99999
            kond    = k0;
        else 
            kr_macr = Psh.fi(ii)*Psh.alfvg(ii)*(1-(1-sef1^(1/em))^em);
            kr_micr = (1-Psh.fi(ii))*Psh.alfvg2(ii)*(1-(1-sef2^(1/em2))^em2);
            kond    = Psh.k0(ii)*(Psh.fi(ii)*sef1 + ...
                         (1-Psh.fi(ii))*sef2)^Psh.bita(ii)* ( (kr_macr+kr_micr)/...
                         (Psh.fi(ii)*Psh.alfvg(ii)+(1-Psh.fi(ii))*Psh.alfvg2(ii)) )^2;
        end
    
    case 4% potential (R&S interacting distributions)
        em2         = 1-1/Psh.en2(ii);
        sef1        = (1+Psh.alfrs(ii)*abs(x))*exp(-Psh.alfrs(ii)*abs(x));
        sef2        = (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^-em2;
        sef         = Psh.fi(ii)*sef1 + (1-Psh.fi(ii))*sef2;
        if sef>0.99999
            kond    = Psh.k0(ii);
        else
            gmacr   = Psh.alfrs(ii) * exp(-Psh.alfrs(ii)*abs(x));        
            gmicr   = Psh.alfvg2(ii)*Psh.en2(ii) * betainc(sef2^(1/em2),1,em2);        
            gmacr_0 = Psh.alfrs(ii);
            gmicr_0 = Psh.alfvg2(ii)*Psh.en2(ii);
            kond    = Psh.k0(ii)*sef^Psh.bita(ii) * ...
                          (Psh.fi(ii)*gmacr+(1-Psh.fi(ii))*gmicr) / ...
                          (Psh.fi(ii)*gmacr_0+(1-Psh.fi(ii))*gmicr_0);
        end

    case 5% potential (R&S independent distributions [see R&S 1993 eq.21])
        em2         = 1-1/Psh.en2(ii);
        sef1        = (1+Psh.alfrs(ii)*abs(x))*exp(-Psh.alfrs(ii)*abs(x));
        sef2        = (1+(Psh.alfvg2(ii)*abs(x))^Psh.en2(ii))^-em2;
        sef         = Psh.fi(ii)*sef1 + (1-Psh.fi(ii))*sef2;
        if sef>0.99999
            kond    = Psh.k0macr(ii);
        else 
            kr_macr = sef1^Psh.bita(ii)*exp(-2*Psh.alfrs(ii)*abs(x));
            kr_micr = sef2^Psh.bita(ii)*(1.-(1.-sef2^(1/em2))^em2)^2;
            kond    = Psh.k0macr(ii)*kr_macr+Psh.k0(ii)*kr_micr;
        end
end
%% end
return