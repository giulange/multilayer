function Ki = multilayer_conductivity_node( x, Psh, i )
% Ki = multilayer_conductivity_node( x, P.sh, ii )
% 
% DESCRIPTION
%   It computes nodal (non-saturated?) hydraulic conductivity.
%   It can work on single node or multiple nodes according to how the i
%   variable is passed in input.
%   This function has similar scope of "hconduc" function in functions.for,
%   SWAP-32 implementation.
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
%   Psh:    The set of hydraulic characteristics of the soil grid at all
%           nodes.
% 
%   i:      Current node(s) of soil grid.
%           Two different usages of the function are possible according to
%           the value assigned to i:
%               *one value  --> Capacity at current i node of the soil
%                               grid.
%               *multiple   --> Capacity at all nodes passed in i.
% 
% OUTPUTs
%   Ki:     Nodal conductivity at compartment(s) i.

%% main
Ki                  = NaN(length(i),1);
if ~length(x)==1
    x               = x(i);
end
for a = 1:length(i)
    ii              = i(a);

    switch Psh.ifc(ii)

        case 1% water content
            m           = 1-1/Psh.en(ii);
            % Eq. 2.7 of SWAP 32 manual, page 28:
            Se          = (x(a)-Psh.tetar(ii))/(Psh.tetas(ii)-Psh.tetar(ii));
            if Se>0.99999
                Ki(a)   = Psh.k0(ii);
            else
                % Eq. 2.6 of SWAP 32 manual, page 28:
                % unsaturated hydraulic conductivity:
                Ki(a)   = Psh.k0(ii) * Se^Psh.bita(ii) * (1-(1-Se^(1/m))^m)^2;
            end

        case 3% water content
            Se          = (x(a)-Psh.tetar(ii))/(Psh.tetas(ii)-Psh.tetar(ii));
            if Se>0.99999
                Ki(a)   = Psh.k0(ii);
            else
                Ki(a)   = Psh.k0(ii)*(Psh.fi(ii) * ...
                              exp(-Psh.bita(ii)*(Psh.tetas(ii)-x(a))) + ...
                             (1.-Psh.fi(ii))*exp(-Psh.bita2(ii)*(Psh.tetas(ii)-x(a))));
            end

        case 2% potential (Durner)
            m           = 1-1/Psh.en(ii);
            em2         = 1-1/Psh.en2(ii);
            sef1        = (1+(Psh.alfvg(ii)*abs(x(a)))^Psh.en(ii))^-m;
            sef2        = (1+(Psh.alfvg2(ii)*abs(x(a)))^Psh.en2(ii))^-em2;
            Se          = Psh.fi(ii)*sef1+(1-Psh.fi(ii))*sef2;
            if Se>0.99999
                Ki(a)   = k0;
            else 
                kr_macr = Psh.fi(ii)*Psh.alfvg(ii)*(1-(1-sef1^(1/m))^m);
                kr_micr = (1-Psh.fi(ii))*Psh.alfvg2(ii)*(1-(1-sef2^(1/em2))^em2);
                Ki(a)   = Psh.k0(ii)*(Psh.fi(ii)*sef1 + ...
                             (1-Psh.fi(ii))*sef2)^Psh.bita(ii)* ( (kr_macr+kr_micr)/...
                             (Psh.fi(ii)*Psh.alfvg(ii)+(1-Psh.fi(ii))*Psh.alfvg2(ii)) )^2;
            end

        case 4% potential (R&S interacting distributions)
            em2         = 1-1/Psh.en2(ii);
            sef1        = (1+Psh.alfrs(ii)*abs(x(a)))*exp(-Psh.alfrs(ii)*abs(x(a)));
            sef2        = (1+(Psh.alfvg2(ii)*abs(x(a)))^Psh.en2(ii))^-em2;
            Se          = Psh.fi(ii)*sef1 + (1-Psh.fi(ii))*sef2;
            if Se>0.99999
                Ki(a)   = Psh.k0(ii);
            else
                gmacr   = Psh.alfrs(ii) * exp(-Psh.alfrs(ii)*abs(x(a)));        
                gmicr   = Psh.alfvg2(ii)*Psh.en2(ii) * betainc(sef2^(1/em2),1,em2);        
                gmacr_0 = Psh.alfrs(ii);
                gmicr_0 = Psh.alfvg2(ii)*Psh.en2(ii);
                Ki(a)   = Psh.k0(ii)*Se^Psh.bita(ii) * ...
                              (Psh.fi(ii)*gmacr+(1-Psh.fi(ii))*gmicr) / ...
                              (Psh.fi(ii)*gmacr_0+(1-Psh.fi(ii))*gmicr_0);
            end

        case 5% potential (R&S independent distributions [see R&S 1993 eq.21])
            em2         = 1-1/Psh.en2(ii);
            sef1        = (1+Psh.alfrs(ii)*abs(x(a)))*exp(-Psh.alfrs(ii)*abs(x(a)));
            sef2        = (1+(Psh.alfvg2(ii)*abs(x(a)))^Psh.en2(ii))^-em2;
            Se          = Psh.fi(ii)*sef1 + (1-Psh.fi(ii))*sef2;
            if Se>0.99999
                Ki(a)   = Psh.k0macr(ii);
            else 
                kr_macr = sef1^Psh.bita(ii)*exp(-2*Psh.alfrs(ii)*abs(x(a)));
                kr_micr = sef2^Psh.bita(ii)*(1.-(1.-sef2^(1/em2))^em2)^2;
                Ki(a)   = Psh.k0macr(ii)*kr_macr+Psh.k0(ii)*kr_micr;
            end
    end
end
%% end
return