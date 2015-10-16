function Ki = multilayer_conductivity_node( theta, Psh, nodes )
% Ki = multilayer_conductivity_node( x, P.sh, ii )
% 
% DESCRIPTION
%   It computes nodal (non-saturated?) hydraulic conductivity as a function
%   of theta.
%   It can work on single node or multiple nodes according to how the i
%   variable is passed in input.
%   This function has similar scope of "hconduc" function in functions.for,
%   SWAP-32 implementation.
% 
% INPUTs
%   teta:	It can be:
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
%   nodes:	Current node(s) of soil grid.
%           Two different usages of the function are possible according to
%           the value assigned to it:
%               *one value  --> Capacity at current i node of the soil
%                               grid.
%               *multiple   --> Capacity at all nodes passed in i.
% 
% OUTPUTs
%   Ki:     Nodal conductivity at compartment(s) i.

%% const
h_crit          = -1.0d-2;
hconode_vsmall  = 1.0d-10;
%% read Psh at nodes
thetar          = Psh.tetar(nodes);
thetas          = Psh.tetas(nodes);
alfamg          = Psh.alfvg(nodes);
ksatfit         = Psh.k0(nodes);
n               = Psh.en(nodes);
m               = 1-1./Psh.en(nodes);
lambda          = Psh.bita(nodes);
h_enpr          = Psh.h_enpr(nodes);
%% output
Ki              = NaN(numel(nodes),1);
%% main :: modified Van Genuchten

% if ksatexm > 0.0d0
%     flksatexm       = true;
%     ksatexm         = cofg(10)
%     relsatthr       = cofg(11)
%     ksatthr         = cofg(12)
% else
%     flksatexm       = false;
% end

for ii = 1:numel(nodes)

    relsat              = (theta(ii)-thetar(ii))/(thetas(ii)-thetar(ii));

    % if flksatexm && relsat > relsatthr
    %     term1   = (relsat-relsatthr)/(1.0d0-relsatthr);
    %     Ki = term1 * ksatexm + (1.0d0-term1) * ksatthr;
    % else
        if h_enpr(ii) > h_crit
            if relsat < 0.001d0
                Ki(ii)  = hconode_vsmall;
            elseif relsat > (1.0d0 - 1.0d-6)
                Ki(ii)  = ksatfit(ii);
            else
                term1   = ( 1.0d0-relsat^(1.0d0/m(ii)) ) ^ m(ii);
                Ki(ii)  = ksatfit(ii)*(relsat^lambda(ii))* ...
                            (1.0d0-term1)*(1.0d0-term1);
            end
        else
    % --- For modified VanGenuchten model:

    % ---   now compute thetam
            thetam      = thetar(ii)+(thetas(ii)-thetar(ii)) * ...
                        ((1.0d0 + (abs(alfamg(ii)*h_enpr(ii))) ^ n(ii) ) ^ m(ii));

            if relsat < 0.001d0
                Ki(ii)  = hconode_vsmall;
            else
                if theta(ii) >= thetam
                    Ki(ii) = ksatfit(ii);
                else
                    relsatm = (theta(ii)-thetar(ii))/(thetam-thetar(ii));
                    relsat1 = (thetas(ii) - thetar(ii))/(thetam-thetar(ii));
                    term1   = ( 1.0d0 - (relsatm)^ (1.0d0/m(ii)) ) ^ m(ii);
                    term2   = ( 1.0d0 - (relsat1)^ (1.0d0/m(ii)) ) ^ m(ii);
                    Ki(ii)  = ksatfit(ii)*(relsat^lambda(ii)) * ...
                                ((1.0d0-term1)/(1.0d0-term2)) ^ 2.d0;
                end
            end
        end
        Ki(ii) = min(Ki(ii),ksatfit(ii));
    % end
end

% --- in case of frost conditions
% if swfrost==1
%     Ki = Ki * rfcp + hconode_vsmall * (1.0d0 - rfcp);
% end
%% cutted code
% m           = 1-1/Psh.en(ii);
% if x(a) < Psh.h_enpr
%     Sc      = 1./(1+abs(Psh.alfvg(ii).*Psh.h_enpr(ii)).^Psh.en(ii)).^m;
%     Se      = ( 1./(1+abs(Psh.alfvg(ii).*x(a)).^Psh.en(ii)).^m  ) ./ Sc;
% elseif x(a) >= Psh.h_enpr
%     Se      = 1;
% end
% if Se < 1
%     Ki(a)   = Psh.k0(ii) * Se.^Psh.bita(ii).* ( (1-(1-(Se*Sc)^(1/m))^m) / (1-(1-(Sc)^(1/m))^m) )^2;
% elseif Se>=1
%     Ki(a)   = Psh.k0(ii);
% end
%% end
return