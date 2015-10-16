function teta = multilayer_fnteta_vgm( h, Psh, nodes )
% teta = multilayer_fnteta_vgm( h, Psh, nodes )
% 
% DESCRIPTION
%   It calculates the moisture water content by means of the modified Van
%   Genuchten function. See the SWAP-32 manual, page 28, §2.3 for an
%   introductory description of the modification, which is based on the
%   introduction of a small minimum capillary height causing a minor shift
%   in the soil retention curve.
%   The modification is necessary to prevent for numerical instabilities of
%   the solution scheme.
%   This function has same scope and same code of the *watcon* function in
%   file functions.for, SWAP-32 implementation.
% 
% INPUTs
%   h:      Pressure Head at the compartments of soil grid defined in
%           nodes.
% 
%   Psh:    Structure array of all hydraulic parameters for each
%           compartment of the whole soil grid available during current
%           simulation.
% 
%   nodes:  The list of nodes in range [1, P.nz] from the soil grid of
%           current simulation.
% 
% OUTPUTs
%   teta:   Soil moisture fraction at each compartment in "nodes" list.

%% pre
h_crit      = -1.0d-2;
%% main
head        = h;%h(nodes);
thetar      = Psh.tetar(nodes);
thetas      = Psh.tetas(nodes);
alfamg      = Psh.alfvg(nodes);
n           = Psh.en(nodes);
m           = 1-1./Psh.en(nodes);
h_enpr      = Psh.h_enpr(nodes);
teta        = NaN(numel(nodes),1);

for ii = 1:numel(nodes)

    if h_enpr(ii) > h_crit

        if head(ii) >= 0.0d0
% ---   saturated moisture content
            teta(ii)    = thetas(ii);

        elseif head(ii) > h_crit
            help        = (abs(alfamg(ii)*h_crit))^n(ii);
            help        = (1.0d0 + help) ^ m(ii);
            help        = thetar(ii) + (thetas(ii)-thetar(ii))/help;
            help        = help + (thetas(ii)-help)/(-h_crit)*(head(ii)-h_crit);
            teta(ii)    = min(help,thetas(ii));
        else
% ---   first compute |alpha * h| ** n
            help        = (abs(alfamg(ii)*head(ii))) ^ n(ii);

% ---   add 1 and raise to the power m
            help        = (1.0d0 + help) ^ m(ii);

% ---   now compute theta
            teta(ii)    = thetar(ii)+(thetas(ii)-thetar(ii))/help;
        end

    else

        h105            = 1.05d0 * h_enpr(ii);
        h095            = 0.95d0 * h_enpr(ii);

        if head(ii) >= h095
% ---   saturated moisture content
            teta(ii)    = thetas(ii);
        elseif head(ii) < h095 && head(ii) > h105

            s_enpr      = (1.0d0 + (abs(alfamg(ii)*h_enpr(ii))) ^ n(ii) ) ^ m(ii);
            alphah      = abs(alfamg(ii)*h105);
            theta105    = thetar(ii)+(thetas(ii)-thetar(ii)) / ...
                            ((1.0d0 + alphah ^ n(ii) ) ^ m(ii)) * s_enpr;
            term1       = alphah ^ (n(ii) -1.0d0);
            term2       = term1 * alphah;
            term2       = (1.0d0 + term2) ^ (m(ii) + 1.0d0);
            term2       = (thetas(ii)-thetar(ii)) / term2;
            C105        = abs(-1.0d0 * n(ii) * m(ii) * alfamg(ii)*term2*term1)*s_enpr;

            term3       = (h095 - h105)^3;
            a           = (C105*h095 - C105*h105 + 2*theta105 - 2*thetas(ii))/term3;
            b           = (C105*(h105^2 + h095*h105 -2*h095^2 ) - ...
                            3*(h095 + h105)*(theta105 - thetas(ii)))/term3;
            c           = (h095*(C105*(h095 - h105)*(h095 + 2*h105) + ...
                            6*h105*(theta105 - thetas(ii))))/term3;
            d           = thetas(ii) - a * h095^3 - b * h095^2 - c * h095;
            teta(ii)    = ((a * head(ii) + b ) * head(ii) + c) * head(ii) + d;

        else

% ---   first compute |alpha * h| ^ n
            help        = (abs(alfamg(ii)*head(ii))) ^ n(ii);

% ---   add 1 and raise to the power m
            help        = (1.0d0 + help) ^ m(ii);

% --- For modified VanGenuchten model:
%     - S_enpr: relative saturation at Entry Pressure h_enpr 
            s_enpr      = (1.0d0 + (abs(alfamg(ii)*h_enpr(ii))) ^ n(ii) ) ^ m(ii);

% ---   now compute theta
            teta(ii)    = thetar(ii)+(thetas(ii)-thetar(ii))/help * s_enpr;
        end
    end
end
%% return
end