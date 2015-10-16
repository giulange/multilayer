function C = multilayer_capacity( h, Psh, nodes, dt )
% C = multilayer_capacity( h, P.sh, nodes, dt )
% 
% NOTE
%   See the latest "moiscap" function in functions.for source code of SWAP
%   to update this function.
% 
% DESCRIPTION
%   It computes capillary capacity.
%   It can work on single node or on multiple nodes according to how the i
%   varible is passed in input.
% 
% INPUTs
%   h:      Pressure head at "nodes" compartments.
%           If i is "multiple" defined, it is the vector of potentials at
%           all soil grid nodes (i.e. at all P.nz nodes!).
% 
%   Psh:    The set of hydraulic characteristics of the soil grid at all
%           nodes.
% 
%   nodes:	Current node(s) of soil grid.
%           Two different usages of the function are possible according to
%           the value assigned to i:
%               *one value  --> Capacity at current i node of the soil
%                               grid.
%               *multiple   --> Capacity at all nodes passed in i.
% 
%   dt:     Simulation timestep.
% 
% OUTPUTs
%   C:      Differential water capacity [cm-1] at node(s) "nodes".

%% const
h_crit              = -1.0d-2;
%% read Psh at nodes
head                = h;%h(nodes);
thetar              = Psh.tetar(nodes);
thetas              = Psh.tetas(nodes);
alfamg              = Psh.alfvg(nodes);
n                   = Psh.en(nodes);
m                   = 1-1./Psh.en(nodes);
h_enpr              = Psh.h_enpr(nodes);
%% output
C                   = NaN(numel(nodes),1);
%% main
for ii = 1:numel(nodes)
    
    if h_enpr(ii) > h_crit
        if head(ii) >= 0.0d0
            C(ii)   = dt * 1.0d-7;

        elseif head(ii) > h_crit
            
            term1   = multilayer_fnteta_vgm( h_crit, Psh, nodes(ii) );
            C(ii)   = (thetas(ii)-term1) / (-h_crit);
        else

% ---   use analytical evaluation of capacity
            alphah  = abs(alfamg(ii)*head(ii));

% ---   compute |alpha * h| to the power n-1
            term1   = alphah ^ (n(ii) -1.0d0);

% ---   compute |alpha*h| to the power n
            term2   = term1 * alphah;

% ---   add one and raise to the power m+1
            term2   = (1.0d0 + term2) ^ (m(ii) + 1.0d0);

% ---   divide theta-s minus theta-r by term2
            term2   = (thetas(ii)-thetar(ii)) / term2;

% ---   calculate the differential moisture capacity
            C(ii)   = abs(-1.0d0 * n(ii) * m(ii) * alfamg(ii) * term2 * term1);
        end
    else

        h105        = 1.05d0 * h_enpr(ii);
        h095        = 0.95d0 * h_enpr(ii);

        if head(ii) >= h095
% ---   saturated moisture content
            C(ii)   = dt * 1.0d-7;

        elseif head(ii) < h095 && head(ii) > h105

            s_enpr  = (1.0d0 + (abs(alfamg(ii)*h_enpr(ii))) ^ n(ii) ) ^ m(ii);
            alphah  = abs(alfamg(ii)*h105);
            theta105= thetar(ii)+(thetas(ii)-thetar(ii)) / ...
                        ((1.0d0 + alphah ^ n(ii) ) ^ m(ii)) * s_enpr;
            term1   = alphah ^ (n(ii) -1.0d0);
            term2   = term1 * alphah;
            term2   = (1.0d0 + term2) ^ (m(ii) + 1.0d0);
            term2   = (thetas(ii)-thetar(ii)) / term2;
            C105    = abs(-1.0d0 * n(ii) * m(ii) * alfamg(ii)*term2*term1)*s_enpr;

            term3   = (h095 - h105)^3;
            a       = (C105*h095 - C105*h105 + 2*theta105 - 2*thetas(ii))/term3;
            b       = (C105*(h105^2 + h095*h105 -2*h095^2 ) - ...
                        3*(h095 + h105)*(theta105 - thetas(ii)))/term3;
            c       = (h095*(C105*(h095 - h105)*(h095 + 2*h105) + ...
                        6*h105*(theta105 - thetas(ii))))/term3;

            C(ii)   = ( 3.0d0 * a * head(ii)  + 2.0d0 * b ) * head(ii) + c;
            C(ii)   = max(C(ii), 1.0d-09);

        else
            alphah  = abs(alfamg(ii)*head(ii));
            term1   = alphah ^ (n(ii) -1.0d0);
            term2   = term1 * alphah;
            term2   = (1.0d0 + term2) ^ (m(ii) + 1.0d0);
            term2   = (thetas(ii)-thetar(ii)) / term2;

% --- For modified VanGenuchten model:
%     - S_enpr: relative saturation at Entry Pressure h_enpr 
            s_enpr  = (1.0d0 + (abs(alfamg(ii)*h_enpr(ii))) ^ n(ii) ) ^ m(ii);

            C(ii)   = abs(-1.0d0 * n(ii) * m(ii) * alfamg(ii)*term2*term1) * s_enpr;

        end

    end  

    if (head(ii) > -1.0d0 && C(ii) < (dt * 1.0d-7))
        C(ii)       = dt * 1.0d-7;
    end
end
%% end
return