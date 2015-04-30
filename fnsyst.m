function h2 = fnsyst( P, W )
% h2 = fnsyst( P, W )
% 
% DESCRIPTION
%   Tridiagonal system
% 
% INPUTs
%   P:      ??
%   W:      ??
% 
% OUTPUTs
%   h2:     Output potentials at all nodes.

%% init
h2              = NaN( P.nz, 1);
AA              = NaN( P.nz, 1);
BB              = NaN( P.nz, 1);
CC              = NaN( P.nz, 1);
kp              = NaN( P.nz, 1);
km              = NaN( P.nz, 1);
alpha           = NaN( P.nz, 1);
beta            = NaN( P.nz, 1);
gamma           = NaN( P.nz, 1);
delta           = NaN( P.nz, 1);
ratio           = NaN( P.nz, 1);
%% internodal conductivities km(i) [i.e. (ki-1/2)] and kp(i) [i.e. (ki+1/2)] 

teta_hbot       = fnteta( W.hbot,  P.sh, P.nz );
teta_hsurf      = fnteta( W.hsurf, P.sh, 1    );

% (ki-1/2) -- top node:
if or(P.sh.ifc(1)==1,P.sh.ifc(1)==3)
%     km(1,1)     = ( P.kond(1)    + multilayer_conductivity_node( teta_hsurf,P.sh,1 )      )/2;
    km(1)       = multilayer_conductivity_internode( ...
                    multilayer_conductivity_node( teta_hsurf,P.sh,1 ), ...
                    P.kond(1), W.Kmeth );
else% 2,4,5
%     km(1,1)     = ( P.kond(1)    + multilayer_conductivity_node( W.hsurf,P.sh,1 )         )/2;
    km(1)       = multilayer_conductivity_internode( ...
                    multilayer_conductivity_node( W.hsurf,P.sh,1 ),...
                    P.kond(1), W.Kmeth );
end
% (ki+1/2) -- bottom node:
if or(P.sh.ifc(P.nz)==1,P.sh.ifc(P.nz)==3)
%     kp(P.nz,1)  = ( P.kond(P.nz) + multilayer_conductivity_node( teta_hbot, P.sh,P.nz )   )/2;
    kp(P.nz)    = multilayer_conductivity_internode( P.kond(P.nz), ...
                    multilayer_conductivity_node(teta_hbot,P.sh,P.nz), W.Kmeth );
else% 2,4,5
%     kp(P.nz,1)  = ( P.kond(P.nz) + multilayer_conductivity_node( W.hbot,P.sh,P.nz )       )/2;
    kp(P.nz)    = multilayer_conductivity_internode( P.kond(P.nz), ...
                    multilayer_conductivity_node(W.hbot,P.sh,P.nz), W.Kmeth );
end

% (ki-1/2) -- bottom node:
% km(P.nz,1)      = ( P.kond(P.nz) + P.kond(P.nz-1)       )/2;
km(P.nz)        = multilayer_conductivity_internode( P.kond(P.nz-1), P.kond(P.nz), W.Kmeth );
% (ki+1/2) -- top node:
% kp(1,1)         = ( P.kond(1)    + P.kond(2)            )/2;
kp(1)           = multilayer_conductivity_internode( P.kond(1), P.kond(2), W.Kmeth );

i               = (2:(P.nz-1));
% (ki-1/2) -- intermediate nodes: (Eq. 2.23 SWAP 32 manual, page 35)
% km(i,1)         = ( P.kond(i)    + P.kond(i-1)          )/2;
km(i)           = multilayer_conductivity_internode( P.kond(i-1), P.kond(i), W.Kmeth );
% (ki+1/2) -- intermediate nodes: (Eq. 2.23 SWAP 32 manual, page 35)
kp(i)           = multilayer_conductivity_internode( P.kond(i), P.kond(i+1), W.Kmeth );

%% risolvere il sistema (vedi swap)
%% intermediate nodes (building the Fi??)
i               = (2:(P.nz-1));
% d(dx/dz)/dz ==> d(dx/dz1)/dz2

% *PROBLEMS with ratio:
%   -should be inverse, that is dz/dt as stated in Eq. 2.28 SWAP 32 manual
%   -^2 is not good for a general implementation using variable compartment
%    heights!
ratio(i,1)      = W.dt ./ ( P.nodes.dz(i).^2 );
% *

alpha(i,1)      = -ratio(i).*km(i);
beta(i,1)       = P.cap(i) + ratio(i).*km(i) + ratio(i).*kp(i);
gamma(i,1)      = -ratio(i).*kp(i);
delta(i,1)      = P.cap(i) .* P.h1(i) + ...
                  P.nodes.dz(i) .* ratio(i) .* ( km(i)-kp(i) )...
                  -W.dt .* P.sink(i);
%% definire W.itbc (top boundary condition) e W.ibbc(bottom boundary condition) (=0 per q e 1 per h) 
%% top node flux (positive upward)
% i               = P.nz-1; % ==> i=1; Is it correct??
i               = 1;
% warning('Antonio wrote i=P.nz-1 (bottom), instead it should be i=1(top)')
ratio(i,1)      = 2*W.dt/(P.nodes.dz(i).^2);
alpha(i,1)      = 0;
gamma(i,1)      = -ratio(i)*kp(i);
if W.itbc==0 % itbc differenzia hsurf da qsurf, visto che adesso abbiamo unificato in hqsurf?
    beta(i,1)   = P.cap(i) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i)*P.h1(i) - ...
                  P.nodes.dz(i)*ratio(i)*(W.qsurf+kp(i)) - W.dt*P.sink(i);
else
    beta(i,1)   = P.cap(i) + ratio(i)*km(i) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i)*P.h1(i) - ...
                  P.nodes.dz(i)*ratio(i)*(km(i)-kp(i)) + ...
                  ratio(i)*km(i)*W.hsurf - W.dt*P.sink(i);
end
%% bottom node flux
ratio(P.nz,1)   = 2*W.dt/(P.nodes.dz(P.nz+1)*P.nodes.dz(P.nz));
alpha(P.nz,1)   = -ratio(P.nz)*km(P.nz);
gamma(P.nz,1)   = 0;
if W.ibbc==0
   beta(P.nz,1) = P.cap(P.nz) + ratio(P.nz)*km(P.nz);
   delta(P.nz,1)= P.cap(P.nz)*P.h1(P.nz) + ...
                  P.nodes.dz(P.nz)*ratio(P.nz)*(km(P.nz)+W.qbot) - ...
                  W.dt*P.sink(P.nz);
else
   beta(P.nz,1) = P.cap(P.nz) + ratio(P.nz)*km(P.nz) + ...
                  ratio(P.nz)*kp(P.nz);
   delta(P.nz,1)= P.cap(P.nz)*P.h1(P.nz) + ...
                  P.nodes.dz(P.nz)*ratio(P.nz)*(km(P.nz)-kp(P.nz)) + ...
                  ratio(P.nz)*kp(P.nz)*W.hbot-W.dt*P.sink(P.nz);
end
%% Thomas algorithm
% NOTES:
%   -see    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
%           It's a simplified version of Gaussian elimination (also called
%           Gauss-Jordan elimination).
%   -see    Appendix 4 of SWAP 32, page 237 to read the calculation of
%           Jacobian coefficients.

BB(1)           = gamma(1)/beta(1);
AA(1)           = delta(1)/beta(1);

% metodo della sostituzione
for i=2:1:P.nz
    CC(i)       = beta(i) - alpha(i).*BB(i-1);
    AA(i)       = ( delta(i) - alpha(i).*AA(i-1) ) ./ CC(i);
    BB(i)       = gamma(i) ./ CC(i);
end

h2(P.nz)        = AA(P.nz);
for i=(P.nz-1):-1:1
    h2(i)       = AA(i) - BB(i).*h2(i+1);
end
%% end
return