function h2 = fnsyst( P, W )
% h2 = fnsyst( P, W )

%% init
h2              = NaN( P.nz, 1);

%% definire km(i) [� il (ki-1/2)] e il kp(i) [� il (ki+1/2)] 

kp(1,1)         = (P.kond(1) + P.kond(2))/2;
teta_hbot       = fnteta( W.hbot,  P.sh, P.nz );
teta_hsurf      = fnteta( W.hsurf, P.sh, 1    );
    
if or(P.sh.ifc==1,P.sh.ifc==3)
    kp(P.nz,1)  = ( P.kond(P.nz) + fncond(teta_hbot, P.sh,P.nz) )/2;
    km(1,1)     = ( P.kond(1)    + fncond(teta_hsurf,P.sh,1   ) )/2;
else
    kp(P.nz,1)  = ( P.kond(P.nz) + fncond(W.hbot,P.sh,P.nz    ) )/2;
    km(1,1)     = ( P.kond(1)    + fncond(W.hsurf,P.sh,1      ) )/2;
end

i               = (2:(P.nz-1));
km(P.nz,1)      = ( P.kond(P.nz) + P.kond(P.nz-1)        )/2;
kp(i,1)         = ( P.kond(i)    + P.kond(i+1)           )/2;
km(i,1)         = ( P.kond(i)    + P.kond(i-1)           )/2;

%% risolvere il sistema (vedi swap)
%% intermediate nodes
i               = (2:(P.nz-1));
ratio(i,1)      = W.dtin ./ ( P.nodes.dz(i).^2 );
alfa(i,1)       = -ratio(i).*km(i);
beta(i,1)       = P.cap(i) + ratio(i).*km(i) + ratio(i).*kp(i);
gamma(i,1)      = -ratio(i).*kp(i);
delta(i,1)      = P.cap(i) .* P.h1(i) + ...
                  P.nodes.dz(i) .* ratio(i) .* ( km(i)-kp(i) )...
                  -W.dtin .* P.sink(i);

%% definire W.itbc (top boundary condition) e W.ibbc(bottom boundary condition) (=0 per q e 1 per h) 
%% top node flux (positive upward)
% i               = P.nz-1; % ==> i=1; Is it correct??
i               = 1;
% warning('Antonio wrote i=P.nz-1 (bottom), instead it should be i=1(top)')
ratio(i,1)      = 2*W.dtin/(P.nodes.dz(i).^2);
alfa(i,1)       = 0;
gamma(i,1)      = -ratio(i)*kp(i);
if W.itbc==0 % itbc differenzia hsurf da qsurf, visto che adesso abbiamo unificato in hqsurf?
    beta(i,1)   = P.cap(i) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i)*P.h1(i) - ...
                  P.nodes.dz(i)*ratio(i)*(W.qsurf+kp(i)) - W.dtin*P.sink(i);
else
    beta(i,1)   = P.cap(i) + ratio(i)*km(i) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i)*P.h1(i) - ...
                  P.nodes.dz(i)*ratio(i)*(km(i)-kp(i)) + ...
                  ratio(i)*km(i)*W.hsurf - W.dtin*P.sink(i);
end

%% bottom node flux
ratio(P.nz,1)   = 2*W.dtin/(P.nodes.dz(P.nz+1)*P.nodes.dz(P.nz));
alfa(P.nz,1)    = -ratio(P.nz)*km(P.nz);
gamma(P.nz,1)   = 0;
if W.ibbc==0
   beta(P.nz,1) = P.cap(P.nz) + ratio(P.nz)*km(P.nz);
   delta(P.nz,1)= P.cap(P.nz)*P.h1(P.nz) + ...
                  P.nodes.dz(P.nz)*ratio(P.nz)*(km(P.nz)+W.qbot) - ...
                  W.dtin*P.sink(P.nz);
else
   beta(P.nz,1) = P.cap(P.nz) + ratio(P.nz)*km(P.nz) + ...
                  ratio(P.nz)*kp(P.nz);
   delta(P.nz,1)= P.cap(P.nz)*P.h1(P.nz) + ...
                  P.nodes.dz(P.nz)*ratio(P.nz)*(km(P.nz)-kp(P.nz)) + ...
                  ratio(P.nz)*kp(P.nz)*W.hbot-W.dtin*P.sink(P.nz);
end

% algoritmo di Thomas
BB(1,1)         = gamma(1)/beta(1);
AA(1,1)         = delta(1)/beta(1);

% metodo della sostituzione
for i=2:1:P.nz
    CC(i,1)     = beta(i) - alfa(i).*BB(i-1);
    AA(i,1)     = ( delta(i) - alfa(i).*AA(i-1) ) ./ CC(i);
    BB(i,1)     = gamma(i) ./ CC(i);
end

h2(P.nz)        = AA(P.nz);
for i=(P.nz-1):-1:1
    h2(i)       = AA(i) - BB(i).*h2(i+1);
end

return