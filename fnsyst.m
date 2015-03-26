function h2 = fnsyst( P, W )
% h2 = fnsyst( P, W )

%% init
h2              = NaN( W.nz, 1);

%% definire km(i) [� il (ki-1/2)] e il kp(i) [� il (ki+1/2)] 

kp(1,1)         = (P.kond(1,P.j) + P.kond(2,P.j))/2;
teta_hbot       = fnteta( W.hbot,   P, W.nz );
teta_hsurf      = fnteta( W.hsurf, P, 1    );
    
if or(P.ifc==1,P.ifc==3)
    kp(W.nz,1)  = ( P.kond(W.nz,P.j) + fncond(teta_hbot, P,W.nz) )/2;
    km(1,1)     = ( P.kond(1,P.j)    + fncond(teta_hsurf,P,1   ) )/2;
else
    kp(W.nz,1)  = ( P.kond(W.nz,P.j) + fncond(W.hbot,P,W.nz    ) )/2;
    km(1,1)     = ( P.kond(1,P.j)    + fncond(W.hsurf,P,1      ) )/2;
end

i               = (2:(W.nz-1));
km(W.nz,1)      = ( P.kond(W.nz,P.j) + P.kond(W.nz-1,P.j)        )/2;
kp(i,1)         = ( P.kond(i,P.j)    + P.kond(i+1,P.j)           )/2;
km(i,1)         = ( P.kond(i,P.j)    + P.kond(i-1,P.j)           )/2;

%% risolvere il sistema (vedi swap)
%% intermediate nodes
i               = (2:(W.nz-1));
ratio(i,1)      = W.dtin ./ ( P.dz(i).^2 );
alfa(i,1)       = -ratio(i).*km(i);
beta(i,1)       = P.cap(i,P.j) + ratio(i).*km(i) + ratio(i).*kp(i);
gamma(i,1)      = -ratio(i).*kp(i);
delta(i,1)      = P.cap(i,P.j) .* P.h1(i,P.j) + ...
                  P.dz(i) .* ratio(i) .* ( km(i)-kp(i) )...
                  -W.dtin .* P.sink(i,P.j);

%% definire W.itbc (top boundary condition) e W.ibbc(bottom boundary condition) (=0 per q e 1 per h) 
%% top node flux (positive upward)
% i               = W.nz-1; % ==> i=1; Is it correct??
i               = 1;
% warning('Antonio wrote i=W.nz-1 (bottom), instead it should be i=1(top)')
ratio(i,1)      = 2*W.dtin/(P.dz(i).^2);
alfa(i,1)       = 0;
gamma(i,1)      = -ratio(i)*kp(i);
if W.itbc==0 % itbc differenzia hsurf da qsurf, visto che adesso abbiamo unificato in hqsurf?
    beta(i,1)   = P.cap(i,P.j) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i,P.j)*P.h1(i,P.j) - ...
                  P.dz(i)*ratio(i)*(W.qsurf+kp(i)) - W.dtin*P.sink(i,P.j);
else
    beta(i,1)   = P.cap(i,P.j) + ratio(i)*km(i) + ratio(i)*kp(i);
    delta(i,1)  = P.cap(i,P.j)*P.h1(i,P.j) - ...
                  P.dz(i)*ratio(i)*(km(i)-kp(i)) + ...
                  ratio(i)*km(i)*W.hsurf - W.dtin*P.sink(i,P.j);
end

%% bottom node flux
ratio(W.nz,1)   = 2*W.dtin/(P.dzbot*P.dz(W.nz));
alfa(W.nz,1)    = -ratio(W.nz)*km(W.nz);
gamma(W.nz,1)   = 0;
if W.ibbc==0
   beta(W.nz,1) = P.cap(W.nz,P.j) + ratio(W.nz)*km(W.nz);
   delta(W.nz,1)= P.cap(W.nz,P.j)*P.h1(W.nz,P.j) + ...
                  P.dz(W.nz)*ratio(W.nz)*(km(W.nz)+W.qbot) - ...
                  W.dtin*P.sink(W.nz,P.j);
else
   beta(W.nz,1) = P.cap(W.nz,P.j) + ratio(W.nz)*km(W.nz) + ...
                  ratio(W.nz)*kp(W.nz);
   delta(W.nz,1)= P.cap(W.nz,P.j)*P.h1(W.nz,P.j) + ...
                  P.dz(W.nz)*ratio(W.nz)*(km(W.nz)-kp(W.nz)) + ...
                  ratio(W.nz)*kp(W.nz)*W.hbot-W.dtin*P.sink(W.nz,P.j);
end

% algoritmo di Thomas
BB(1,1)         = gamma(1)/beta(1);
AA(1,1)         = delta(1)/beta(1);

% metodo della sostituzione
for i=2:1:W.nz
    CC(i,1)     = beta(i) - alfa(i).*BB(i-1);
    AA(i,1)     = ( delta(i) - alfa(i).*AA(i-1) ) ./ CC(i);
    BB(i,1)     = gamma(i) ./ CC(i);
end

h2(W.nz)        = AA(W.nz);
for i=(W.nz-1):-1:1
    h2(i)       = AA(i) - BB(i).*h2(i+1);
end

return