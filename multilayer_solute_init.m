%% Initialize SOLUTE rate/state variables -- Swap::Solute(1)

% === initialize Solute  ==================================================
% I should add to "cml" solute quantities coming from fertilizer applications!!!  
% --- determine initial solute profile from input concentrations
P.cml(:,1)      = P.CDECinNH;
P.cml(:,2)      = P.CDECinNO;
% --- determine derived solute concentrations
% Freundlich adsorption coefficient, [0..100 cm3 g-1] --> or [cm3 mg-1]
KF              = S.CDE.NX.Kf1;
% Freundlich exponent, [0..10 -]
FREXP           = S.CDE.NX.Kf2;
% TO BE PARAMETERIZED:
CREF            = 1.0;   % Reference solute concentration for Freundlich adsorption, [0..1000 mg cm-3]
% solute balance:
for sl=1:2
    P.cmsy(:,sl)= (P.teta.*P.cml(:,sl) + P.sh.dap.*KF(sl).*CREF.*(P.cml(:,sl)./CREF).^FREXP(sl));
    P.samini(sl)= sum( P.cmsy(:,sl) .* P.nodes.dz(1:P.nz) );
end

% Flag indicating that cumulative fluxes should be reset to zero:
flZeroCumu      = true;
% Following SWAP, I initialize with "flZeroCumu = true" and when I enter
% the solute module this switch – after the first dt iteration on water
% transport – would set "P.samini = sampro". But sampro was never defined
% before (the first time will be defined at the end of solute module).
% This means that I have – differently from SWAP – to initialize sampro
% on a correct basis:
P.sampro        = zeros(1,2);

% === in itinere Solute rate/state updates ================================
% % I should delete this line, because P.cml already stores
% % the latest concentration computed at j-1.
% P.cml            = squeeze(O.C2(:,P.j-1,mm,:));