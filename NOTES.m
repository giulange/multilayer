% W.hsurf / W.qsurf
%   I cannot change the two as hqsurf, because they both coexist in
%   multilayer_prog.m (lines 325-326):
%                     W.qsurf         = P.Ep;
%                     W.hsurf         = (3*P.h1(1,P.j)-P.h1(2,P.j))/2;
% Either we change that, or we leave all as it is now.



%% stats on montecarlo: cumulative std during Monte Carlo Simulation ??
% We should compute on the mm repetitions the following statistics for the
% Ti output of the i simulation, where i=[1, ..., M.nvp]:
%   ++mean++    Tave = sum( Ti ) / M.nvp
%   ++bias++    NULL
%   ++ std++    Tstd = sqrt( sum( (Ti - Tave)^2 ) / (M.nvp-1) )
%   ++ mse++    NULL
%   other?
% For instance, for P.C2 we have:
%   Tave = sum( Ti ) / M.nvp
%       where Ti = squeeze( P.C2(:,:,:,mm) )
%   Tstd = sqrt( sum( (Ti - Tave)^2 ) / (M.nvp-1) )
%       where Tave = Ti / M.nvp at current i, and finally the equation is:
%   (Ti - Ti/M.nvp)
%     O.C2 = O.C2 + squeeze( P.C2(:,:,:,mm) );
%     
%     % 
%     squeeze( P.C2(:,:,:,mm) ) - squeeze( P.C2(:,:,:,mm) )/M.nvp

% An example follows:
Ti          = rand(10,1);

Tave        = sum( Ti ) / length( Ti );

Tave_i = 0;
for ii=1:length(Ti)
    Tave_i  = Tave_i + Ti(ii) / length( Ti );
end
fprintf('mean:\n reference=%.4f, manual=%.4f, iterative=%.4f\n',mean(Ti), Tave, Tave_i);

Tstd        = sqrt( sum( (Ti - Tave).^2 ) / (length(Ti)-1) );
% Tstd_i      = zeros(size(Ti));
% for ii=1:length(Ti)
%     Tstd_i  = Tstd_i - Ti ./ (length( Ti )-1);
% end
% iterative is not possible due to the power 2 in expression!!
fprintf('std:\n reference=%.4f, manual=%.4f, iterative=%.4f\n',std(Ti), Tstd, NaN);



