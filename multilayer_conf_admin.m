W.tolle1                        = 0.0001;
W.dtmin                         = 0.00000001;
W.dtmax                         = 0.5;
W.multmin                       = 1.2;
W.multmax                       = 0.7;

% Nj_shp:                    Set how number of iterations grows
%                               according to increasing W.tmax values. For
%                               small values of W.tmax or for fast
%                               convergence case studies W.Nj_shp can
%                               converge to 2.6, otherwise it must decrease
%                               towards 1.0 (rarely below 1.0). A value of
%                               1.5 is advised for most applications.
W.Nj_shp                     = 1.5;
