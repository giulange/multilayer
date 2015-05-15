%% ponding & runoff & hsurf
if fl_isRunoff || W.SwMacro
    % [boundtop.for,PONDRUNOFF]
    % --- check whether h0max, the max value of pond, yields a runoff
    if h0max <= W.pondmax
        runoff      = 0.0d0;
        pond        = h0max;
        hsurf       = pond;

    else
        runoff = multilayer_runoff(h0max,W.pondmax,W.rsro,W.rsroExp,P.dt);
        % if no runoff occurs: first estimation of pond is OK
        if runoff <= 1.0d-06
            pond	= h0max;
            hsurf   = pond;

        % if exponent of empirical equation (4.2 SWAP-32) is less then 1:
        elseif (runoff >= 1.0d-06) && (W.rsroExp-1 < 1.0d-6)
            p1      = K1max/0.5*P.nodes.dz(1) * P.dt;
            p2      = 1 / ( p1 + 1 + P.dt/W.rsro );
            pond    = p2 * ( pond_jm1 + q0*P.dt - K1max*P.dt + p1*P.h(1) + P.dt*W.rsro*W.pondmax );
            hsurf   = pond;
            runoff  = multilayer_runoff(pond,W.pondmax,W.rsro,W.rsroExp,P.dt);

        % if runoff occurs: find values for pond and runoff iteratively
        else
            p1      = K1max/0.5*P.nodes.dz(1) * P.dt;
            p2      = 1 / ( p1 + 1 );
            % estimation of maximum ponding: ignoring runoff
            h0max   = p2 * (pond_jm1 + q0*P.dt - K1max*P.dt + p1*P.h(1));
            h0min   = 0;
            fl_pond_conv = false;
            for ii = 1:30
                pond    = 0.5*(h0max+h0min);
                runoff  = multilayer_runoff(pond,W.pondmax,W.rsro,W.rsroExp,P.dt);
                h0      = p2 * (pond_jm1 + q0*P.dt - K1max*P.dt + p1*P.h(1) - runoff);
                if pond-h0 < 1.0d-06
                    hsurf = pond;
                    fl_pond_conv = true;
                    break
                else
                    if h0 > pond
                        h0min = pond;
                    else
                        h0max = pond;
                    end
                end
            end
            % if convergence has not been reached: proceed with final value
            if ~fl_pond_conv
                pond    = 0.5*(h0max+h0min);
                hsurf   = pond;
                runoff  = multilayer_runoff(pond,W.pondmax,W.rsro,W.rsroExp,P.dt);
            end
        end
    end
end