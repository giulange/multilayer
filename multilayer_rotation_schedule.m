%% INFO
% This script runs twice per crop:
%   1.\ the first time when the crop STARTS the cycle;
%   2.\ the second time when the crop ENDS  the cycle;
%% PRE
err_msg_wrong_par_set = @(PAR) sprintf('This parameter was set in a wrong way: %s', PAR);
len_units_time_1      = @(ed,sd) (datenum(ed)-datenum(sd)+2)/W.timestep;
len_units_time        = @(ed,sd) (datenum(ed)-datenum(sd)+1)/W.timestep;% +1 to deal with both extremes of interval 
len_units_time_0      = @(ed,sd) (datenum(ed)-datenum(sd)+0)/W.timestep;
%% rotation
ExitMe              = true;% flag for bare soil (ExitMe=true)
%% crop presence/absence
if isempty(V.rotation)
    % I set W.iveg=0 even if user didn't, this way program will run faster:
    W.iveg          = 0;
    init_crop_vars
    return
end
%% unload crop
for ii = 1:size(V.rotation,1)
    if P.tidx == len_units_time_1(V.rotation(ii,2), W.sdate)
        init_crop_vars
        ExitMe      = true;
        break
    end
end
%% load crop
% current day vs start of (any) crop
for ii = 1:size(V.rotation,1)
    if P.tidx == len_units_time(V.rotation(ii,1), W.sdate)
        % load current crop, if any
        eval( ['run ', fullfile( proj.ipath, 'Crop', V.rotation{ii,4} )] )
        ExitMe      = false;
        V.currcrop  = ii;% current selected crop    val={NaN, ii}
        V.lastcrop  = ii;% latest selected crop     val={NaN, ii}
        break
    end
end
%% exit point
if ExitMe
    return
end
%% check parameters
if V.swLAI==1, error('The value of this flag is not implemeted yet!'), end
%% adjust all dates in CROP parameters
% V.LAI(:,1)          = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.LAI(:,1) );
V.Kc(:,1)           = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.Kc(:,1) );
if ~isempty(V.Ke)
    V.Ke(:,1)       = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.Ke(:,1) );
end
if V.swDroot==0
    V.Droot(:,1)    = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.Droot(:,1) );
end
V.irrintervals(:,1) = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.irrintervals(:,1) );
V.irrintervals(:,2) = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.irrintervals(:,2) );
V.irri(:,1)         = dates_from_MD_to_YMD( V.rotation(V.currcrop,1), V.irri(:,1) );
%% check :: interval-time vs irrigation-time
% It is assumed that irrigation-time cannot be INSIDE interval-time!
% Since just the first row in V.irri can define the value for V.irri
% parameters till the end of interval-time (even more rows in V.irri would
% properly define interval-time till its end), it is only necessary to
% check that the date in the first row is not INSIDE the interval-time:

% |---------simulation time--------------------------------|
%       |-------|      |------interval time------------|
%      |-->        |-->      |-->   irrigation-time (with three rows) GOOD
%        |-->      |-->      |-->   irrigation-time (with three rows) BAD
if datenum(V.irri{1,1}) > datenum( V.irrintervals{1,1} )
    error('Check config!\n Irrigation-time cannot be INSIDE interval-time.')
end
%% irrigation [V.irri]
if V.flirri~=0 && size(V.irri,1)==0
    err_msg_wrong_par_set('V.irri')
elseif V.flirri~=0 && isempty(V.irri{1})
    err_msg_wrong_par_set('V.irri')
end

switch V.flirri
    case 0 % no irrigation schedule for current crop
%         P.irri( tRot )      = 0;
%         P.cirr( :,tRot )    = 0;% [NH4+, NO3-]
%         P.irTYPE( tRot )    = 0;
%         P.irDEPTH( tRot )   = 0;
%         P.h_from( tRot )    = 0;
%         P.h_to( tRot )      = 0;
    case 1 % discrete irrigation dates with iDEPTH volumes
        % The schedule can be written in advance with discrete
        % irrigations.
        % When V.flirri=1, the last two columns of the irrigation table
        % are useless!
        % Useful fields are: { Time,iDEPTH,cirr,Type,rDEPTH }
        for ii = 1:size(V.irri,1)
            dat             = len_units_time(V.irri{ii,1},W.sdate);
            P.irri(dat)     = V.irri{ii,2};
            P.cirr(:,dat)   = V.irri{ii,3:4}; % [NH4+, NO3-]
            P.irTYPE(dat)   = V.irri{ii,5};
            P.irDEPTH(dat)  = V.irri{ii,6};
            P.h_from(dat)   = NaN;% 7
            P.h_to(dat)     = NaN;% 8
        end
    case 2 % discrete irrigation dates with h_to thresholds
        % schedule: known,  volumes:unknown
        % Volumes must be calculated in multilayer_irrigation using the
        % h_to threshold.
        % Useful fields are: { Time,iDEPTH,cirr,Type,rDEPTH,h_to }
        for ii = 1:size(V.irri,1)
            dat             = len_units_time(V.irri{ii,1},W.sdate);
            % skip current date if it is outside interval-time:
            if ~intersect(dat,tInt), continue, end
            % write parameters for current date:
            P.irri(dat)     = 0;
            P.cirr(:,dat)   = V.irri{ii,3:4}; % [NH4+, NO3-]
            P.irTYPE(dat)   = V.irri{ii,5};
            P.irDEPTH(dat)  = V.irri{ii,6};
            P.h_from(dat)   = NaN;% 7
            P.h_to(dat)     = V.irri{ii,8};
        end
    case 3 % continuous dates defining h intervals [from,to]
        % schedule:unknown, volumes:unknown
        % Dates and volumes are calculated in multilayer_irrigation!
        % When V.flirri=3 the irrigation schedule and volumes are not
        % known in advance, but the program calculates them according
        % to pressure head critical thresholds set in h_from and h_to
        % of V.irri!
        % This means that during the simulation at a certain point in
        % time(j) and when the day starts (flStartOfDay), the volume of
        % irrigation is updated according to the difference between the
        % two pressure heads thresholds provided.
        % It is important to highlight the interaction of the different
        % "times" within the program:
        %   |------ simulation time ------------------------|
        %       |------ rotation time -----------------|
        %      |------ interval time -|     |------------|
        %        |---- irrigation time -------------|
        % 
        % Useful fields are: { Time,iDEPTH,cirr,Type,rDEPTH,h_from,h_to }

        % rotation-time:
        %   (In multilayer_check_and_load I already check that rotation-time is not
        %    outside simulation-time).
        tRot    = len_units_time( V.rotation(V.currcrop,1), W.sdate ) : ...
                  len_units_time( V.rotation(V.currcrop,2), W.sdate );

        % interval-time:
        %   (interval-time is allowed to be outside rotation-time (to avoid to
        %    modify CROP file for any modification of main config file), but only
        %    dates within rotation-time are activated).
        tInt = [];
        for ii = 1:size(V.irrintervals,1)
            tInt= [tInt,len_units_time( V.irrintervals(ii,1), W.sdate ) : ...
                        len_units_time( V.irrintervals(ii,2), W.sdate ) ]; %#ok<AGROW>
        end
        
        % irrigation-time:
        %   (the dates in V.irri can be: (i) discrete, and crisp dates are defined
        %    by the values of V.irri parameters; (ii) continuous, and intervals od
        %    dates are considered in tIrr).
        % tIrr = [];
        if size(V.irri,1)==1
            tIrr            = len_units_time( V.irri(1,1), W.sdate ) : ...
                              len_units_time( V.rotation(V.currcrop,2), W.sdate );
            % intersect rotation-time vs interval-time vs irrigation-time:
            dat             = intersect( tIrr, intersect(tRot,tInt) );
            P.irri(dat)     = 0;
            P.cirr(:,dat)   = V.irri{1,3:4};% [NH4+, NO3-]
            P.irTYPE(dat)   = V.irri{1,5};
            P.irDEPTH(dat)  = V.irri{1,6};
            P.h_from(dat)   = V.irri{1,7};
            P.h_to(dat)     = V.irri{1,8};
        else
            for ii = 1:size(V.irri,1)-1
                tIrr            = len_units_time( V.irri(ii,1),  W.sdate ) : ...
                              len_units_time( V.irri(ii+1,1),W.sdate );
                dat             = intersect( tIrr, intersect(tRot,tInt) );
                P.irri(dat)     = 0;
                P.cirr(:,dat)   = V.irri{ii,3:4};% [NH4+, NO3-]
                P.irTYPE(dat)   = V.irri{ii,5};
                P.irDEPTH(dat)  = V.irri{ii,6};
                P.h_from(dat)   = V.irri{ii,7};
                P.h_to(dat)     = V.irri{ii,8};
            end
            tIrr            = len_units_time( V.irri(end,1), W.sdate ) : ...
                              len_units_time( V.rotation(V.currcrop,2), W.sdate );
            dat             = intersect( tIrr, intersect(tRot,tInt) );
            P.irri(dat)     = 0;
            P.cirr(:,dat)   = V.irri{end,3:4};% [NH4+, NO3-]
            P.irTYPE(dat)   = V.irri{end,5};
            P.irDEPTH(dat)  = V.irri{end,6};
            P.h_from(dat)   = V.irri{end,7};
            P.h_to(dat)     = V.irri{end,8};
        end
    otherwise
        error('Wrong definition of V.flirri!')
end
%% V.Kc (reduction coeff. of potential evapotraspiration)
if ~numel(V.Kc)
    err_msg_wrong_par_set('V.Kc')
else
    P.Kc = multilayer_sub_valorization_time_cropars( V.Kc, V.rotation(V.currcrop,2), W.sdate, W.timestep, P.Kc );
end
%% V.Ke (reduction coeff. of potential soil evaporation)
if ~numel(V.Ke)
    % calculate following FAO paper 56, Chapter 7
    warning('I substituted V.Ke with B.top.cfbs following SWAP-32 manual!')
else
    P.Ke = multilayer_sub_valorization_time_cropars( V.Ke, V.rotation(V.currcrop,2), W.sdate, W.timestep, P.Ke );
end
%% V.Droot (Rooting depth)
if ~numel(V.Droot)
    err_msg_wrong_par_set('V.Droot')
end

% Verhulst-Pearl logistic growth function (Eq. 2.22, HYDRUS-1D manual):
%   Lr = Lm * ( L0 / (L0+(Lm-L0)*exp(-rt)) )
%   see this paper to include also the "heat unit concept/model":
%       "Modeling of Carbon Dioxide Transport and Production in Soil. 1.
%        Model Development", Simunek & Suarez, Water Resources Research,
%        vol.29, No.2, pages 487-497, 1993.
% VPL                   = @(L0,Lm,r,t) Lm .* (L0 ./ (L0+(Lm-L0).*exp(-r.*t)));
VPL                     = @(Dr,t) Dr(2) .* (Dr(1) ./ (Dr(1)+(Dr(2)-Dr(1)).*exp(-Dr(3).*t)));

if V.swDroot==1%        parameters are given to calculate root depths
    dat                 = len_units_time( V.rotation(V.currcrop,1), W.sdate ) : ...
                          len_units_time( V.rotation(V.currcrop,2), W.sdate );
    P.Droot(dat)        = VPL( V.Droot, 1:numel(dat) );
elseif V.swDroot==0%    values of root depth are prescribed at certain dates
    if ischar(V.Droot{1})% I'm using the MMDD definition instead of DVS
        error('%s\n  %s',...
            'I need to convert definition of time in Droot from MMDD dates to DVS times.',...
            'See how LAI was implemented and apply the same for Droot.')
    end
    if size(V.Droot,1)==1
% --- the same rooting depth is applied over the whole simulation time
        P.Droot         = repmat(V.Droot{1,2},1,W.tmax);
    else
% --- a linear trend of rooting depth is applied (till supplied timestep)
        X1              = len_units_time(V.Droot(:,1),W.sdate);
        for ii = 1:W.tmax
            if isempty( find(X1>ii,1,'first') )
                break
            end
            F_x1        = find(X1<=ii,1,'last');
            y1          = V.Droot{F_x1,2};
            slope       = (V.Droot{F_x1+1,2}-V.Droot{F_x1,2})/(len_units_time(V.Droot{F_x1+1,1},V.Droot{F_x1,1}) -1);
            P.Droot(ii) = y1 + (ii-X1(F_x1))*slope;
        end
% --- a plateau of rooting depth is applied (for timesteps exceding V.Droot definition)
        P.Droot(ii:end) = V.Droot{end,2};
    end
end
%% LAI
if V.swLAI==0
    % Crop Start:
    cropS               = len_units_time( V.rotation(V.currcrop,1), W.sdate );
    % Crop Duration:
    cropDur             = diff(datenum( V.rotation(V.currcrop,1:2) ))+1;
    X                   = round(cropDur * V.LAI(:,1) / 2  *10)/10;
    for ii = 2:length(X)
        % linearization in each DVS interval:
        int             = V.LAI(ii-1,2): (V.LAI(ii,2)-V.LAI(ii-1,2))/(X(ii)-X(ii-1)) :V.LAI(ii,2);
        P.LAI(dat(1)+X(ii-1):cropS+X(ii)-1) = int(2:end);
    end
elseif V.swLAI==1
    % ...to be implemented yet!
end
%% clear workspace
clear dat ii F_x1 X1s slope y1 VPL len_units_time* err_msg_wrong_*
clear tIrr tInt tRot ExitMe cropS cropDur X int
return