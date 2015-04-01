%% define internal pars:
err_msg_wrong_par_set       = @(PAR) sprintf('This parameter was set in a wrong way: %s', PAR);
%% check dirs
if ~exist( proj.ipath, 'dir' )
    error('The Project directory does not exist!')
end
if ~exist( proj.opath, 'dir' )
    mkdir( proj.opath )
elseif length(dir(proj.opath))-2>0
    ListOfFiles = dir(proj.opath);
    ListOfFiles = {ListOfFiles(3:end).name};
    warning('The output folder\n\t%s\nalready has this %d files:\n', ...
            proj.opath,length(ListOfFiles) )
    for ii=1:length(ListOfFiles)
        fprintf('\t> %s\n',ListOfFiles{ii})
    end
    
    STOP = input('Type 0 to continue OR 1 to stop execution!\n');
    if STOP, error('You decided to stop the program'), end
    clear ListOfFiles ii STOP
end
%% n-iterations
% Nj:                        Number of iterations as a function 
%                               f(W.tmax,W.dtin,W.dtmax) according to the
%                               W.Nj_shp shape parameter.
P.Nj                        = round( W.tmax / (10^(W.Nj_shp + ...
                                  mean([log10(W.dtin),log10(W.dtmax)]))) );
%% soil geometry
W.sg.sublayers_names        = {'SoilLay','SubLay','hSubLay','hNode','nNodes'};
if W.sg.type == 1
    P.nz = W.sg.regular(1);
elseif W.sg.type == 2
    P.nz = sum(W.sg.sublayers(:,5));
else
    error('Type of soil geometry not properly set!')
end
% other?
% ...
%% W.dtin<W.dtmin
if W.dtin<W.dtmin
    warning('W.dtin=%f < W.dtmin=%f',W.dtin,W.dtmin)
    W.dtin = W.dtmin;
end
%% condizioni al contorno superiore e inferiore
% if (W.itopvar==0 && W.ibotvar==0) || (W.itopvar==1 && W.ibotvar==1)
%     error('W.itopvar=%d cannot be equal to W.ibotvar=%d',W.itopvar,W.ibotvar)
% end
%% trasporto soluti
switch W.isol
    case 0
        % do nothing!

    case 1
        if ~W.iCtopvar==0, error('Something wrong!'), end
        
    case 2
        if ~W.iCtopvar==1, error('Something wrong!'), end
        %Con questo controllo si assume che il primo strato del profilo
        %corrisponda allo strato di interramento del concime azotato
        if and(B.Ctop.dL<W.zint(2),B.Ctop.dL>W.zint(2)) %# --> never TRUE!!!
            error('check B.Ctop.dL-W.zint(2)')
        end
end
%% vegetazione
if ~W.iveg==1, V = NaN; end
%% controllo congruenza W.itopvar, W.iCtopvar, W.iveg)
if or(W.iveg==1,(and(W.isol==2,W.iCtopvar==1)))
    if W.itopvar==0
       error('check (W.itopvar, W.iCtopvar, W.iveg)=1')
    end
end
%% montecarlo

if W.MTCL==0
    M.nvp               = 1;
elseif W.MTCL==1
    if M.nvp<2
        warning('\n\t%s (W.MTCL=%d),\n\t%s M.nvp=%d times!',  ...
            'You set stochastic simulation',W.MTCL, ...
            'but multilayer will run',M.nvp )
    end
    if sum(M.nlay) ~= size(M.data,1)
        error('%s\n  (M.nlay=%d, while M.data size is %d)', ...
            'The number of stochastic simulations is not properly defined',...
            sum(M.nlay), size(M.data,1))
    end
    if M.combinatorial==0
        if sum( min( M.nlay , M.nvp )~=M.nvp )
            M.bound     = min( M.nlay , M.nvp )~=M.nvp;
            fprintf('\nM.nvp = %d\n', M.nvp)
            warning('\nlayer={%d} is limited according to the value of M.nvp!\n', ...
                    find(M.bound)' )
            M.nvp_old   = M.nvp;
            M.nvp       = min( M.nlay );
            warning('A new value is set for M.nvp=%d (instead of %d)', ...
                    M.nvp,M.nvp_old)
        end

        % -----------------------------------------------------------------
        % compute the combinations using the first M.nvp stochastic
        % simulations from each layer in sequence:
        % -----------------------------------------------------------------
        M.combinations = NaN(M.nvp,W.nlay);
        for ii = 1:W.nlay
            if ii==1
                M.combinations(:,1) = 1:M.nvp;
            else
                s=sum(M.nlay(1:ii-1))+1;
                M.combinations(:,ii) = s:s+M.nvp-1;
            end
            
        end
        clear ii s
        % -----------------------------------------------------------------
        
    elseif M.combinatorial==1
        M.nvp_old       = M.nvp;
        M.nvp           = prod( M.nlay );
        
        % -----------------------------------------------------------------
        % compute all possible combinations from W.nlay sets of stochastic
        % simulations:
        % -----------------------------------------------------------------
        varin = cell( 1,W.nlay );
        for ii = 1:W.nlay
            if ii==1, start=1; else start=sum(M.nlay(1:ii-1))+1; end
            varin{ii}=start:start+M.nlay(ii)-1;
        end
        % replicate the grid vectors to compute the full grid:
        [tmp{1:W.nlay}]     = ndgrid( varin{1:3} );
        % concatenate
        M.combinations      = reshape(cat(W.nlay+1,tmp{:}),[],W.nlay);
        clear varin ii start tmp
        % -----------------------------------------------------------------

    else error( err_msg_wrong_par_set( 'W.combinatorial' ) );
    end
else
    error( err_msg_wrong_par_set( 'W.MTCL' ) );
end
if (W.MTCL==0 && M.nvp>1) || (W.MTCL==1 && M.nvp<2)
    error( err_msg_wrong_par_set( 'M.nvp' ) );
end
%% cutted code:

% %% condizioni al contorno superiore e inferiore ==> mettere in _conf.m
% if W.itopvar==1 && W.ibotvar==0
%     B.sel.name          = 'top';
%     B.sel.tqstar        = B.top.thqstar;
%     B.sel.qstar         = B.top.hqstar;
% elseif W.itopvar==0 && W.ibotvar==1
%     B.sel.name          = 'bot';
%     B.sel.tqstar        = B.bot.thqstar;
%     B.sel.qstar         = B.bot.hqstar;
% else % if they are [0,0] OR [1,1]
%     % see multilayer_checkpars.m
% end
% 
% %% trasporto soluti ==> mettere in _conf.m
% switch W.isol
%     case 0
%         % do nothing!
%     case 1
%         % do nothing!
%     case 2 % ==> check if that is what you were thinking about!!
%         B.sel.name      = 'Ctop';
%         B.sel.tqstar    = B.Ctop.tqstar;
%         B.sel.qstar     = NaN; % B.Ctop.qstar does not exist!
% end
% 
% %% vegetazione ==> mettere in _conf.m, solo se non � inutile!!
% if W.iveg==1 % ==> inutile, meglio definire un unico tqstar per tutto il programma!! 
%     B.sel.name          = 'veg';
%     B.sel.tqstar        = V.tqstar;
% %     B.sel.qstar         = B.bot.hqstar; % ??non c'�??
% end
%% ...something else?
%% include checks on EC (in particular EC.matrix!!)
%% end
clear err_msg_wrong_par_set