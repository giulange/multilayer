function O = multilayer_sub_valorization_depth( I, nz, z )
% O = multilayer_sub_valorization_depth( I, nz, z )
% 
% DESCRIPTION
%   This function translates the definition of DEPTH-dependent parameters
%   made in multilayer_conf.m in MatLab environment variables that can be
%   handled to properly execute the program.
%   In multilayer_conf.m parameters are defined using the couples z and the
%   value of the parameter.
%   Here this definition is draped on the soil geometry defined by user.
% 
% INPUTs
%   I:          Any DEPTH-dependent parameter defined in multilayer_conf.m
%               Example: S.CDE.Cin.NO
% 
%   nz:         The total number of nodes within the soil grid (excluding
%               any other ficticious node outside the pedon).
% 
%   z:          The depth of each node starting from ground. The node is
%               the center of gravity of the compartment.
% 
% OUTPUT
%   O:          It is the new variable assigned to the P structure array
%               when calling this function.
%               Example: P.CDECinNO

%% main
if isempty( I )
    % read from I_depth.txt file, and then applies the same procedure used
    % for user defined with depth!
    % P.hin = VAR(:,col);
else
    if size( I , 1 )==1
        O   = repmat( I(1,2), nz, 1 );
    else
        O   = NaN( nz, 1 );
        for ii = 1:size( I, 1 )
            O( z(1:end-1) > I(ii,1) ) = I( ii, 2 );
        end
    end
end
% [z,O]
%% end
return