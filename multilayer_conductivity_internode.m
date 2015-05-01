function Kis2 = multilayer_conductivity_internode( Ki_, Ki, method, dzi_, dzi )
% Kis2 = multilayer_conductivity_internode( Ki_, Ki, method, dzi_, dzi )
% 
% DESCRIPTION
%   This function computes the internodal conductivity according to the
%   passed method amongst the four methods implemented (see SWAP-32
%   manual, page 35, Eqs. 2.23, 2.24, 2.25 and 2.26).
%   An equivalent method is implemented in "hcomean" function,
%   functions.for, SWAP-32.
%   The available methods are:
%       *1  arithmic mean
%       *2  weighted arithmic mean
%       *3  geometric mean
%       *4  weighted geometric mean
%   The values of Ki_ and Ki can be passed only according to the following
%   two sets:
%     -------------------------------
%           Ki_         Ki
%     -------------------------------
%       {   K(i-1),     K(i)    }
%       {   K(i),       K(i+1)	}
%     -------------------------------
% 
%   That is Ki_ is always the upper compartment of i-th compartment defined
%   by Ki.
% 
% INPUTs
%   Ki_:        Nodal conductivity at upper compartment (i.e. at "i-1").
%               [mandatory]
%                   *K(i-1)     if Ki=K(i)
%                   *K(i)       if Ki=K(i+1)
% 
%   Ki:         Nodal conductivity at current "i" compartment.
%               [mandatory]
%                   *K(i)       Ki_=K(i-1)
%                   *K(i+1)     Ki_=K(i)
% 
%   method:     One the above explained methods.
%               [mandatory]
%                   *1  --> arithmic mean
%                   *2  --> weighted arithmic mean
%                   *3  --> geometric mean
%                   *4  --> weighted geometric mean
% 
%   dzi_:       Thickness of compartment "i-1" ("i").
%               This parameter is needed only for methods {2,4}
%               [optional]
% 
%   dzi:        Thickness of compartment "i" ("i+1").
%               This parameter is needed only for methods {2,4}
%               [optional]
% 
% OUTPUTs
%   Kis2:       Internodal conductivity according to the selected method.

%% check
if (method==2 || method==4) && nargin<5
    error('To run weighted methods you need to pass also the thickness of compartments')
end
if (nargin==5) && (dzi_<=0 || dzi<=0)
    error('You cannot define negative thickness for compartments')
end
%% main
switch method
    case 1
        Kis2 = (Ki_ + Ki) /2;
    case 2
        Kis2 = (dzi_*Ki_ + dzi*Ki) / (dzi_+dzi);
    case 3
        Kis2 = (Ki_^0.5) * (Ki^0.5);
    case 4
        Kis2 = (Ki_)^(dzi_/(dzi_+dzi)) * (Ki)^((dzi/(dzi_+dzi)));
    otherwise
        error('You defined a wrong method to calculate the internodal conductivity')
end
%% end
return