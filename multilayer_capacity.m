function C = multilayer_capacity( x, Psh, i )
% C = multilayer_capacity( x, P.sh, i )
% 
% NOTE
%   See the latest "moiscap" function in functions.for source code of SWAP
%   to update this function.
% 
% DESCRIPTION
%   It computes capillary capacity.
%   It can work on single node or on multiple nodes according to how the i
%   varible is passed in input.
% 
% INPUTs
%   x:      Matrix potential (e.g. pressure head "h").
%           It can be:
%               -teta   contenuto d'acqua, with ifc = { 1 OR 3 }
%               -h      potenziale
%                        *ifc=2 (DURNER);
%                        *ifc=4 (R&S interacting distributions);
%                        *ifc=5 (R&S independent distributions [see R&S
%                               1993 eq.21]).
%           If i is "multiple" defined, it is the vector of potentials at
%           all soil grid nodes (i.e. at all P.nz nodes!).
%   Psh:    The set of hydraulic characteristics of the soil grid at all
%           nodes.
%   i:      Current node(s) of soil grid.
%           Two different usages of the function are possible according to
%           the value assigned to i:
%               *one value  --> Capacity at current i node of the soil
%                               grid.
%               *multiple   --> Capacity at all nodes passed in i.
% 
% OUTPUTs
%   C:      Capacity at node(s) i.

%% main
C                   = NaN(length(i),1);
x                   = x(i);
for a = 1:length(i)
    ii              = i(a);
    df              = Psh.tetas(ii)-Psh.tetar(ii);
    m               = 1-1./Psh.en(ii);
    m2              = 1-1./Psh.en2(ii);
    
    % adjust Capacity according to the value of x:
    if x(a)>=0 
        C(a)        = 0;
        continue
    end

    switch Psh.ifr(ii)
        case 1 
            % see Eq. 2.8, SWAP 32 manual page 28:
            C(a)	= Psh.alfvg(ii).*m.*Psh.en(ii).* abs(Psh.alfvg(ii).*x(a)).^(Psh.en(ii)-1) .* ...
                          df .* (1+abs(Psh.alfvg(ii).*x(a)).^Psh.en(ii)).^(-m-1);
        case 2
            piece1  = -Psh.alfrs(ii).^2.*abs(x(a)).*exp(-Psh.alfrs(ii).*abs(x(a)));
            piece2  = +Psh.alfvg2(ii).*Psh.en2(ii).*m2 .* ...
                          (Psh.alfvg2(ii).*abs(x(a))).^(Psh.en2(ii)-1) .* ...
                          (1+(Psh.alfvg2(ii).*abs(x(a))).^Psh.en2(ii)).^(-1-m2);
            C(a)	= df.*(Psh.fi(ii).*piece1 + (1-Psh.fi(ii)).*piece2);
        case 3
            piece1  = Psh.alfvg(ii).*Psh.en(ii).*m .* ...
                          (Psh.alfvg(ii).*abs(x(a))).^(Psh.en(ii)-1) .* ...
                          (1+(Psh.alfvg(ii).*abs(x(a))).^Psh.en(ii)).^(-1-m);
            piece2  = Psh.alfvg2(ii).*Psh.en2(ii).*m2 .* ...
                          (Psh.alfvg2(ii).*abs(x(a))).^(Psh.en2(ii)-1) .* ...
                          (1+(Psh.alfvg2(ii).*abs(x(a))).^Psh.en2(ii)).^(-1-m2);
            C(a)    = df.*(Psh.fi(ii).*piece1 + (1-Psh.fi(ii)).*piece2);
        otherwise
            error('Bad setting of P.sh.ifr parameter at node %d!',ii)
    end
end
%% end
return