function teta = fnteta( x, Psh, i )
% teta = fnteta( x, Psh, i )
% 
% NOTE
%   Add hysteresis!!
% 
% DESCRIPTION
%   Computes the moisture fraction at node(s) i according to the passed
%   pressure head in x (which can also account for other forms of
%   potential...).
%   A similar function is watcon in functions.for, SWAP-32.
% 
% INPUTs
%   x:      Matrix potential (e.g. pressure head "h"), at i.
%   Psh:    The set of hydraulic characteristics of the soil grid at all
%           nodes.
%   i:      Two different usages of the function are possible according to
%           the value assigned to i:
%               *one value  --> Moisture at current i node of the soil
%                               grid.
%               *multiple   --> Moisture at all nodes passed in i.
% 
% OUTPUTs
%   teta:   Moisture fraction at node(s) i.

%% main
teta            = NaN(length(i),1);
if numel(x)~=1% at top/bottom nodes calls I don't need this adjustment:
    x           = x(i);
end
for a = 1:length(i)
    ii          = i(a);
    df          = Psh.tetas(ii)-Psh.tetar(ii);
    m           = 1-1./Psh.en(ii);
    m2          = 1-1./Psh.en2(ii);
    
    % adjust moistures according to the value of x:
    if x(a)>0
        teta(a) = Psh.tetas(ii);
        continue
    end
    
    switch Psh.ifr(ii)
        case 1
            % SWAP-32, Eq. 2.4, page 27:
            teta(a) = Psh.tetar(ii) + df./(1+(Psh.alfvg(ii).*abs(x(a))).^Psh.en(ii)).^m;%#ERROR: abs on alfa!
        case 2
            teta(a) = Psh.tetar(ii) + df .* ( Psh.fi(ii) .* ...
                      ( (1+Psh.alfrs(ii).*abs(x(a))).*exp(-Psh.alfrs(ii).*abs(x(a))) ) + ...
                      (1-Psh.fi(ii)).*(1+(Psh.alfvg2(ii).*abs(x(a))).^Psh.en2(ii)).^(-m2)    );
        case 3
            sef1    = (1+(Psh.alfvg (ii).*abs(x(a))).^Psh.en(ii)).^-m;
            sef2    = (1+(Psh.alfvg2(ii).*abs(x(a))).^Psh.en2(ii)).^-m2;
            sef     = Psh.fi(ii).*sef1 + (1-Psh.fi(ii)).*sef2;
            teta(a) = Psh.tetar(ii) + df.*sef;
        case 4% modified Van Genuchten
            if x(a) < Psh.h_enpr(ii)
                %#modification: abs outside alfa!
                Sc      = 1./(1+abs(Psh.alfvg(ii).*Psh.h_enpr(ii)).^Psh.en(ii)).^m;
                teta(a) = Psh.tetar(ii) + ( df./(1+(Psh.alfvg(ii).*abs(x(a))).^Psh.en(ii)).^m  ) ./ Sc;
            elseif x(a) >= Psh.h_enpr(ii)
                teta(a) = Psh.tetas(ii);
            end
        otherwise
            error('Bad setting of P.sh.ifr parameter at node %d!',ii)
    end
end
%% end
return