function ics_fin = mln_bisection(teta,dap,S1,Kf1,Kf2,C2phys,SsSk,Wdt)

f_bis_sol       = @(yps) ...
                    teta*yps        - ...
                    teta*C2phys     + ...
                    dap*Kf1*yps^Kf2 - ...
                    dap*S1          - ...
                    SsSk*Wdt;
a_bs            = -100;
b_bs            = 100;
espon           = -30;
eps_bs          = 10^espon;
fa              = feval(f_bis_sol,a_bs);
fb              = feval(f_bis_sol,b_bs);
if fa*fb>0
    error('Intervallo iniziale non accettabile')
end
Nit_sol         = (log(b_bs-a_bs)-log(eps_bs))/log(2);
for u=3:Nit_sol+2
    ics         = (a_bs+b_bs)/2;
    fics        = feval(f_bis_sol,ics);
    if fa*fics<=0
        b_bs    = ics;
    else
        a_bs    = ics;
        fa      = fics;
    end
end
ics_fin         = (a_bs+b_bs)/2;

return