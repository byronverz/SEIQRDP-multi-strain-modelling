function Kappa = kappa_func(t,kappas)
kappa1 = kappas(1);
kappa2 = kappas(2);
Kappa = kappa1*exp(-kappa2*t);
end