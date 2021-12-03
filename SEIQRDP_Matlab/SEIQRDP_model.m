function X = SEIQRDP_model(t, x0, params,wave)
X = zeros(size(x0));
S = x0(1);
E = x0(2);
I = x0(3);
Q = x0(4);
R = x0(5);
D = x0(6);
P = x0(7);

alpha = params(1);
beta1 = params(2);
beta2 = params(3);
beta3 = params(4);
gamma = params(5);
delta = params(6);
lambda1 = params(7);
lambda2 = params(8);
lambda3 = params(9);
kappa1 = params(10);
kappa2 = params(11);

beta = beta_func(t,beta1,beta2,beta3,wave);
% lambda = params(5); % 0.039 - 0.039
lambda = lambda_func(t,[lambda1,lambda2,lambda3]);

% kappa = params(6); % 0.002
kappa = kappa_func(t,[kappa1,kappa2]);
N = 57.78e6;

X(1) = -alpha*S - (beta/N)*S*I;      % Susceptible
X(2) = -gamma*E + (beta/N)*S*I;      % Exposed
X(3) = gamma*E - delta*I;            % Infectious
X(4) = delta*I - (lambda + kappa)*Q; % Quarantined
X(5) = lambda*Q;                     % Recovered
X(6) = kappa*Q;                      % Deceased
X(7) = alpha*S;                      % Insusceptible
end