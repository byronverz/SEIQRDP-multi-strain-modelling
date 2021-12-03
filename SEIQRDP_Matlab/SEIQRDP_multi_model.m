function X = SEIQRDP_multi_model(t, x0, params,wave)
X = zeros(size(x0));
S = x0(1);
E_s1 = x0(2);
E_s2 = x0(3);
E_s3 = x0(4);
I_s1 = x0(5);
I_s2 = x0(6);
I_s3 = x0(7);
Q_s1 = x0(8);
Q_s2 = x0(9);
Q_s3 = x0(10);
R = x0(11);
D = x0(12);
P = x0(13);

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

gamma_s1 = gamma;
gamma_s2 = gamma;
gamma_s3 = gamma;

beta = beta_func(t,beta1,beta2,beta3,wave);
beta_s1 = beta;
if wave == 1
    t_true = t;
elseif wave == 2
    t_true = t + days(datetime([2020,10,8])-datetime([2020,4,18]));
elseif wave == 3
    t_true = t + days(datetime([2021,3,30])-datetime([2020,4,18]));
end
if t_true>160
    beta_s2 = 1.08*beta_func(t-160,beta1,beta2,beta3,2);
    beta = beta_s2;
else
    beta_s2 = 0.00*beta;
end
if t_true>346
    beta_s3 = beta_func(t-346,beta1,beta2,beta3,3);
    beta = beta_s3;
else
    beta_s3 = 0.00*beta;
end
% lambda = params(5); % 0.039 - 0.039
lambda = lambda_func(t,[lambda1,lambda2,lambda3]);

% kappa = params(6); % 0.002
kappa = kappa_func(t,[kappa1,kappa2]);
N = 57.78e6;

X(1) = -alpha*S - (beta/N)*S*(I_s1+I_s2+I_s3);      % Susceptible
X(2) = -gamma_s1*E_s1 + (beta_s1/N)*S*I_s1;      % Exposed S1
X(3) = -gamma_s2*E_s2 + (beta_s2/N)*S*I_s2;      % Exposed S2
X(4) = -gamma_s3*E_s3 + (beta_s3/N)*S*I_s3;      % Exposed S3
X(5) = gamma_s1*E_s1 - delta*I_s1;            % Infectious S1
X(6) = gamma_s2*E_s2 - delta*I_s2;            % Infectious S2
X(7) = gamma_s3*E_s3 - delta*I_s3;            % Infectious S3
X(8) = delta*I_s1 - (lambda + kappa)*Q_s1; % Quarantined S1
X(9) = delta*I_s2 - (lambda + kappa)*Q_s2; % Quarantined S1
X(10) = delta*I_s3 - (lambda + kappa)*Q_s3; % Quarantined S1
X(11) = lambda*(Q_s1+Q_s2+Q_s3);                     % Recovered
X(12) = kappa*(Q_s1+Q_s2+Q_s3);                      % Deceased
X(13) = alpha*S;                      % Insusceptible
end