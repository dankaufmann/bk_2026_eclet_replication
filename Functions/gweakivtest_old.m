
function output = gweakivtest(y,Y,X,Z,varargin)

% Reference: Daniel Lewis and Karel Mertens, 
% A Robust Test for Weak Instruments with Multiple Endogenous Regressors 
% This version 14/10/2022

%--------------------------------------------------------------------------
% Model: 
% y       = Y*beta+X*By+u;

% Required Inputs
%----------------
% y  : Regressand (T x 1)
% Y  : Endogenous Regressors (T x N)
% X  : Exogenous Regressors (T x Nx)
% Z  : Instrumental Variables (T x K)

% Optional Inputs 
%----------------
% cov_type: HAR variance estimator, 'NW' or 'EHW', where 'NW': Newey West, 'EHW': Eicker-Huber-White, default is EHW  
% alfa:     confidence level, default is 0.05
% tau:      bias tolerance, default is 0.10
% points:   number of starting points for the optimization step, default is 1000 

% Outputs:
%---------
% nobs : Sample size
% beta_2SLS: Point estimates
% gmin_generalized: generalized test statistic
% gmin_generalized_critical_value: critical value for gmin_generalized
% gmin_generalized_simplified_critical_value: simplified critical value for gmin_generalized
% stock_yogo_test_statistic: test statistic for the Stock Yogo (2005) test
% stock_yogo_critical_value_nagar: critical value for stock_yogo_test_statistic based on the Nagar approximation (not the same as the value in the Stock Yogo tables which are based numerical integration)

%--------------------------------------------------------------------------

if nargin >4 
    if ~isempty(varargin{1})
        cov_type = varargin{1};
    else 
        cov_type = [];
    end
else 
    cov_type = [];
end

if nargin >5
    if ~isempty(varargin{2})
        alfa = varargin{2}; 
    else 
        alfa = 0.05; 
    end
else 
    alfa = 0.05; 
end

if nargin >6
    if ~isempty(varargin{3})
        tau = varargin{3}; 
    else 
        tau  = 0.10; 
    end
else 
    tau  = 0.10; 
end

if nargin >7
    if ~isempty(varargin{4})
        points = varargin{4}; 
    else 
        points = 1000; 
    end
else 
    points = 1000; 
end

 
% Force Dimensions 
if size(y,1)<size(y,2); y=y'; end
if size(Y,1)<size(Y,2); Y=Y'; end
if size(Z,1)<size(Z,2); Z=Z'; end
if size(X,1)<size(X,2); X=X'; end

% Drop Missing Observations
sel_sample = ~isnan(sum([y Y Z X],2));
y = y(sel_sample,:);
Y = Y(sel_sample,:);
Z = Z(sel_sample,:);

T = length(y);

% Add constant to X if absent
if ~isempty(X)
    X = X(sel_sample,:);
    X = [X(:,var(X)~=0) ones(T,1)];
else
    X = ones(T,1);
end

[~,Nx] = size(X);
[~,N]  = size(Y);
[~,K]  = size(Z);

% 
Zo = Z - X*(X\Z);
Zo = Zo-mean(Zo); 
Zo = Zo*(Zo'*Zo/T)^-0.5;
Yo = Y - X*(X\Y);
yo = y - X*(X\y);

PYo = Zo*(Zo\Yo);
Pyo = Zo*(Zo\yo);

betahat = PYo\Pyo;

v1 = yo-Pyo;
v2 = Yo-PYo;

ZV= repmat(Zo,1,N+1).*repelem([v1 v2] ,1,K);

if strcmp(cov_type,'NW')
    L=ceil(1.3*T^(1/2)); % Newey-West (1987) with truncation parameter recommended by Lazarus, Lewis, Stock, Watson (JBES 2018)
    for j=0:L
        if j>0
            acv(:,:,j+1)=ZV(j+1:end,:)'*ZV(1:end-j,:)/T + ZV(1:end-j,:)'*ZV(j+1:end,:)/T;
            w_l=1-j/(L);
            W=W+w_l*acv(:,:,j+1);
        else
            acv(:,:,j+1)=ZV'*ZV/T;
            W=acv(:,:,j+1);
        end
    end
else
    W=ZV'*ZV/T; % Eicker–Huber–White
end

W       = W*T/(T-K-Nx);

RNK     = kron(eye(N),reshape(eye(K),K*K,1));
W2      = W(K+1:end,K+1:end);
Phi     = RNK'*kron(W2,eye(K))*RNK;
gmin_generalized    = min(eig(Phi^-0.5*(PYo'*PYo)*Phi^-0.5));
[gmin_generalized_critical_value,gmin_generalized_critical_value_simplified,stock_yogo_critical_values_nagar] = gweakivtest_critical_values_old(W,K,alfa,tau,points);


Svv     = v2'*v2/(T-K-Nx);
gmin_stock_yogo = min(eig((Svv^-.5*Yo'*Zo*(Zo'*Zo)^-1*Zo'*Yo*Svv^-.5)/K));


output.nobs                                         = T;
output.beta_2SLS                                    = betahat';
output.gmin_generalized                             = gmin_generalized;
output.gmin_generalized_critical_value              = gmin_generalized_critical_value;
output.gmin_generalized_critical_value_simplified   = gmin_generalized_critical_value_simplified;
output.stock_yogo_test_statistic                    = gmin_stock_yogo;
output.stock_yogo_critical_value_nagar              = stock_yogo_critical_values_nagar;
