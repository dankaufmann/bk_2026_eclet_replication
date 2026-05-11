
function output = gweakivtest_LRR1(y,Y,X,Z,ind,varargin)

% Reference: Daniel Lewis and Karel Mertens,
% A Robust Test for Weak Instruments with Multiple Endogenous Regressors
% First version 9/22/2024
% This version 26/06/2025

%--------------------------------------------------------------------------
% Model:
% y       = Y*beta+X*By+u;

% Required Inputs
%----------------
% y  : Regressand (T x 1)
% Y  : Endogenous Regressors (T x N)
% X  : Exogenous Regressors (T x Nx)
% Z  : Instrumental Variables (T x K)
% ind: Integer less than or equal to N specifying retained regressor in
% transformation

% Optional Inputs
%----------------
% cov_type: HAR variance estimator, 'NW' or 'EHW', where 'NW': Newey West, 'EHW': Eicker-Huber-White, default is EHW
% alfa:     confidence level, default is 0.05
% tau:      bias tolerance, default is 0.10
% points:   number of starting points for the optimization step, default is 1000
% target:   either 'beta' or equal to the integer input 'ind'. Default is
% beta.

% Outputs:
%---------
% nobs : Sample size
% beta_2SLS: Point estimates
% gmin_generalized: generalized test statistic
% gmin_generalized_critical_value: critical value for gmin_generalized
% gmin_generalized_simplified_critical_value: simplified critical value for gmin_generalized
% stock_yogo_test_statistic: test statistic for the Stock Yogo (2005) test,
% as in Sanderson and Windmeijer (2018)
% stock_yogo_critical_value_nagar: critical value for stock_yogo_test_statistic based on the Nagar approximation (not the same as the value in the Stock Yogo tables which are based numerical integration)

%--------------------------------------------------------------------------

if nargin >5
    if ~isempty(varargin{1})
        cov_type = varargin{1};
    else
        cov_type = [];
    end
else
    cov_type = [];
end

if nargin >6
    if ~isempty(varargin{2})
        alfa = varargin{2};
    else
        alfa = 0.05;
    end
else
    alfa = 0.05;
end

if nargin >7
    if ~isempty(varargin{3})
        tau = varargin{3};
    else
        tau  = 0.10;
    end
else
    tau  = 0.10;
end

if nargin >8
    if ~isempty(varargin{4})
        points = varargin{4};
    else
        points = 1000;
    end
else
    points = 1000;
end
if nargin >9
    if ~isempty(varargin{5})
        target = varargin{5};
        if isnumeric(target) && target~=ind
            fprintf('Error: Target coefficient ''target'' must be the same as retained regressor, ''ind''')
            return
        end
    else
        target = 'beta';
    end
else
    target = 'beta';
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
nind=1:N; nind(ind)=[];

%
Zo = Z - X*(X\Z);
Zo = Zo-mean(Zo);
Zo = Zo*(Zo'*Zo/T)^-0.5;
Yo = Y - X*(X\Y);
yo = y - X*(X\y);

PYo = Zo*(Zo\Yo);
Pyo = Zo*(Zo\yo);
betahat = PYo\Pyo;
what = yo-Pyo;
vhat = Yo-PYo;
ehat=[what vhat];

dtilde=PYo(:,nind)\Yo(:,ind);
Ystar=Yo(:,ind)-PYo(:,nind)*dtilde;
ystar=yo(:,1)-PYo(:,nind)*(PYo(:,nind)\yo(:,1));
Z2=Zo(:,N:end)-PYo(:,nind)*(PYo(:,nind)\Zo(:,N:end));
shat=cov(vhat,1);
dbar=[-dtilde(1:ind-1); 1; -dtilde(ind:end)];
denom=dbar'*shat*dbar;
gmin_stock_yogo=Ystar'*Z2*(Z2'*Z2)^-1*Z2'*Ystar/((K-N+1)*denom);

Yhs=Z2*(Z2\Ystar);
Z2=Z2*(Z2'*Z2/T)^-.5;
estar=[what vhat(:,ind)-vhat(:,nind)*dtilde];
Kstar=size(Z2,2);
ZV= [repmat(Z2,1,2).*repelem(estar,1,Kstar)];

if strcmp(cov_type,'NW')
    L=ceil(1.3*T^(1/2)); % Newey-West (1987) with truncation parameter recommended by Lazarus, Lewis, Stock, Watson (JBES 2018)
    for j=0:L
        if j>0
            acv(:,:,j+1)=ZV(j+1:end,:)'*ZV(1:end-j,:)/T + ZV(1:end-j,:)'*ZV(j+1:end,:)/T;
            w_l=1-j/(L);
            Wstar=Wstar+w_l*acv(:,:,j+1);
            acvS(:,:,j+1)=estar(j+1:end,:)'*estar(1:end-j,:)/T + estar(1:end-j,:)'*estar(j+1:end,:)/T;            
            Sstar=Sstar+w_l*acvS(:,:,j+1);
            acvs(:,:,j+1)=vhat(j+1:end,:)'*vhat(1:end-j,:)/T + vhat(1:end-j,:)'*vhat(j+1:end,:)/T;
            S_full=S_full+w_l*acvs(:,:,j+1);
        else
            acv(:,:,j+1)=ZV'*ZV/T;
            Wstar=acv(:,:,j+1);
            acvS(:,:,j+1)=estar'*estar/T;
            Sstar=acvS(:,:,j+1);
            acvs(:,:,j+1)=vhat'*vhat/T;
            S_full=acvs(:,:,j+1);
        end
    end
else
    Wstar=ZV'*ZV/T; % Eicker–Huber–White \
    Sstar=estar'*estar/T;
    S_full=vhat'*vhat/T;
end

Wstar       = Wstar*T/(T-K-Nx);
Sstar   = Sstar*T/(T-K-Nx);
Phi=trace(Wstar(Kstar+1:end,Kstar+1:end));
gmin_generalized=Yhs'*Yhs/Phi;
[gmin_generalized_critical_value,gmin_generalized_critical_value_simplified,stock_yogo_critical_values_nagar] = gweakivtest_critical_valuesLRR1(Wstar,Kstar,Sstar,alfa,tau,points,target,S_full);

output.nobs                                         = T;
output.beta_2SLS                                    = betahat';
if strcmp(target,'beta')
    output.target='beta vector';
else
    output.target=strcat('beta_',num2str(target));
end
output.gmin_generalized                             = gmin_generalized;
output.criterion                                    = 'abs';
output.gmin_generalized_critical_value              = gmin_generalized_critical_value;
output.gmin_generalized_critical_value_simplified   = gmin_generalized_critical_value_simplified;
output.sanderson_windmeijer_test_statistic          = gmin_stock_yogo;
output.sanderson_windmeijer_critical_value_nagar    = stock_yogo_critical_values_nagar;

