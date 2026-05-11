function [wiv_cv,wiv_cv_simplified,wiv_cv_sy,Bmax] = gweakivtest_critical_valuesLRR1(W,K,Sig,varargin)

% Reference: Daniel Lewis and Karel Mertens,
% A Robust Test for Weak Instruments with Multiple Endogenous Regressors
% This version 22/09/2024

% Note : The 'OptStiefelGBB' function is by Z. Wen and W. Yin, A feasible method for optimization with orthogonality constraints

% Required:
% W: (N+1)*K x (N+1)*K HAR covariance matrix of score z*(w, v) in TRANSFORMED regression
% K:        number of instruments in TRANSFORMED regression
% Sig:      2 x 2 HAR covariance matrix of errors (w, v) in TRANSFORMED regression

% Optional:
% alfa:     confidence level, default is 0.05
% tau:      bias tolerance, default is 0.10
% points:   number of starting points for the optimization step, default is 1000
% target:   either 'beta' for entire vector or an integer j<=N corresponding to Y_j's position in beta
% Sigv:     N x N HAR covariance matrix of errors (w, v) in ORIGINAL regression; required if target is not 'beta'

if nargin >3
    if ~isempty(varargin{1})
        alfa = varargin{1};
    else
        alfa = 0.05;
    end
else
    alfa = 0.05;
end

if nargin >4
    if ~isempty(varargin{2})
        tau = varargin{2};
    else
        tau  = 0.10;
    end
else
    tau  = 0.10;
end

if nargin >5
    if ~isempty(varargin{3})
        points = varargin{3};
    else
        points = 1000;
    end
else
    points = 1000;
end

if nargin >6
    if ~isempty(varargin{4})
        target = varargin{4};
    else
        target = 'beta';
    end
else
    target = 'beta';
end

if nargin >7 && isnumeric(target)
    if ~isempty(varargin{5})
        Sigv = varargin{5};
    else
        fprintf('Error: Sigv required for single coefficient test')
        return
    end
elseif isnumeric(target)
    fprintf('Error: Sigv required for single coefficient test')
    return
end


N      = length(W)/K-1; % Number of endogenous regressors;


% Optimization Settings
opts.record = 0; %
opts.mxitr  = 1000;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1.e-7;
%options_fmincon = optimoptions('fmincon','Display','off');
%optimset('fmincon')
options_fmincon = optimset('Display', 'off');



% Construct some matrices
RNK     = kron(eye(N),reshape(eye(K),K*K,1));
RNN     = kron(eye(N),reshape(eye(N),N^2,1));
RNpK    = kron(eye(N+1),reshape(eye(K),K*K,1));
M1      = RNN'*(eye(N^3)+kron(spKgen(N,N),eye(N)));
M2      = (RNK*RNK'/(1+N)-eye(N*K^2));

W1      = W(1:K,1:K);
W2      = W(K+1:end,K+1:end);
W12     = W(1:K,K+1:end);
Phi=RNK'*kron(W2,eye(K))*RNK;
S       = kron((Phi/K)^-0.5,eye(K))*W2^0.5;
Sigma   = S*S';
Psibar     = kron(kron((Phi/K)^-0.5,eye(K))*[W12;W2]',eye(K))*RNpK;
Psi     = Psibar*Sig^-.5*norm(Phi^-.5*Sig(2:end,2:end)^.5);

X1      = kron(kron(speye(N),spKgen(K^2,N)),speye(N^2))*kron(reshape(speye(N),N^2,1),speye(K^2*N^2))*kron(kron(speye(K),spKgen(K,N)),speye(N))*(speye(N^2*K^2)+spKgen(N*K,N*K));
M2PsiM2 = M2*(Psi*Psi')*M2';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K < N
    fprintf('Error: not identified')
    return
end


if K>N+1
    Bmax(2) = min((2*(N+1)/K)^0.5*norm(M2*Psi),norm(Psi));
else
    Bmax(2) = max((2*(N+1)/K)^0.5*norm(M2*Psi),norm(Psi));
end

if K>N+1
    % Sharp upper bound
    for iter = 1:points
        [X,~]     = qr(randn(K,K));
        L0        = X(:,1:N)';
        [X, out1] = OptStiefelGBB(L0',@(x) objL0(x,M1,M2PsiM2,X1,N,K),opts);
        %   L0(:,:,iter) = X';
        Bmax_iters(iter) = sqrt(-out1.fval);
    end
    Bmax(1) = max(Bmax_iters);
else
    Bmax(1) = Bmax(2);
end

% Stock-Yogo under Nagar Approximation
if K>N+1
    Bmax(3) = (K-(1+N))/K;
else
    Bmax(3) = NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescale tau if necessary for median bias
if K==N
    fprintf('Model is just-identified, test is for median bias')
    if N==1
        tau=tau/0.455; % sharper bound for median bias
    end
end

% Rescale tau for single-coefficient test
if isnumeric(target)
    tau=tau*(sqrt(Sig(2,2))/sqrt(Sigv(target,target)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get critical value based on Imhof Approximation
for j = 1:3
    lmin = Bmax(j)/tau;
    if j <3 % Imhof approximation
        for n = 1:3
            k(n) = 2^(n-1)*factorial(n-1)*(norm(RNK'*(kron(Sigma^n,eye(K)))*RNK)+n*K*lmin*norm(Sigma)^(n-1));
        end
        ome = k(2)/k(3);
        nu  = 8*k(2)*ome^2;
        cc  = chi2inv(1-alfa,nu);
        % Check Kuhn-Tucker Conditions at the corner solution
        warning off
        fun_phiz  = @(z) ome*(1+(z-k(1))/(2*k(2)*ome)).^(nu/2-1).*exp(-nu/2*(1+(z-k(1))/(2*k(2)*ome)))*nu.^(nu/2-1)/(2^(nu/2-2))/gamma(nu/2);
        G1fun = @(q) -1/2*(q-2*nu*(nu-2)./q+nu) +3*nu/2*((log(q/2))-psi(0,nu/2));
        G2fun = @(q)  1/2*(q-nu*(nu-2)./q)  -nu*((log(q/2))-psi(0,nu/2));

        D1fun =  @(q) (1+(q-k(1))*2*ome)/(2*k(2)*ome).*(1+(q-k(1))/(2*k(2)*ome)).^(-1).*fun_phiz(q);
        D2fun =  @(q) fun_phiz(q)/k(2).*G1fun(nu+(q-k(1))*4*ome);
        D3fun =  @(q) G2fun(nu+(q-k(1))*4*ome)/k(3).*fun_phiz(q);

        ID1fun= @(xx) integral(D1fun,xx,Inf);
        ID2fun= @(xx) integral(D2fun,xx,Inf);
        ID3fun= @(xx) integral(D3fun,xx,Inf);

        kt_cond1 = ID1fun(((cc-nu)/4/ome+k(1)));
        kt_cond2 = ID2fun(((cc-nu)/4/ome+k(1)));
        kt_cond3 = ID3fun(((cc-nu)/4/ome+k(1)));
        kt_cond = (kt_cond1>=0)&&(kt_cond2>=0)&&(kt_cond3>=0);
        warning on
        if kt_cond~=1 % If Kuhn-Tucker Conditions fail, find cumulants that maximize the critical value at alfa numerically
            k_old = k;
            if N>1
                fun = @(x) -((chi2inv(1-alfa,8*x(2)*(x(2)/x(3))^2)-8*x(2)*(x(2)/x(3))^2)/4/(x(2)/x(3))+x(1));
                [k,fval] = fmincon(fun,k,eye(3),k,[],[],.01*ones(3,1),[],[],options_fmincon);
                ome  = k(2)/k(3);
                nu   = 8*k(2)*ome^2;
            else
                fun = @(x) -((chi2inv(1-alfa,8*x(1)*(x(1)/x(2))^2)-8*x(1)*(x(1)/x(2))^2)/4/(x(1)/x(2))+k(1));
                [knew,fval] = fmincon(fun,k(2:3),eye(2),k(2:3),[],[],.01*ones(2,1),[],[],options_fmincon);
                ome  = knew(1)/knew(2);
                nu   = 8*knew(1)*ome^2;
                k(2:3)=knew;
            end
            cc  = chi2inv(1-alfa,nu);
        end

        cv(j)     = ((cc-nu)/4/ome+k(1))/K;

    elseif j==3
        cv(j) = ncx2inv(1-alfa,K,K*lmin)/K;
    end
end

wiv_cv            = cv(1);
wiv_cv_simplified = cv(2);
wiv_cv_sy         = cv(3);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions

function [fval,gradient] = objL0(x,M1,M2PsiM2,X1,N,K) % Objective function and gradient
L0      = x';
vecL0   = reshape(L0,N*K,1);
QLL     = kron(kron(eye(N),L0),L0);
Mobj    = M1*QLL*M2PsiM2*QLL'*M1'/K;
Mobj    = 0.5*(Mobj+Mobj');
Mobj    = nearestSPD(Mobj);
[Qobj,Dobj] = eig(Mobj);
[~,ind]     = sort(diag(Dobj),'descend');
%Dobj        = Dobj(ind,ind);
Qobj        = Qobj(:,ind);
ev          = Qobj(:,1);
fval        = -ev'*Mobj*ev;
gradient    = 2*kron(ev'*M1*QLL*M2PsiM2,ev'*M1)*X1*kron(eye(N*K),vecL0);
gradient    =-reshape(gradient,N,K);
gradient    = gradient';
end

function [K] = spKgen(m,n)

%[m, n] = size(A);

I = reshape(1:m*n, [m, n]);
I = I';
I = I(:);
K = speye(m*n);
K = K(I,:);

end

function [X, out]= OptStiefelGBB(X, fun, opts, varargin)
%-------------------------------------------------------------------------
% curvilinear search algorithm for optimization on Stiefel manifold
%
%   min F(X), S.t., X'*X = I_k, where X \in R^{n,k}
%
%   H = [G, X]*[X -G]'
%   U = 0.5*tau*[G, X];    V = [X -G]
%   X(tau) = X - 2*U * inv( I + V'*U ) * V'*X
%
%   -------------------------------------
%   U = -[G,X];  V = [X -G];  VU = V'*U;
%   X(tau) = X - tau*U * inv( I + 0.5*tau*VU ) * V'*X
%
%
% Input:
%           X --- n by k matrix such that X'*X = I
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= OptStiefelGBB(X0, @fun, opts, data1, data2);
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 mxitr       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
%
% Output:
%           X --- solution
%         Out --- output information
%
% -------------------------------------
% For example, consider the eigenvalue problem F(X) = -0.5*Tr(X'*A*X);
%
% function demo
%
% function [F, G] = fun(X,  A)
%   G = -(A*X);
%   F = 0.5*sum(dot(G,X,1));
% end
%
% n = 1000; k = 6;
% A = randn(n); A = A'*A;
% opts.record = 0; %
% opts.mxitr  = 1000;
% opts.xtol = 1e-5;
% opts.gtol = 1e-5;
% opts.ftol = 1e-8;
%
% X0 = randn(n,k);    X0 = orth(X0);
% tic; [X, out]= OptStiefelGBB(X0, @fun, opts, A); tsolve = toc;
% out.fval = -2*out.fval; % convert the function value to the sum of eigenvalues
% fprintf('\nOptM: obj: %7.6e, itr: %d, nfe: %d, cpu: %f, norm(XT*X-I): %3.2e \n', ...
%             out.fval, out.itr, out.nfe, tsolve, norm(X'*X - eye(k), 'fro') );
%
% end
% -------------------------------------
%
% Reference:
%  Z. Wen and W. Yin
%  A feasible method for optimization with orthogonality constraints
%
% Author: Zaiwen Wen, Wotao Yin
%   Version 0.1 .... 2010/10
%   Version 0.5 .... 2013/10
%-------------------------------------------------------------------------


%% Size information
if isempty(X)
    error('input X is an empty matrix');
else
    [n, k] = size(X);
end

if nargin < 2; error('[X, out]= OptStiefelGBB(X0, @fun, opts)'); end
if nargin < 3; opts = [];   end

if ~isfield(opts, 'X0');        opts.X0 = [];  end
if ~isfield(opts, 'xtol');      opts.xtol = 1e-6; end
if ~isfield(opts, 'gtol');      opts.gtol = 1e-6; end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-12; end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'tau');       opts.tau  = 1e-3; end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.1; end
if ~isfield(opts, 'retr');      opts.retr = 0; end
if ~isfield(opts, 'gamma');     opts.gamma  = 0.85; end
if ~isfield(opts, 'STPEPS');    opts.STPEPS  = 1e-10; end
if ~isfield(opts, 'nt');        opts.nt  = 5; end
if ~isfield(opts, 'mxitr');     opts.mxitr  = 1000; end
if ~isfield(opts, 'record');    opts.record = 0; end
if ~isfield(opts, 'tiny');      opts.tiny = 1e-13; end

%-------------------------------------------------------------------------------
% copy parameters
xtol    = opts.xtol;
gtol    = opts.gtol;
ftol    = opts.ftol;
rhols   = opts.rhols;
STPEPS  = opts.STPEPS;
eta     = opts.eta;
gamma   = opts.gamma;
retr    = opts.retr;
record  = opts.record;
nt      = opts.nt;
crit    = ones(nt, 3);
tiny    = opts.tiny;
%-------------------------------------------------------------------------------

%% Initial function value and gradient
% prepare for iterations
[F,  G] = feval(fun, X , varargin{:});  out.nfe = 1;
GX = G'*X;

if retr == 1
    invH = true; if k < n/2; invH = false;  eye2k = eye(2*k); end
    if invH
        GXT = G*X';  H = 0.5*(GXT - GXT');  RX = H*X;
    else
        U =  [G, X];    V = [X, -G];       VU = V'*U;
        %U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
        %VX = VU(:,k+1:end); %VX = V'*X;
        VX = V'*X;
    end
end
dtX = G - X*GX;     nrmG  = norm(dtX, 'fro');

Q = 1; Cval = F;  tau = opts.tau;

%% Print iteration header if debug == 1
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '----------- Gradient Method with Line search ----------- \n');
    fprintf(fid, '%4s %8s %8s %10s %10s\n', 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff');
    %fprintf(fid, '%4d \t %3.2e \t %3.2e \t %5d \t %5d	\t %6d	\n', 0, 0, F, 0, 0, 0);
end

%% main iteration
for itr = 1 : opts.mxitr
    XP = X;     FP = F;   GP = G;   dtXP = dtX;
    % scale step size

    nls = 1; deriv = rhols*nrmG^2; %deriv
    while 1
        % calculate G, F,
        if retr == 1
            if invH
                [X, infX] = linsolve(eye(n) + tau*H, XP - tau*RX);
            else
                [aa, infR] = linsolve(eye2k + (0.5*tau)*VU, VX);
                X = XP - U*(tau*aa);
            end
        else
            [X, RR] = myQR(XP - tau*dtX, k);
        end

        if norm(X'*X - eye(k),'fro') > tiny; X = myQR(X,k); end

        [F,G] = feval(fun, X, varargin{:});
        out.nfe = out.nfe + 1;

        if F <= Cval - tau*deriv || nls >= 5
            break;
        end
        tau = eta*tau;          nls = nls+1;
    end

    GX = G'*X;
    if retr == 1
        if invH
            GXT = G*X';  H = 0.5*(GXT - GXT');  RX = H*X;
        else
            U =  [G, X];    V = [X, -G];       VU = V'*U;
            %U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
            %VX = VU(:,k+1:end); % VX = V'*X;
            VX = V'*X;
        end
    end
    dtX = G - X*GX;     nrmG  = norm(dtX, 'fro');
    S = X - XP;         XDiff = norm(S,'fro')/sqrt(n);
    tau = opts.tau;     FDiff = abs(FP-F)/(abs(FP)+1);

    %Y = G - GP;     SY = abs(iprod(S,Y));
    Y = dtX - dtXP;     SY = abs(iprod(S,Y));
    if mod(itr,2)==0; tau = (norm(S,'fro')^2)/SY;
    else tau  = SY/(norm(Y,'fro')^2); end
    tau = max(min(tau, 1e20), 1e-20);

    if (record >= 1)
        fprintf('%4d  %3.2e  %4.3e  %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr, tau, F, nrmG, XDiff, FDiff, nls);
        %fprintf('%4d  %3.2e  %4.3e  %3.2e  %3.2e (%3.2e, %3.2e)\n', ...
        %    itr, tau, F, nrmG, XDiff, alpha1, alpha2);
    end

    crit(itr,:) = [nrmG, XDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    %if (XDiff < xtol && nrmG < gtol ) || FDiff < ftol
    %if (XDiff < xtol || nrmG < gtol ) || FDiff < ftol
    %if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol
    %if ( XDiff < xtol || FDiff < ftol ) || nrmG < gtol
    %if any(mcrit < [gtol, xtol, ftol])
    if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol || all(mcrit(2:3) < 10*[xtol, ftol])
        out.msg = 'converge';
        break;
    end

    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + F)/Q;
end

if itr >= opts.mxitr
    out.msg = 'exceed max iteration';
end

out.feasi = norm(X'*X-eye(k),'fro');
if  out.feasi > 1e-13
    %X = MGramSchmidt(X);
    X = myQR(X,k);
    [F,G] = feval(fun, X, varargin{:});
    out.nfe = out.nfe + 1;
    out.feasi = norm(X'*X-eye(k),'fro');
end

out.nrmG = nrmG;
out.fval = F;
out.itr = itr;
end

function a = iprod(x,y)
%a = real(sum(sum(x.*y)));
a = real(sum(sum(conj(x).*y)));
end



function [Q, RR] = myQR(XX,k)
[Q, RR] = qr(XX, 0);
diagRR = sign(diag(RR)); ndr = diagRR < 0;
if nnz(ndr) > 0
    Q = Q*spdiags(diagRR,0,k,k);
    %Q(:,ndr) = Q(:,ndr)*(-1);
end
end

function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
    error('Exactly one argument must be provided.')
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
    error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
    % A was scalar and non-positive, so just return eps
    Ahat = eps;
    return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[U,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
    [R,p] = chol(Ahat);
    k = k + 1;
    if p ~= 0
        % Ahat failed the chol test. It must have been just a hair off,
        % due to floating point trash, so it is simplest now just to
        % tweak by adding a tiny multiple of an identity matrix.
        mineig = min(eig(Ahat));
        Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
    end
end


end
