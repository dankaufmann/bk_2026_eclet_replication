clear all; close all;

 rng(0)

% Octave packages
pkg load optim
pkg load statistics

 % True Model Parameters
 %%%%%%%%%%%%%%%%%%%%%%%
 K  = 2; % Number of instruments
 N  = 2; % Number of endogenous variables
 Nx = 2; % Number of Exogenous Regressors

 beta = ones(N,1); % Structural Parameters

 DL      = 0.08*diag(ones(N,1));
 QL      = eye(N);
 [X,~]   = qr(randn(K,K));
 L0      = X(:,1:N)';
 Pi      = sqrt(K)*QL*DL^0.5*L0;
 Pi      = Pi'; % First Stage Parameters

 By       = randn(Nx+1,1); % Parameters on Exogenous Regressors
 BY       = randn(Nx+1,N);

 % Generate Sample
 %%%%%%%%%%%%%%%%%
 T       = 200;        % Sample Size
 u       = randn(T,1); % True Structural Equation Errors
 v       = randn(T,N); % True First-Stage Errors

 CONST   = ones(T,1);
 X       = [randn(T,Nx) CONST];
 Z       = randn(T,K);
 Y       = Z*Pi+X*BY+v;
 y       = Y*beta+X*By+u;


% Test for bias:

output = gweakivtest(y,Y,X,Z,'NW')

% With Optional Inputs
% alfa = 0.05;
% tau  = 0.10;
% points = 1000;
% target = 'beta'
% output = gweakivtest(y,Y,X,Z,'NW',alfa,tau,points,target)

% Test for bias in first coefficient only:

output_1 = gweakivtest(y,Y,X,Z,'NW',[],[],[],1)

% LRR1 tests (projecting out Y2hat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_LRR1 = gweakivtest_LRR1(y,Y,X,Z,1,'NW')

% Test for bias in first coefficient only:

output_LRR1_1 = gweakivtest_LRR1(y,Y,X,Z,1,'NW',[],[],[],1)

