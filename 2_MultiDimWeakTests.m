% 2_MultiDimWeakTests.m
% Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach
% by Marc Burri and Daniel Kaufmann
%
% Runs weak-instrument tests (Lewis and Mertens, 2025) for two-dimensional
% monetary policy shocks across all combinations of short- and medium/long-term
% interest rate specifications.  Tests are conducted for:
%   (a) HET-IV instruments (heteroskedasticity-based)
%   (b) HF-IV instruments (high-frequency surprises, non-recursive)
%   (c) HF-IV instruments (recursive, using both instruments jointly)
%
% Input:  ./Results/WeakData2Dim_*.mat     (written by 1_MultiDimShocks.R)
%         ./Results/WeakData2Dim_HF_*.mat  (written by 1_MultiDimShocks.R)
% Output: ./Results/WeakIVTest2Dim_Summary.txt     (HET-IV 2D test results)
%         ./Results/WeakIVTest1Dim_Summary.txt     (HET-IV 1D test results)
%         ./Results/WeakIVTest2Dim_HF_Summary.txt  (HF-IV non-recursive 2D)
%         ./Results/WeakIVTest2Dim_HFRec_Summary.txt (HF-IV recursive 2D)
%         ./Results/WeakIVTest1Dim_HF_Summary.txt  (HF-IV 1D test results)
%
% Run after: 1_MultiDimShocks.R
% Requires:  Octave with optim and statistics packages, or Matlab with Optimization Toolbox
%            gweakivtest_critical_values.m and gweakivtest.m by Lewis and Mertens (2025)
%
% Reference: Lewis, D.J. and Mertens, K. (2025). A robust test for weak instruments
%            with multiple endogenous regressors. Review of Economic Studies.

% Add path so that it finds the functions by Lewis and Mertens (2025)
addpath('./Functions/')

clear all; close all; clc;
pkg load optim
pkg load statistics

%--------------------------------------------------------------------------
% Settings
%--------------------------------------------------------------------------
% Significance level and bias tolerance for Lewis-Mertens (2025) test
% gmin_generalized statistic is compared against gmin_generalized_critical_value
cov_type   = 'NW';
alfa = 0.05;
tau  = 0.10;
points = 1000;
target1 = 1;
target2 = 2;
crit = 'abs';
code = 'old';   % old, new

% output = gweakivtest(y,Y,X,Z,cov_type,alfa,tau,points,target,crit)

%--------------------------------------------------------------------------
% Section 1: HET-IV weak instrument tests
% Tests heteroscedasticity-based instruments across all combinations of
% short-term (firstShock) and medium/long-term (secondShock) interest rates.
% For E=1: tests the first instrument alone (target shock).
% For E=2: tests both instruments jointly (target + path shock).
%--------------------------------------------------------------------------
firstShock  = {'FFR', 'IR3Mfed', 'IR6Mfed', 'IRSTfed'};
secondShock = {'IR2Yfed', 'IR3Yfed', 'IR5Yfed', 'IRMTfed'};

myResults = struct();
myResultsE1 = struct();
for i = 1:length(firstShock)



  for j = 1:length(secondShock)

    name1 = firstShock{i};
    name2 = secondShock{j};
    myFile = strcat(name1, "_", name2);

    load(strcat("./Results/WeakData2Dim_", myFile, ".mat"));

    %% Set up the data and perform the test for every shock (adapt if control variables)
    T = size(myTable.y1, 1);
    C = ones(T, 1);
    X = [C];

    %% Add control variables if any
    myNames = fieldnames(myTable);
    startsWithO = strncmp(myNames, 'o', 1);
    if(sum(startsWithO)>0)
      myNames = myNames(startsWithO);

      for o = 1:numel(myNames)

        varname = myNames{o};
        X = [X, myTable.(varname)];

      end
    end

    if j == 1

      % Variables for e = 2, without exogenous regressors
      Z1 = [myTable.Z1];
      Y1 = [myTable.y1];
      y1 = myTable.y3;

      disp('Results e = 1');
      if(strcmp(code, 'old'))
        resE1 = gweakivtest_old(y1,Y1,X,Z1,cov_type,alfa,tau,points)
      end
      if(strcmp(code, 'new'))
        resE1 = gweakivtest(y1,Y1,X,Z1,cov_type,alfa,tau,points,target1,crit)
      end

      % Save results for later use
      myResultsE1.(myFile) = [resE1.beta_2SLS, resE1.gmin_generalized, resE1.gmin_generalized_critical_value];

    endif

    % Variables for e = 2, without exogenous regressors
    Z2 = [myTable.Z1, myTable.Z2];
    Y2 = [myTable.y1, myTable.y2];
    y2 = myTable.y3;

    disp(myFile)
    disp('Results e = 2');
    if(strcmp(code, 'old'))
      resE2 = gweakivtest_old(y2,Y2,X,Z2,cov_type,alfa,tau,points)
    end
    if(strcmp(code, 'new'))
      resE2 = gweakivtest(y2,Y2,X,Z2,cov_type,alfa,tau,points,target2,crit)
    end

    % Save results for later use
    myResults.(myFile) = [resE2.beta_2SLS, resE2.gmin_generalized, resE2.gmin_generalized_critical_value];
  end
end

%Save results to text file
outfile = fopen("./Results/WeakIVTest2Dim_Summary.txt", "w");
fprintf(outfile, "\n %s \t %s \t %s \t %s \t %s", "Spec.", "1st Beta", "F-stat", "Crit. val.");

for CellStr = fieldnames(myResults).'
  Name = CellStr{1};
  fprintf(outfile, "\n %s \t %.2f \t %.2f %.2f %.2f", Name, myResults.(Name));

end

 fflush(outfile);
 fclose(outfile);


%Save results to text file
outfile = fopen("./Results/WeakIVTest1Dim_Summary.txt", "w");
fprintf(outfile, "\n %s \t %s \t %s \t %s \t %s", "Spec.", "1st Beta", "F-stat", "Crit. val.");

for CellStr = fieldnames(myResultsE1).'
  Name = CellStr{1};
  Name2 =  strsplit(CellStr{1}, '_'){1};
  fprintf(outfile, "\n %s \t %.2f \t %.2f %.2f %.2f", Name2, myResultsE1.(Name));

end

 fflush(outfile);
 fclose(outfile);



%--------------------------------------------------------------------------
% Section 2: HF-IV weak instrument tests
% Tests high-frequency surprise instruments (Swanson, 2021).
% Non-recursive: each instrument tested separately (one instrument per equation).
% Recursive: both instruments tested jointly (Z1 and Z2 together for second shock).
%--------------------------------------------------------------------------
myResults = struct();
myResultsE1 = struct();
myResultsRec = struct();

for i = 1:length(firstShock)

  for j = 1:length(secondShock)

    name1 = firstShock{i};
    name2 = secondShock{j};
    myFile = strcat(name1, "_", name2);

    load(strcat("./Results/WeakData2Dim_HF_", myFile, ".mat"));

    %% Set up the data and perform the test for every shock (adapt if control variables)
    T = size(myTable.y1, 1);
    C = ones(T, 1);
    X = [C];

    %% Add control variables if any
    myNames = fieldnames(myTable);
    startsWithO = strncmp(myNames, 'o', 1);
    if(sum(startsWithO)>0)
      myNames = myNames(startsWithO);

      for o = 1:numel(myNames)

        varname = myNames{o};
        X = [X, myTable.(varname)];

      end
    end

    if j == 1

      % Variables for e = 1, without exogenous regressors
      Z1 = [myTable.Z1];
      Y1 = [myTable.y1];
      y1 = myTable.y3;

      disp('Results e = 1');
      if(strcmp(code, 'old'))
        resE1 = gweakivtest_old(y1,Y1,X,Z1,cov_type,alfa,tau,points)
      end
      if(strcmp(code, 'new'))
        resE1 = gweakivtest(y1,Y1,X,Z1,cov_type,alfa,tau,points,target1,crit)
      end

      % Save results for later use
      myResultsE1.(myFile) = [resE1.beta_2SLS, resE1.gmin_generalized, resE1.gmin_generalized_critical_value];

    endif

    if i == 1
      % Variables for e = 2, without exogenous regressors
      Z2 = [myTable.Z2];
      Y2 = [myTable.y2];
      y2 = myTable.y3;

      disp(myFile)
      disp('Results e = 2');

      if(strcmp(code, 'old'))
        resE2 = gweakivtest_old(y2,Y2,X,Z2,cov_type,alfa,tau,points)
      end
      if(strcmp(code, 'new'))
        resE2 = gweakivtest(y2,Y2,X,Z2,cov_type,alfa,tau,points,target1,crit)
      end

      % Save results for later use
      myResults.(myFile) = [resE2.beta_2SLS, resE2.gmin_generalized, resE2.gmin_generalized_critical_value];
    endif

      % Variables for e = 2, without exogenous regressors (recursive)
    Z2 = [myTable.Z1, myTable.Z2];
    Y2 = [myTable.y1, myTable.y2];
    y2 = myTable.y3;

    disp(myFile)
    disp('Results e = 2');
    if(strcmp(code, 'old'))
     resE2 = gweakivtest_old(y2,Y2,X,Z2,cov_type,alfa,tau,points)
    end
    if(strcmp(code, 'new'))
      resE2 = gweakivtest(y2,Y2,X,Z2,cov_type,alfa,tau,points,target2,crit)
    end


    % Save results for later use
    myResultsRec.(myFile) = [resE2.beta_2SLS, resE2.gmin_generalized, resE2.gmin_generalized_critical_value];

  end
end

%Save results to text file
outfile = fopen("./Results/WeakIVTest2Dim_HF_Summary.txt", "w");
fprintf(outfile, "\n %s \t %s \t %s \t %s \t %s", "Spec.", "1st Beta", "F-stat", "Crit. val.");

for CellStr = fieldnames(myResults).'
  Name = CellStr{1};
  Name2 =  strsplit(CellStr{1}, '_'){2};
  fprintf(outfile, "\n %s \t %.2f \t %.2f %.2f %.2f", Name2, myResults.(Name));

end

 fflush(outfile);
 fclose(outfile);


%Save results to text file
outfile = fopen("./Results/WeakIVTest2Dim_HFRec_Summary.txt", "w");
fprintf(outfile, "\n %s \t %s \t %s \t %s \t %s", "Spec.", "1st Beta", "F-stat", "Crit. val.");

for CellStr = fieldnames(myResultsRec).'
  Name = CellStr{1};
  fprintf(outfile, "\n %s \t %.2f \t %.2f %.2f %.2f", Name, myResultsRec.(Name));

end

 fflush(outfile);
 fclose(outfile);



%Save results to text file
outfile = fopen("./Results/WeakIVTest1Dim_HF_Summary.txt", "w");
fprintf(outfile, "\n %s \t %s \t %s \t %s \t %s", "Spec.", "1st Beta", "F-stat", "Crit. val.");

for CellStr = fieldnames(myResultsE1).'
  Name = CellStr{1};
  Name2 =  strsplit(CellStr{1}, '_'){1};
  fprintf(outfile, "\n %s \t %.2f \t %.2f %.2f %.2f", Name2, myResultsE1.(Name));

end

 fflush(outfile);
 fclose(outfile);




