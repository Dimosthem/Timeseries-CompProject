%% Timeseries analysis project

%Bounarelis Dimosthenis
%Giachoudis Christos


%% Import data and plot

%******Import a .dat file******
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["e03", "Var2"];
opts.SelectedVariableNames = "e03";
opts.VariableTypes = ["double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Var2", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var2", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "e03", "TrimNonNumeric", true);
opts = setvaropts(opts, "e03", "ThousandsSeparator", ",");

% Import the data
tbl = readtable("C:\Users\Christos\chris\auth\Timeseries\project\normalEEG.dat", opts);

% Convert to output type
normEEG = tbl.e03;

% Clear temporary variables
clear opts tbl

% Plot the data
x1 = 1:length(normEEG);
y1 = normEEG;
figure("Name", "Timeserie", "NumberTitle", "off");
plot(x1,y1, "color", "#A2142F");

% or

% plot the data in parts
figure("Name", "History", "NumberTitle", "off");
tittxt = 'normEEG';
nseg = 3;
type = 'b';
tau_s = 1;
pser(normEEG, tittxt, nseg, type, tau_s);

%% Linear analysis

% _____Stationarity_____
% (no need to make the normal EEG stationary, as it already is)
% maorder = 100;
% stationary_normEEG = movingaveragesmooth(normEEG, maorder);
% Plot the stationary data
% xS = 1:length(stationary_normEEG);
% yS = stationary_normEEG;
% figure("Name", "Stationary Timeserie of normal EEG", "NumberTitle", "off");
% plot(xS,yS, "color", "#A2142F");


% _____Mean value_____and_____variance_____
normEEG_mean = mean(normEEG);
formatSpec = '******The mean value of the normal EEG is %4.4f\n';
fprintf(formatSpec,normEEG_mean);
normEEG_var = var(normEEG);
formatSpec = '******The variance of the normal EEG is %4.4f\n';
fprintf(formatSpec,normEEG_var);


% _____Box-Cox method for stabilizing standard deviation/variance_____

% !!!!!!if you want to make the best transformation you can
% [transdat,lambda] = boxcox(normEEG);
% if lambda == 0
%     formatSpec = 'We use the following: transdat = log(normEEG)\n';
%     fprintf(formatSpec);
% else
%     formatSpec = 'We use the following: transdat = (normEEG^lambda - 1)/lambda \n';
%     fprintf(formatSpec);
% end

% !!!!!!if you just want a logarithmic transformation
% lambda = 0;
% transdat = boxcox(lambda,normEEG);

% Plot the transformed data
% xBC = 1:length(normEEG);
% yBC = transdat;
% figure("Name", "BOX-COX-transform", "NumberTitle", "off");
% plot(xBC, yBC, "color", "#EDB120");


%_____Periodicity_____
% per = input('Give the periodicity order of timeserie >');
% non_per = seasonalcomponents(normEEG, per);
% Plot the non-periodical data
% xnp = 1:length(normEEG);
% ynp = transdat;
% figure("Name", "Non-periodical data", "NumberTitle", "off");
% plot(xnp, ynp, "color", "#EDB120");

%_____White Noise or more information ?_____

% autocorrelation
%method 1
lags = input('Give how many lags you want to compute for autocorrelation and partial autocorrelation>');
figure("Name", "Autocorrelation", "NumberTitle", "off");
autocorr(normEEG, 'NumLags', lags);
%method 2
% tmax = lags;
% tittxt = 'Autocorrelation for normEEG';
% type = 'cd';
% [acM] = autocorrelation(normEEG, tmax, tittxt, type);

% Portmanteau test to be sure
thesize = lags;
alpha = 0.05;
tittxt = 'normal EEG';
[hV, pV, QV, xautV] = portmanteauLB(normEEG,thesize,alpha,tittxt);
% formatSpec = 'The h value of the Portmanteau test for the normal EEG is %4.4f\n';
% fprintf(formatSpec, hV);

% Partial autocorrelation
figure("Name", "Partial Autocorrelation", "NumberTitle", "off");
parcorr(normEEG, 'NumLags', lags);


%_____Linear Models Simulations_____

n = length(normEEG);

% ***AR model***
% q = 0;
% p = input('Give the order p of the AR model >');
% Tmax = input('How many steps ahead you want the function to compute the fitting error ? >');
% [ARnrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(normEEG, p, q, Tmax);
% fprintf('===== AR model ===== \n');
% fprintf('Estimated coefficients of phi(B):\n');
% disp(phiallV');
% fprintf('Estimated coefficients of theta(B):\n');
% disp(thetaallV');
% fprintf('SD of noise: %f \n',SDz);
% fprintf('AIC: %f \n',aicS);
% fprintf('FPE: %f \n',fpeS);
% fprintf('\t T \t\t NRMSE \n');
% disp([[1:Tmax]' ARnrmseV]);
% figure('Name', 'Normalized Root Mean Square Error for AR', 'NumberTitle', 'off');
% plot([1:Tmax]',ARnrmseV,'.-k')
% hold on
% plot([1 Tmax],[1 1],'y')
% xlabel('T')
% ylabel('NRMSE')
% title(sprintf('AR(%d,%d), fitting error',p,q))
% % Further prediction error analysis
% nlast = 0.2 * n;
% figure('Name', 'Normalized Root Mean Square Error for AR v2', 'NumberTitle', 'off');
% [ARnrmseV,ARpreM] = predictARMAnrmse(normEEG,p,q,Tmax,nlast,'AR');
% figure("Name", "Comparison for AR", "NumberTitle", "off");
% plot([n-nlast+1:n]',normEEG(n-nlast+1:n),'.-')
% hold on
% plot([n-nlast+1:n]',ARpreM(:,1),'.-r')
% if Tmax>1
%     plot([n-nlast+1:n]',ARpreM(:,2),'.-c')
% 	if Tmax>2
%         plot([n-nlast+1:n]',ARpreM(:,3),'.-k')
%     end
% end
% switch Tmax
%     case 1
%         legend('true','T=1','Location','Best')
%     case 2
%         legend('true','T=1','T=2','Location','Best')
%     otherwise
%         legend('true','T=1','T=2','T=3','Location','Best')
% end

% ***MA model***
% p = 0;
% q = input('Give the order p of the MA model >');
% Tmax = input('How many steps ahead you want the function to compute the fitting error ? >');
% [MAnrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(normEEG, p, q, Tmax);
% fprintf('===== MA model ===== \n');
% fprintf('Estimated coefficients of phi(B):\n');
% disp(phiallV');
% fprintf('Estimated coefficients of theta(B):\n');
% disp(thetaallV');
% fprintf('SD of noise: %f \n',SDz);
% fprintf('AIC: %f \n',aicS);
% fprintf('FPE: %f \n',fpeS);
% fprintf('\t T \t\t NRMSE \n');
% disp([[1:Tmax]' MAnrmseV]);
% figure('Name', 'Normalized Root Mean Square Error for MA', 'NumberTitle', 'off');
% plot([1:Tmax]',MAnrmseV,'.-k')
% hold on
% plot([1 Tmax],[1 1],'y')
% xlabel('T')
% ylabel('NRMSE')
% title(sprintf('MA(%d,%d), fitting error',p,q))
% % Further prediction error analysis
% nlast = 0.2 * length(normEEG);
% figure('Name', 'Normalized Root Mean Square Error for MA v2', 'NumberTitle', 'off');
% [MAnrmseV,MApreM] = predictARMAnrmse(normEEG,p,q,Tmax,nlast,'MA');
% figure("Name", "Comparison for MA", "NumberTitle", "off");
% plot([n-nlast+1:n]',normEEG(n-nlast+1:n),'.-')
% hold on
% plot([n-nlast+1:n]',MApreM(:,1),'.-r')
% if Tmax>1
%     plot([n-nlast+1:n]',MApreM(:,2),'.-c')
% 	if Tmax>2
%         plot([n-nlast+1:n]',MApreM(:,3),'.-k')
%     end
% end
% switch Tmax
%     case 1
%         legend('true','T=1','Location','Best')
%     case 2
%         legend('true','T=1','T=2','Location','Best')
%     otherwise
%         legend('true','T=1','T=2','T=3','Location','Best')
% end


% ***ARMA model***
p = input('Give the order p of the AR part of ARMA >');
q = input('Give the order p of the MA part of ARMA >');
Tmax = input('How many steps ahead you want the function to compute the fitting error ? >');
[ARMAnrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(normEEG, p, q, Tmax);
fprintf('===== ARMA model ===== \n');
fprintf('Estimated coefficients of phi(B):\n');
disp(phiallV');
fprintf('Estimated coefficients of theta(B):\n');
disp(thetaallV');
fprintf('SD of noise: %f \n',SDz);
fprintf('AIC: %f \n',aicS);
fprintf('FPE: %f \n',fpeS);
fprintf('\t T \t\t NRMSE \n');
disp([[1:Tmax]' ARMAnrmseV]);
figure('Name', 'Normalized Root Mean Square Error for ARMA', 'NumberTitle', 'off');
plot([1:Tmax]',ARMAnrmseV,'.-k')
hold on
plot([1 Tmax],[1 1],'y')
xlabel('T')
ylabel('NRMSE')
title(sprintf('ARMA(%d,%d), fitting error',p,q))
% Further prediction error analysis
nlast = 0.5 * length(normEEG);
figure('Name', 'Normalized Root Mean Square Error for ARMA v2', 'NumberTitle', 'off');
[ARMAnrmseV,ARMApreM] = predictARMAnrmse(normEEG,p,q,Tmax,nlast,'ARMA');
figure("Name", "Comparison", "NumberTitle", "off");
plot([n-nlast+1:n]',normEEG(n-nlast+1:n),'.-')
hold on
plot([n-nlast+1:n]',ARMApreM(:,1),'.-r')
if Tmax>1
    plot([n-nlast+1:n]',ARMApreM(:,2),'.-c')
	if Tmax>2
        plot([n-nlast+1:n]',ARMApreM(:,3),'.-k')
    end
end
switch Tmax
    case 1
        legend('true','T=1','Location','Best')
    case 2
        legend('true','T=1','T=2','Location','Best')
    otherwise
        legend('true','T=1','T=2','T=3','Location','Best')
end
% multistep predictions
n1 = length(normEEG) - Tmax;
figure("Name", "multistep figure", "NumberTitle", "off");
[preV] = predictARMAmultistep(normEEG,n1,p,q,Tmax,'normEEG');
%% Non-linear analysis

% Permutation process to create 40 different iid datasets based on the
% stationary timeserie
dataset40 = zeros(length(normEEG), 40); % the new iid data holder
for i = 1:40
    n = randperm(length(normEEG));
    for j = 1:length(n)
        dataset40(j, i) = normEEG(n(j));
    end
end

% Linear characteristic 1: autocorrelation

autocorrelations40 = zeros(lags, 40);
figure("Name", "Autocorrelations for non-linear analysis", "NumberTitle", "off");
for i = 1:40
    autocorrelations40(:, i) = autocorr(dataset40(:, i), 'NumLags', lags);
    hold on;
end
% fprintf('Estimated autocorrelations of 40 iid samples:\n');
% disp(autocorrelations40');


% Linear characteristic 2: partial autocorrelation
partialautocorrelations40 = zeros(lags, 40);
figure("Name", "Partial autocorrelations for non-linear analysis", "NumberTitle", "off");
for i = 1:40
    partialautocorrelations40(:, i) = parcorr(dataset40(:, i), 'NumLags', lags);
    hold on;
end
% fprintf('Estimated partial autocorrelations of 40 iid samples:\n');
% disp(partialautocorrelations40');


% Non-linear characteristic 1: correlation dimension
tau = 1;
mmax = 5;
tittxt = 'normEEG corellation dimension';
fac = 0.5;
logrmin = 1;
logrmax = 10;
[rcM,cM,rdM,dM,nuM] = correlationdimension(normEEG,tau,mmax,tittxt,fac,logrmin,logrmax);



% Non-linear characteristic 2: Lyapunov exponents
tau = 1;
mmin = 2;
mmax = 5;
tittxt = 'Lyapunov exponents for normEEG';
nitecal = 20;
itewidth = 5;
[l1M,sdl1V] = maxlyapunov(normEEG,tau,mmin,mmax,tittxt,nitecal,itewidth);






