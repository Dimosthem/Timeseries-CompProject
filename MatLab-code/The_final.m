%% Timeseries analysis project
 
%Bounarelis Dimosthenis
%Giachoudis Christos
 
 
%% Import data and plot

% Change the filename so that you can load the right file. The file must
% exist at the same folder as the script file
filename = 'normalEEG.dat';
data = load (filename);
 
 
 
% Plot the data
x1 = 1:length(data);
y1 = data;
figure('Name', 'normalEEG', 'NumberTitle', 'off');
% figure('Name', 'epilepticEEG', 'NumberTitle', 'off');
plot(x1, y1, 'color', '#A2142F');
 
%or
 
% Plot the data in parts(just for better visual)
figure('Name', 'normalEEG', 'NumberTitle', 'off');
% figure('Name', 'epilepticEEG', 'NumberTitle', 'off');
tittxt = 'data';
nseg = 3;
type = 'b';
tau_s = 1;
% Be sure that pser function file is at the same folder as the script file
pser(data, tittxt, nseg, type, tau_s);
 
% _____Mean value_____and_____variance(and standard deviation)_____
data_mean = mean(data);
formatSpec = '******The mean value of the normal EEG is %4.4f\n';
fprintf(formatSpec, data_mean);
data_var = var(data);
formatSpec = '******The variance of the normal EEG is %4.4f\n';
fprintf(formatSpec, data_var);
data_std = std(data);
formatSpec = '******The variance of the normal EEG is %4.4f\n';
fprintf(formatSpec, data_std);


%% Linear analysis


% _____Stationarity_____(question 1)
stationary_data = data;
figure('Name', 'stationarity before', 'NumberTitle', 'off');
movMean = movmean(data, 100);
plot(movMean);
choice = input('Do you want to remove the timeseries trend? Choose 1(to detrend) or 0(to continue without detrending) \n');
if choice == 1
    stationary_data = data(2:end) - data(1:end-1);
    %Plot the stationary data
    xS = 1:length(stationary_data);
    yS = stationary_data;
    figure('Name', 'Stationary Timeseries of normal EEG', 'NumberTitle', 'off');
    plot(xS,yS, 'color', '#A2142F');
end
figure('Name', 'stationarity after', 'NumberTitle', 'off');
movMean = movmean(stationary_data, 100);
plot(movMean);

% _____Check for unstable variance_____(question 2)
% Negative inputs for box-cox are invalid. Thus, we will fix that problem
% with a constant.
figure('Name', 'Check for unstable variance', 'NumberTitle', 'off');
histfit(stationary_data, 100, 'Normal')
transdat = stationary_data;
choice = input('Do you want to stabilize the timeseries variance? Choose 1(for stabilizing) or 0(for not-stabilizing) \n');
if choice == 1
    constant = abs(min(stationary_data)) + 0.0000001;
    stationary_data = stationary_data + ones(size(stationary_data))*constant;
    % Box-Cox method for stabilizing standard deviation/variance
    [transdat,lambda] = boxcox(stationary_data); % make the best BoxCox transformation
    % lets see whether or not boxcox was 'useless' 
    figure('Name', 'After stabilizing variance', 'NumberTitle', 'off');
    histfit (transdat, 100, 'Normal')
    if lambda == 0
        formatSpec = 'We use the following transformation: transdat = log(stationary_data)\n';
        fprintf(formatSpec);
    else
        formatSpec = 'We use the following transformation: transdat = (stationary_data^lambda - 1)/lambda \n';
        fprintf(formatSpec);
    end
    %Plot the transformed data
    stationary_data = stationary_data - ones(size(stationary_data))*constant;
    xBC = 1:length(transdat);
    yBC = transdat;
    figure('Name', 'BOX-COX-transformed normalEEG', 'NumberTitle', 'off');
    % figure('Name', 'BOX-COX-transformed epilepticEEG', 'NumberTitle', 'off');
    plot(xBC, yBC, 'color', '#EDB120');
end


%_____Periodicity_____(question 3)
non_per = transdat;
figure('Name', 'Perriodogram of normalEEG', 'NumberTitle', 'off');
periodogram(transdat);
choice = input('Do you want to remove the timeseries periodicity? Choose 1(to remove periodicity) or 0(to not remove periodicity) \n');
if choice == 1
    per = input('Give the periodicity order of timeseries >');
    % Be sure that seasonalcomponents function file is at the same folder as the script file
    non_per = seasonalcomponents(transdat, per);
    %Plot the non-periodical data
    xnp = 1:length(non_per);
    ynp = non_per;
    figure('Name', 'Non-periodical normalEEG', 'NumberTitle', 'off');
    % figure('Name', 'Non-periodical epilepticEEG', 'NumberTitle', 'off');
    plot(xnp, ynp, 'color', '#EDB120');
end


%_____White Noise or more information ?_____(question 4)
% Check autocorrelations
lags = input('Give how many lags you want to compute for autocorrelation and partial autocorrelation>');
figure('Name', 'Autocorrelations of normalEEG', 'NumberTitle', 'off');
% figure('Name', 'Autocorrelations of epilepticEEG', 'NumberTitle', 'off');
autocorr(non_per, 'NumLags', lags);

% Check partial autocorrelations
figure('Name', 'Partial Autocorrelation of normalEEG', 'NumberTitle', 'off');
% figure('Name', 'Partial Autocorrelation of epilepticEEG', 'NumberTitle', 'off');
parcorr(non_per, 'NumLags', lags);
 
 
%_____Linear Models Simulations_____(question 5)
n = length(non_per);
% ***ARMA model***
p = input('Give the order p of the AR part of ARMA >'); % if you want to create an MA model set p to zero
q = input('Give the order q of the MA part of ARMA >'); % if you want to create an AR model set p to zero
Tmax = input('How many steps ahead you want the function to predict and then compute the errors ? >');
[ARMAnrmseV,phiallV,thetaallV,SDz,aicS,fpeS,model]=fitARMA(non_per, p, q, Tmax);
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
% Further prediction res analysis
nlast = 0.5 * length(non_per);
figure('Name', 'Normalized Root Mean Square Error for ARMA v2', 'NumberTitle', 'off');
[ARMAnrmseV,ARMApreM] = predictARMAnrmse(non_per,p,q,Tmax,nlast,'ARMA');
figure('Name', 'Comparison', 'NumberTitle', 'off');
plot([n-nlast+1:n]',non_per(n-nlast+1:n),'.-')
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
n1 = length(non_per) - Tmax;
figure('Name', 'multistep figure', 'NumberTitle', 'off');
[preV] = predictARMAmultistep(non_per,n1,p,q,Tmax,'data');


%% Non-linear analysis
mdl = arima(p,0,q);
estMdl = estimate(mdl, non_per);
 
[res,V] = infer (estMdl, non_per);

%Permutation process to create 40 different iid datasets based on the
%stationary timeseries
dataset40 = zeros(length(res), 40); % the new iid data holder
for i = 1:40
    n = randperm(length(res));
    for j = 1:length(n)
        dataset40(j, i) = res(n(j));
    end
end
 
% Linear characteristic 1: autocorrelation
 
figure('Name', 'Original autocorrelation', 'NumberTitle', 'off');
autocorr(res, 'NumLags', lags);
figure('Name', 'Autocorrelations for 1-10 datasets', 'NumberTitle', 'off');
for i = 1:10
    nexttile
    autocorr(dataset40(:, i), 'NumLags', lags);
 
end
figure('Name', 'Autocorrelations for 11-20 datasets', 'NumberTitle', 'off');
for i = 11:20
    nexttile
    autocorr(dataset40(:, i), 'NumLags', lags);
 
end
figure('Name', 'Autocorrelations for 21-30 datasets', 'NumberTitle', 'off');
for i = 21:30
    nexttile
    autocorr(dataset40(:, i), 'NumLags', lags);
 
end
figure('Name', 'Autocorrelations for 31-40 datasets', 'NumberTitle', 'off');
for i = 31:40
    nexttile
    autocorr(dataset40(:, i), 'NumLags', lags);
 
end
 
 
 
 
 
% Linear characteristic 2: partial autocorrelation
figure('Name', 'Original partial autocorrelations', 'NumberTitle', 'off');
parcorr(res, 'NumLags', lags);
partialautocorrelations40 = zeros(lags + 1, 40);
figure('Name', 'Partial autocorrelations for non-linear analysis', 'NumberTitle', 'off');
hold on;
for i = 1:40
    parcorr(dataset40(:, i), 'NumLags', lags);
 
end
hold off;
fprintf('Estimated partial autocorrelations of 40 iid samples:\n');
disp(partialautocorrelations40');
 
% Linear characteristic 3: Mean values
m = [mean(res) mean(dataset40)];
figure ('Name', 'Mean values')
plot (m)
% Linear characteristic 4: Variance
v = [var(res) var(dataset40)];
figure('Name', 'Variance values')
plot(v)


% Non-linear characteristic 1: correlation dimension
tau = 2;
mmax = 40;

n = calcCorDims(res, tau, mmax);

Ns = [];
for i = 1:40
    Ns = [Ns; calcCorDims(dataset40(:,i), tau, mmax)];
end

figure('Name' , 'Correlation dimension');


hold on
plot (n, 'Color', 'black');
for i = 1:40
    plot (Ns(i, :), 'Color', 'red')
end
hold off



% Non-linear characteristic 2: Lyapunov exponents
l = calcLyapExp(res, tau, mmax);
Ls = [];
for i = 1:40
    Ls = [Ls; calcLyapExp(dataset40(:,i), tau, mmax)];
end


figure('Name','Lyapunov Exponents')
hold on
plot (l, 'Color', 'black');
for i = 1:40
    plot (Ls(i, :), 'Color', 'red')
end
hold off
 
function n = calcCorDims(data, tau, mmax)
    n = [];
    for i = 1:mmax
        i
        n = [n correlationDimension(data, tau, i)];
    end
 
end
 
function n = calcLyapExp(data, tau, mmax)
    n = [];
    for i = 1:mmax
        i
        n = [n lyapunovExponent(data,'Lag',tau,'Dimension',i)];
    end
 
end