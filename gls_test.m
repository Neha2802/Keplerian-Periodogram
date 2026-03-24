clc 
clearvars
close all

fileID = 'C:/Users/vadin/Desktop/RV_nomp/harps_rvs_HD10180.xlsx';
% fileID = 'C:/Users/vadin/DesktopResearch\Codes';

opt = detectImportOptions(fileID);
data = readtable(fileID, opt);

bjd = table2array(data(:,4));
rv = table2array(data(:,5));
e_rv = table2array(data(:,6));
bis = table2array(data(:,11));
fwhm = table2array(data(:,12));
s_hk = table2array(data(:,13));
e_s_hk = table2array(data(:,14));
snr = table2array(data(:,15));

% rv_mean = mean(rv);
% bis_mean = mean(bis);
% fwhm_mean = mean(fwhm);
% s_hk_mean = mean(s_hk);
% 
% rv = rv - rv_mean.*ones(length(bjd),1);
% bis = bis - bis_mean.*ones(length(bjd),1);
% fwhm = fwhm - fwhm_mean.*ones(length(bjd),1);
% s_hk = s_hk - s_hk_mean.*ones(length(bjd),1);

R = 600;
f = (0:1/R:(R-1)/R)*2*pi;
for i = 1
    % [freq, pow, fap, P_est,w_est,ecc_est,K_est, T0_est] = keplerian_GLS(bjd, rv, f);
    % % [freq, pow, fap, P_est,w_est,ecc_est,K_est, T0_est] = keplerian_GLS(bjd, rv, f, e_rv);
    % % [freq, pow, fap] = gls(bjd,rv);
    % S.freq = freq;
    % S.pow = pow;
    % S.fap = fap;
    % S.K_est = K_est;
    % S.ecc_est = ecc_est;
    % S.T0_est = T0_est;
    % S.P_est = P_est;
    % S.w_est = w_est;
    % save("info.mat",'-struct','S');
    
    S = load("info.mat");
    % 
    figure('Position', [100, 100, 900, 400]);

    % subplot(1,2,1);
    % plot(S.freq, S.pow, 'b-');
    % xlabel('Frequency');
    % ylabel('Power');
    % title('GLS Periodogram');
    % grid on;

    % subplot(1,2,2);
    period = 1./S.freq;
    plot(period, S.pow, 'b-');
    xlabel('Period');
    ylabel('Power');
    title('Period Domain');
    set(gca, 'XScale', 'log');
    grid on;

    % Mark best period
    [~, idx] = max(S.pow);
    best_p = period(idx);
    xline(best_p, 'r--', sprintf('P=%.3f', best_p), 'LabelOrientation', 'horizontal');


end