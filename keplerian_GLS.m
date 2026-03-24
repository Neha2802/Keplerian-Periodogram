function [f,pow,fap,P_est,w_est,ecc_est,K_est, T0_est] = keplerian_GLS(t, y, f, e, fap_vals)
% [f,pow,fap] = gls(t, y, f, e, fap_vals)
% Calculates The generalized Lomb-Scargle periodogram (Zechmeister 2009).
% If errors are not given, then unweighted means are used.
%
% Input: t - time vecotor of the data
%        y - the data
%        f - frequencies vector to calculate GLS (optional)
%        e - errors of the data (optional)
%        fap_vals - GLS power values to calculate the fap (default = [0.1 0.01 0.001])
%
% Output: f - frequencies.
%         pow - normalized GLS power for each f.
%         fap - selected false-alarm probability GLS values.
%
% Created: 2015, Shay Zucker
% Modified: 2018, Lev Tal-Or & Mathias Zechmeister

nslots = 12;
mypool = parpool(nslots);

% data
N = length(t);
% center time
tmin = min(t);
t = t - tmin;
tbase = max(t);

% default f:
if nargin<3
    ofac = 10;
    hifac = 1;
    fstep = 1 / tbase / ofac ; % frequency sampling depends on the time span, default for start frequency
    fbeg = fstep;
    fnyq = 0.5 / tbase * N;     % Nyquist frequency
    fend = fnyq * hifac;
    f = fbeg:fstep:fend;
end

% Allocate output array.
%
pow = 0 * f;
ecc_est = 0*f;
w_est = 0*f;
K_est = 0*f;
P_est = 0*f;
T0_est = 0*f;

% Convert frequencies to angular frequencies.
%
w = 2 * pi * f;

% P = w/(2*pi);
r = length(f); % length of frequency spectrum
disp(r);

T0 = t(N)*(0:1/r:1);
T0 = T0(2:end);

% M = 2*pi*(t - T0)/P;

ecc = 1*(0:1/r:(r-1)/r);


%--------------------------------------------------------------------------
% calculate the FAP lines (analytic aproximation):
%
frange = max(f) - min(f);
M = frange * tbase;
if nargin<5
    fap_vals = [0.1 0.01 0.001];
end
P = 1 - (1 - fap_vals).^(1/M);
fap = 1 - P.^(2/(N-3));

%--------------------------------------------------------------------------
if nargin<4 % if erors are not given use means
%
% Calculate values that are independent of frequency
%
Y = mean(y);
YY_hat = mean(y.^2);
YY = YY_hat - Y^2;


parfor k = 1:length(f)
    pow_max = -inf;
    P = w(k)/(2*pi);
    
    for j = 1:r
        M = 2*pi*(t - T0(j))/P;
        for i = 1:r
            E = kepSolver(M,ecc(i));

            nu = 2*atan(sqrt((1+ecc(i))/(1-ecc(i)))*tan(E/2));
            % x = w(k) * t;   % phases
            x = nu; % true anomaly
        
            cosx = cos(x);
            sinx = sin(x);
            cos2x = cos(2*x);
            sin2x = sin(2*x);
        
            omegatau = atan2(mean(sin2x)-2*mean(cosx)*mean(sinx),...
                mean(cos2x)-(mean(cosx)^2-mean(sinx)^2))/2;
        
            cosx_tau = cos(x-omegatau);
            sinx_tau = sin(x-omegatau);
        
            C = mean(cosx_tau);
            S = mean(sinx_tau);
        
            YC_hat = mean(y.*cosx_tau);
            YS_hat = mean(y.*sinx_tau);
        
            YC = YC_hat - Y*C;
            YS = YS_hat - Y*S;
        
            CC_hat = mean(cosx_tau.^2);
            SS_hat = mean(sinx_tau.^2);
            CS_hat = mean(cosx_tau.*sinx_tau);
        
            CC = CC_hat - C^2;
            SS = SS_hat - S^2;
            CS = CS_hat - C*S;
        
            D = CC*SS - CS^2;
            
            power = 1/YY * (YC^2/CC+YS^2/SS);
            % pow(k) = 1/YY * (YC^2/CC+YS^2/SS);
        
            

            if power > pow_max
                % pow_temp(i) = 1/YY * (YC^2/CC+YS^2/SS);
                pow_max = power;
                pow(k) = power;
                a = (YC*SS - YS*CS)/D;
                b = (YS*CC - YC*CS)/D;
            
                w_est(k) = atan(-b/a);
                K_est(k) = a/cos(w_est(k));
                ecc_est(k) = ecc(i);
                P_est(k) = P;
                T0_est(k) = T0(j);
            end

        end
        
    end
end

%--------------------------------------------------------------------------
else % if erors are given use wmeans
%
% Calculate values that are independent of frequency
%
wei = 1./e.^2;
Y = wmean(y,wei);
YY_hat = wmean(y.^2,wei);
YY = YY_hat - Y^2;

for k = 1:length(f)
    pow_max = -inf;
    P = w(k)/(2*pi);
    
    for j = 1:r
        M = 2*pi*(t - T0(j))/P;
        for i = 1:r
            E = kepSolver(M,ecc(i));

            nu = 2*atan(sqrt((1+ecc(i))/(1-ecc(i)))*tan(E/2));
            % x = w(k) * t;   % phases
            x = nu; % true anomaly

        
            cosx = cos(x);
            sinx = sin(x);
            cos2x = cos(2*x);
            sin2x = sin(2*x);
        
            omegatau = atan2(wmean(sin2x,wei)-2*wmean(cosx,wei)*wmean(sinx,wei),...
                wmean(cos2x,wei)-(wmean(cosx,wei)^2-wmean(sinx,wei)^2))/2;
        
            cosx_tau = cos(x-omegatau);
            sinx_tau = sin(x-omegatau);
        
            C = wmean(cosx_tau,wei);
            S = wmean(sinx_tau,wei);
        
            YC_hat = wmean(y.*cosx_tau,wei);
            YS_hat = wmean(y.*sinx_tau,wei);
        
            YC = YC_hat - Y*C;
            YS = YS_hat - Y*S;
        
            CC_hat = wmean(cosx_tau.^2,wei);
            SS_hat = wmean(sinx_tau.^2,wei);
            CS_hat = wmean(cosx_tau.*sinx_tau,wei);
        
            CC = CC_hat - C^2;
            SS = SS_hat - S^2;
            CS = CS_hat - C*S;
        
            D = CC*SS - CS^2;
        
            power = 1/YY * (YC^2/CC+YS^2/SS);
            % pow(k) = 1/YY * (YC^2/CC+YS^2/SS);
        
            if power > pow_max
                pow_max = power;
                pow(k) = power;
                a = (YC*SS - YS*CS)/D;
                b = (YS*CC - YC*CS)/D;
            
                w_est(k) = atan(-b/a);
                K_est(k) = a/cos(w_est(k));
                ecc_est(k) = ecc(i);
                P_est(k) = P;
                T0_est(k) = T0(j);
            end
        end
    end
end

disp("Periodogram done");

delete(mypool)

figure('Position', [100, 100, 900, 400]);
    
    subplot(1,2,1);
    plot(f, pow, 'b-');
    xlabel('Frequency');
    ylabel('Power');
    title('GLS Periodogram');
    grid on;
    
    subplot(1,2,2);
    period = 1./f;
    plot(period, pow, 'b-');
    xlabel('Period');
    ylabel('Power');
    title('Period Domain');
    set(gca, 'XScale', 'log');
    grid on;
    
    % Mark best period
    [~, idx] = max(pow);
    best_p = period(idx);
    xline(best_p, 'r--', sprintf('P=%.3f', best_p), 'LabelOrientation', 'horizontal');

end


% % % Plot the result
% % figure;
% % plot(theta, f, 'LineWidth', 1.5);
% % grid on;
% % xlabel('\theta (radians)');
% % ylabel('f(\theta)');
% % title(sprintf('Kapteyn Fourier Series (N=%d, z=%.2f)', N, z));
