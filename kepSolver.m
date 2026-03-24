% Function to solve Kepler's equation for a single set of parameters
function E = keplerSolver(M, e)
    % t: time vector
    % T0: time of passage of periastron
    % M: mean anomaly
    % e: eccentricity
    
    E = 0; % Initial guess
    tolerance = 1e-10;
    maxIter = 100;
    
    for iter = 1:maxIter
        f = E - e * sin(E) - M;
               
        f_prime = 1 - e * cos(E);
   
        % delta = f / f_prime;
        delta = f_prime\f;
     
        E = E - delta;
        
        if abs(delta) < tolerance
            break;
        end
    end
end

% % E - e*sinE = M;
% 
% % E = M + sum_{n=1} (2/n)Jn(nz)sin(n*M) 
% % Kapteyn Fourier Series Example in MATLAB
% % f(theta) = sum_{n=1}^N a_n * J_n(n*z) * cos(n*theta)
% 
% % Parameters
% n = 20;              % Number of terms in the series
% z = 0.5;             % Parameter for Bessel function
% 
% % Coefficients a_n (example: all ones)
% a = (2/n)*ones(1, n);
% 
% % Initialize series sum
% E = zeros(size(M));
% 
% % Compute the series
% for i = 1:n
%     % Bessel function of the first kind J_n(n*z)
%     Jn = besselj(i, i*z);
% 
%     % Add term to the sum
%     E = E + a(i) * Jn * sin(i * M);
% end