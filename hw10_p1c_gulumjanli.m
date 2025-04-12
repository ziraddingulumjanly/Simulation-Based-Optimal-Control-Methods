clear; clc; close all;

%  System parameters 
alpha = 7.14;    % rad/s
beta  = 286.3;   % rad/(V*s^2)
A = [0    1;
     0 -alpha];
B = [0; beta];

Q = [1 0;
     0 0];

%  Initial conditions & time horizon 
x0 = [0.3; 0];         
simulationTmax = 0.2;  
targetRiseTime = 0.05; 

rhoRange = logspace(-2, 2, 100);
bestRho       = NaN;
bestError     = Inf;
bestRiseTime  = NaN;

for rhoCandidate = rhoRange
    % Compute the LQR gain & closed-loop A
    K = lqr(A, B, Q, rhoCandidate);
    Acl = A - B*K;

    % Simulate closed-loop for 200 ms
    tSpan = linspace(0, simulationTmax, 1001);
    [tSol, xSol] = ode45(@(t,x) Acl*x, tSpan, x0);

    % Measure 10-90% "rise time" for theta(t)
    thetaSol = xSol(:,1);
    t90 = firstCrossingTime(tSol, thetaSol, 0.3 * 0.9); 
    t10 = firstCrossingTime(tSol, thetaSol, 0.3 * 0.1); 

    if ~isnan(t90) && ~isnan(t10) && (t10>t90)
        thisRiseTime = t10 - t90;      
        thisError    = abs(thisRiseTime - targetRiseTime);
        if thisError < bestError
            bestError     = thisError;
            bestRho       = rhoCandidate;
            bestRiseTime  = thisRiseTime;
        end
    end
end

fprintf('Best rho = %g\n', bestRho);
fprintf('Achieved 10--90%% rise time = %.3f ms (target 50 ms)\n',...
    bestRiseTime*1000);

%  Final simulation 
Kfinal   = lqr(A, B, Q, bestRho);
AclFinal = A - B*Kfinal;
[tFinal, xFinal] = ode45(@(t,x) AclFinal*x, [0 simulationTmax], x0);
thetaFinal = xFinal(:,1);

% Helper function 
function tC = firstCrossingTime(tVals, signal, threshold)
    idx = find(signal <= threshold, 1, 'first');
    if isempty(idx)
        tC = NaN;
    else
        tC = tVals(idx);
    end
end
