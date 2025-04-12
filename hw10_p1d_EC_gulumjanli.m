function p1d_finite_time_LQR()
    clear; clc; close all;

    % System parameters
    alpha = 7.14;     
    beta  = 286.3;    
    A = [0  1;
         0 -alpha];
    B = [0; beta];

    % LQR weighting
    Q = [1  0;
         0  0];
    rho = 0.5;
    R = rho;

    x0 = [0.3; 0];

    [tInf, xInf, uInf] = infiniteHorizonResponse(A,B,Q,R,x0,0.20);

    horizon_list = [0.05, 0.20];  % 50 ms and 200 ms

    solStruct = struct();

    for iH = 1:length(horizon_list)
        Tf = horizon_list(iH);

        % Backward Riccati
        tSpanBack = [Tf 0];
        P0 = zeros(4,1);
        [tBackward, Psol] = ode45(@(t,Pvec) riccatiBackward(t,Pvec,A,B,Q,R), tSpanBack, P0);

        tBackwardFlipped = flipud(tBackward);
        PsolFlipped      = flipud(Psol);

        % Forward state
        tSpanFwd = [0 Tf];
        [tForward, xForward] = ode45(@(t,x) forwardStateODE(t,x,A,B,R,tBackwardFlipped,PsolFlipped), ...
                                     tSpanFwd, x0);

        % Control u(t)
        uForward = zeros(size(tForward));
        for k = 1:length(tForward)
            Pmat_k = getPmatrix(tForward(k), tBackwardFlipped, PsolFlipped);
            K_k = (1/rho)*B'*Pmat_k;
            uForward(k) = - K_k * xForward(k,:)';
        end

        % Store
        solStruct(iH).Tf        = Tf;
        solStruct(iH).tBackward = tBackwardFlipped;
        solStruct(iH).Psol      = PsolFlipped;
        solStruct(iH).tForward  = tForward;
        solStruct(iH).xForward  = xForward;
        solStruct(iH).uForward  = uForward;
    end

    %  Position Plot 
    hf1 = figure('Position',[100 100 800 600]);
    hf1.Color = 'w'; hold on;

    plot(tInf*1000, xInf(:,1), 'k--', 'LineWidth',2, 'DisplayName','Infinite-Horizon $\theta(t)$');

    for iH = 1:length(horizon_list)
        Tf  = solStruct(iH).Tf;
        tF  = solStruct(iH).tForward;
        xF  = solStruct(iH).xForward;
        lbl = sprintf('Finite-Horizon (T=%.0f ms) $\\theta(t)$', Tf*1000);

        % Custom color
        lineColor = [0.1+0.4*iH, 0.2, 0.8-0.3*iH];
        plot(tF*1000, xF(:,1), '-o', ...
             'LineWidth',2, 'MarkerSize',6, ...
             'Color', lineColor, ...
             'MarkerFaceColor', lineColor, ...
             'DisplayName', lbl);

        if abs(Tf - 0.05) < 1e-6
            xline(Tf*1000, ':', 'LineWidth',1, 'Color', [0.3 0.3 0.3]);
        end
    end

    xlabel('Time [ms]', 'Interpreter','latex', 'FontSize',14);
    ylabel('$\theta(t)$ [rad]', 'Interpreter','latex', 'FontSize',14);
    legend('Interpreter','latex', 'FontSize',12, 'Location','best');
    % %title('Position $\theta(t)$ Comparison: Finite vs. Infinite Horizon', ...
    %       'Interpreter','latex', 'FontSize',16);
    box on

    %  Control Plot 
    hf2 = figure('Position',[130 130 800 600]);
    hf2.Color = 'w'; hold on;

    plot(tInf*1000, uInf, 'k--', 'LineWidth',2, 'DisplayName','Infinite-Horizon $u(t)$');

    for iH = 1:length(horizon_list)
        Tf  = solStruct(iH).Tf;
        tF  = solStruct(iH).tForward;
        uF  = solStruct(iH).uForward;
        lbl = sprintf('Finite-Horizon (T=%.0f ms) $u(t)$', Tf*1000);

        lineColor = [0.1+0.4*iH, 0.2, 0.8-0.3*iH];
        plot(tF*1000, uF, '-o', ...
             'LineWidth',2, 'MarkerSize',6, ...
             'Color', lineColor, ...
             'MarkerFaceColor', lineColor, ...
             'DisplayName', lbl);

        if abs(Tf - 0.05) < 1e-6
            xline(Tf*1000, ':', 'LineWidth',1, 'Color', [0.3 0.3 0.3]);
        end
    end

    xlabel('Time [ms]', 'Interpreter','latex', 'FontSize',14);
    ylabel('$u(t)$ [Volts]', 'Interpreter','latex', 'FontSize',14);
    legend('Interpreter','latex', 'FontSize',12, 'Location','best');
    % title('Control $u(t)$ Comparison: Finite vs. Infinite Horizon', ...
    %       'Interpreter','latex', 'FontSize',16);
    box on

    exportgraphics(hf1,'p1d-finite-horizon-state.pdf','ContentType','vector');
    exportgraphics(hf2,'p1d-finite-horizon-control.pdf','ContentType','vector');
end

function dPvecdt = riccatiBackward(t, Pvec, A, B, Q, R)
    P = [Pvec(1) Pvec(2);
         Pvec(3) Pvec(4)];
    invR = 1/R;
    dP = -(Q + A'*P + P*A - P*B*(invR)*(B')*P);
    dPvecdt = [dP(1,1); dP(1,2); dP(2,1); dP(2,2)];
end

function dxdt = forwardStateODE(t, x, A, B, R, tP, Psol)
    Pmat_t = getPmatrix(t, tP, Psol);
    Kt = (1/R)*(B')*Pmat_t;
    dxdt = (A - B*Kt)*x;
end

function Pmat = getPmatrix(tquery, tGrid, PmatGrid)
    p1 = interp1(tGrid, PmatGrid(:,1), tquery, 'spline');
    p2 = interp1(tGrid, PmatGrid(:,2), tquery, 'spline');
    p3 = interp1(tGrid, PmatGrid(:,3), tquery, 'spline');
    p4 = interp1(tGrid, PmatGrid(:,4), tquery, 'spline');
    Pmat = [p1 p2; p3 p4];
end

function [tSol, xSol, uSol] = infiniteHorizonResponse(A,B,Q,R,x0,Tsim)
    Kinf = lqr(A,B,Q,R);
    Acl = A - B*Kinf;
    tSpan = linspace(0, Tsim, 400);
    [tSol,xSol] = ode45(@(t,x) Acl*x, tSpan, x0);

    uSol = zeros(size(tSol));
    for k=1:length(tSol)
       uSol(k) = -Kinf*xSol(k,:)';
    end
end
