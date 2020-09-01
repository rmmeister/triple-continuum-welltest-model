%% INITIALIZE MATLAB
clc
clear all
close all
format long
%% DASHBOARD
beta = 0 ; % friction coefficient,
gamma = .01; % fluctuation coefficient
alpha = 1; % damping coefficient

kD = .99; % kappa: dimensionless permeability; kf/(kf + kv)

wf = 1e-3; % fracture storage contrast
wv = .9; % vug storage contrast
lmf = 1e-5; % lambda_mf: matrix-fracture interporosity
lmv = 1e-7; % lambda_mv: matrix-vug interporosity
lvf = 1e-7; % lamvda_vf: vug-fracture interporosity
wm = 1 - wf - wv; % matrix storage contrast
Sw = [1 1 1 1]; % wellbore skin
CwD =[0 1 10 20]; % wellbore storage
Slv = 0; % large vug skin
ClvD = 0; % large vug storage
rlvD = 1.5; % large vug radius
reD1 = 700; % external reservoir radius 1
reD2 = 2000; % external reservoir radius 2
hD = .75; % thickness ratio

% TIME VECTOR
tD = [1e-3:2e-4:1e-2-2e-4 .01:.002:.1-.002 .1:0.02:1-0.02 ...
    1:.2:10-.2 10:2:98 100:20:1e3-20 1e3:200:1e4-200 1e4:2e3:1e5-2e3 1e5:2e4:1e6-2e4 1e6:2e6:1e9-2e6];
intHandle = @(tD) (tD.^(alpha-1)).*exp(-tD);
Gamma = integral ( intHandle, 0, inf );

%% KEYS
key = 3; % (3) for Camacho-Velazquez (2005), (1) for Du et al. (2019)
keyPlot = 1; % (1) for log-log plot, (2) for semi-log plot
WR = 'off';
LT = 'off';
stehfestIndex = 6;

%% CALCULATE LAPLACE OF PD AND PD
figure('Position', [200 200 900 600])
% for jj = 1:2        
for j = 1:length(Sw)
    if j == 1 || j == 2
        WR = 'off';
    else
        WR = 'off';
    end
    PD = zeros(length(tD), 1);
    for i = 1:length(tD)
        PDbar = @(u) laplaceP ( Gamma, alpha, beta, gamma, lmf, lvf, lmv, wf, wv, ...
            wm, kD, Sw(j), Slv, u, rlvD, reD1, reD2, CwD(j), ClvD, hD, key, WR, LT );
        PD(i) = stehfestAlgorithm(PDbar, tD(i), stehfestIndex ); 
    end
    
    PDPrime = diff(PD)'./diff(tD).*tD(2:end);
    
    %% RESULTS
    col = {'k', 'r', 'b', 'g', 'y', 'c', 'm', [0.9290, 0.6940, 0.1250]};
    marker = ['*', 'd', '^', 'o', 's', 'o', 's'];
    markerCol = ['k', 'k', 'r', 'r', 'r', 'w', 'w'];
    if keyPlot == 1
        pressPlot(j) = loglog(tD, PD, 'LineWidth', 2, ...
            'DisplayName',...
            ['S_w = ', num2str(Sw(j)), '  C_{wD} = ', num2str(CwD(j))],...
            'Color', col{j} );
        hold on
        derPlot(j) = loglog(tD(2:end), PDPrime, 'LineWidth', 2, ...
            'HandleVisibility', 'off',...
            'Color', col{j}, 'LineStyle', '--');
        hold on
    elseif keyPlot == 2
        semilogPlot(j) = semilogx(tD, PD, 'LineWidth', 2, 'Color', col{j}, ...
            'DisplayName', ['\lambda_{mf} = ', num2str(lmf(j)), '  \kappa = ', num2str(kD(j))]);
        hold on
    end
end

legend('show', 'Location', 'Best');
grid on
xlabel('Dimensionless Time, t_D');
if keyPlot == 1
    ylim([1e-2 1e2]);
    xlim([1e-2 1e8]);
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'FontSize', 12, ...
        'Interpreter', 'latex', 'BackgroundColor', 'w' ...
        ,'String', ['$\lambda_{mf} = $', num2str(lmf),...
        ', $\lambda_{mv} = $', num2str(lmv),...
        ', $\lambda_{vf} = $', num2str(lmf)])
    annotation('textbox', [0.18, 0.72, 0.1, 0.1], 'FontSize', 12, ...
        'Interpreter', 'latex', 'BackgroundColor', 'w' ...
        ,'String', ['$\omega_f = $', num2str(wf),...
        ', $\omega_v = $', num2str(wv)])
    ylabel(['Dimensionless Pressure, P_D', char(10)' 'Pressure Derrivative, P^''_D']);
    
elseif keyPlot == 2
    ylim([0 10])
    xlim([1e-1 1e8]);
    ylabel('Dimensionless Pressure, P_D');
    annotation('textbox', [0.55, 0.8, 0.1, 0.1], 'FontSize', 12, ...
        'Interpreter', 'latex', 'BackgroundColor', 'w' ...
        ,'String', ['$w_f = $', num2str(wf),...
                    ', $w_v = $', num2str(wv),...
                    ', $\lambda_{mv} = $', num2str(lmv),...
                      '$\lambda_{vf} = $', num2str(lvf)])
end
