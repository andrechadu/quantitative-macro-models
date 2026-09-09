// Based on: "Fiscal Backing, Inflation and US Business Cycles"
// Authors: Frank Smets and Raf Wouters
// Dynare implementation: André Chadú

// 2.2. A New Keynesian model with partial fiscal backing
var 
    y pi R Rb b tau ystar Pb r         // monetary-led economy
    yF piF RF RbF bF tauF yFstar PbF  // fiscal-led shadow economy
    eps_a eps_d;                      // AR(1) shocks modeled as endogenous


varexo e_a e_d eps_mp eps_tau;



// === Parameter block ===
parameters 
    kappa beta r_par rho_m b_par
    psi_pi psi_piF psi_y rho_R
    delta_b delta_bF delta_y delta_dy rho_tau
    lambda;


// === CALIBRAÇÃO PARA O BRASIL (Baseada na Tabela 9) ===
kappa     = 0.0075;    // Phillips Curve Slope (IPEA)
beta      = 0.9879;    // Discount Factor (BCB - Neutral real rate ~5%)
r_par     = 1/0.9806;  // Gross Nominal Rate (CMN - Inflation target 3%)
rho_m     = 0.939;     // Debt Maturity (National Treasury ~4.1 years)
b_par     = 3.2;       // Debt/GDP Ratio (IFI Projections ~80%)
psi_pi    = 1.40;      // Taylor Rule: Inflation Response (BCB WP)
psi_piF   = 0.0;       // Theoretical assumption for fiscal-led regime
psi_y     = 0.10;      // Taylor Rule: Output Gap Response (BCB WP)
rho_R     = 0.83;      // Taylor Rule: Interest Rate Smoothing (BCB WP)
delta_b   = 0.02;      // Fiscal Rule: Debt Response (IPEA)
delta_bF  = 0.0;       // Theoretical assumption for fiscal-led regime
delta_y   = 0.4;       // Fiscal Rule: Output Level Response (IPEA)
delta_dy  = 0.3;       // Fiscal Rule: Output Growth Response (IPEA)
rho_tau   = 0.9;       // Fiscal Rule: Smoothing (SW, 2024)
lambda    = 0.5;       // Experiment parameter, not a calibrated value

// === Linearized model block ===

model(linear);

// Definição da taxa de juros real
r = R - pi(+1);

//// === Regime principal (monetary-led): (11) a (16) ===

// (11) IS curve
y = y(+1) - (R - pi(+1)) + eps_d;

// (12) Phillips curve
pi = kappa * (y - ystar) + beta * pi(+1);

// (13) Output target
ystar = eps_a;

// (14) No arbitrage
R = Rb(+1);

// (15) Return on long-term bond — linearized version
Rb = (rho_m / r_par) * Pb - Pb(-1);

// (16) Government budget constraint
b = (1 / beta) * b(-1) + b_par * (1 / beta) * (Rb - y + y(-1) - pi) - tau;

//// === Policy rules (monetary-led): Equações (17) e (18)

// (17) Monetary policy rule
R = rho_R * R(-1)
  + (1 - rho_R) * (psi_pi * (pi - piF) + psi_piF * piF + psi_y * (y - ystar))
  + eps_mp;

// (18) Fiscal policy rule
tau = rho_tau * tau(-1)
   + (1 - rho_tau) * (delta_b * (b(-1) - bF(-1)) + delta_bF * bF(-1) 
   + delta_y * y) + delta_dy * (y - y(-1)) - eps_tau;

//// === Shadow economy (fiscal-led): (11F) a (16F)

// (11F) IS curve
yF = yF(+1) - (RF - piF(+1)) + (1 - lambda) * eps_d;

// (12F) Phillips curve
piF = kappa * (yF - yFstar) + beta * piF(+1);

// (13F) Output target
yFstar = (1 - lambda) * eps_a;

// (14F) No arbitrage
RF = RbF(+1);

// (15F) Return on long-term bond — linearized
RbF = (rho_m / r_par) * PbF - PbF(-1);

// (16F) Government budget constraint
bF = (1 / beta) * bF(-1) + b_par * (1 / beta) * (RbF - yF + yF(-1) - piF) - tauF;

//// === Policy rules (fiscal-led): Equações (19) e (20)

// (19) Monetary policy rule in fiscal-led regime
RF = rho_R * RF(-1)
   + (1 - rho_R) * (psi_piF * piF + psi_y * (yF - yFstar))
   + (1 - lambda) * eps_mp;

// (20) Fiscal policy rule in fiscal-led regime
tauF = rho_tau * tauF(-1)
     + (1 - rho_tau) * (delta_bF * bF(-1) + delta_y * yF)
     + delta_dy * (yF - yF(-1))
     - (1 - lambda) * eps_tau;


//// === AR(1) shocks ===
eps_a = 0.9 * eps_a(-1) - e_a;
eps_d = 0.9 * eps_d(-1) + e_d;
end;

initval;
    y     = 0; pi    = 0;
    R     = 0; Rb    = 0;
    b     = 0; tau   = 0;
    ystar = 0; Pb    = 0;

    yF     = 0; piF    = 0;
    RF     = 0; RbF    = 0;
    bF     = 0; tauF   = 0;
    yFstar = 0; PbF    = 0;

    eps_a  = 0;
    eps_d  = 0;
end;

//steady;

resid;

check;

//model_diagnostics;


shocks;
    var e_a;     stderr 1.0;
    var e_d;     stderr 1.0;
    var eps_mp;  stderr 0.01; // Choques de política monetária são usualmente menores
    var eps_tau; stderr 1.0;
end;

%==========================================================================
% SCRIPT DE ANÁLISE FLEXÍVEL
%==========================================================================

% --- CONTROLE DO EXPERIMENTO ---
% Mude esta variável para escolher qual figura/análise gerar:
% 1 = Varia LAMBDA para todos os choques (Figuras 2, 5, 6, 7)
% 2 = Varia KAPPA para o choque fiscal (Figura 3)
% 3 = Varia PSI_PI_F para o choque fiscal (Figura 4)
experimento_a_rodar = 3; % Mude para 1, 2 ou 3

% --- CONFIGURAÇÕES GERAIS ---
var_list_to_plot = {'y', 'pi', 'R', 'tau', 'r', 'b'};
title_list = {'Output', 'Inflation', 'Nominal rate', 'Primary surplus', 'Real rate', 'Debt ratio'};
irf_periods = 40; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%==========================================================================
% INÍCIO DOS EXPERIMENTOS
%==========================================================================

if experimento_a_rodar == 1
    % --- EXPERIMENTO 1: Variação de Lambda (estilo gradiente vermelho) ---
    disp('Rodando Experimento 1: Variando lambda para todos os choques...');

    shocks_to_plot = {'e_a', 'e_d', 'eps_mp', 'eps_tau'};

    % Mais resolucao de lambda para produzir gradiente suave
    lambda_values = linspace(1.0, 0.0, 11);  % 1.00, 0.90, ..., 0.00
    results = struct();

    % --- Simulacoes Dynare ---
    for i = 1:length(lambda_values)
        set_param_value('lambda', lambda_values(i));
        stoch_simul(order=1, irf=40, nograph, noprint, irf_shocks=(e_a, e_d, eps_mp, eps_tau));

        for shock_idx = 1:length(shocks_to_plot)
            current_shock_name = shocks_to_plot{shock_idx};
            for var_idx = 1:length(var_list_to_plot)
                current_var_name = var_list_to_plot{var_idx};
                irf_field_name = [current_var_name, '_', current_shock_name];
                results.(current_shock_name).(current_var_name)(:, i) = oo_.irfs.(irf_field_name)'; %#ok<AGROW>
            end
        end
    end

    % --- Paleta: gradiente vermelho claro -> vermelho escuro (neutro) ---
    nL = length(lambda_values);
    startRGB = [1.00, 0.85, 0.85];   % rosa bem claro
    stopRGB  = [0.55, 0.00, 0.00];   % vermelho escuro
    cmap = [linspace(startRGB(1), stopRGB(1), nL)', ...
            linspace(startRGB(2), stopRGB(2), nL)', ...
            linspace(startRGB(3), stopRGB(3), nL)'];
    cmap = flipud(cmap);   % << agora λ=1 (primeira coluna) fica mais escuro

    line_w = 1.25;

    % --- Plots ---
    for shock_idx = 1:length(shocks_to_plot)
        current_shock_name = shocks_to_plot{shock_idx};
        fig = figure('Name', ['IRFs to shock: ', strrep(current_shock_name, '_', ' ')], 'Color', 'w');

        % guardar os eixos para escolher onde colocar a legenda
        ax_array = gobjects(1, length(var_list_to_plot));

        for var_idx = 1:length(var_list_to_plot)
            current_var_name = var_list_to_plot{var_idx};
            ax_array(var_idx) = subplot(2, 3, var_idx);
            ax = ax_array(var_idx);

            data_to_plot = results.(current_shock_name).(current_var_name);
            irf_periods = size(data_to_plot, 1);

            % aplicar gradiente como ordem de cores do eixo
            set(ax, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
            hold(ax, 'on');

            % linhas reais NÃO entram na legenda
            h_lines = plot(ax, 1:irf_periods, data_to_plot, 'LineWidth', line_w);
            set(h_lines, 'HandleVisibility', 'off');

            yline(ax, 0, 'k-');
            xlim(ax, [0, irf_periods]);
            grid(ax, 'on'); box(ax, 'off');
            title(ax, title_list{var_idx}, 'Interpreter', 'latex');
            xlabel(ax, 'Periods', 'Interpreter', 'latex');
            ylabel(ax, 'PP dev from SS', 'Interpreter', 'latex');
            hold(ax, 'off');
        end

        % ===== LEGENDA (posição dinâmica dentro do Real rate) =====
        ax_for_legend = ax_array(5);   % Real rate

        target_lambdas = [1.0 0.7 0.5 0.2 0.0];
        idx_leg = arrayfun(@(v) find(abs(lambda_values - v) == min(abs(lambda_values - v)), 1, 'first'), target_lambdas);
        idx_leg = unique(idx_leg, 'stable');

        h_dummy = gobjects(numel(idx_leg), 1);
        hold(ax_for_legend, 'on');
        for k = 1:numel(idx_leg)
            j = idx_leg(k);
            h_dummy(k) = plot(ax_for_legend, NaN, NaN, '-', ...
                              'Color', cmap(j,:), 'LineWidth', line_w);
        end
        lgd_labels = arrayfun(@(j) sprintf('$\\lambda=%.1f$', lambda_values(j)), ...
                              idx_leg, 'UniformOutput', false);

        lgd = legend(ax_for_legend, h_dummy, lgd_labels, ...
                     'Location', 'best', ...     % << posição dinâmica
                     'Interpreter', 'latex', ...
                     'FontSize', 8, ...
                     'Box', 'off', ...
                     'NumColumns', 1);           % opcional: 1 coluna
        lgd.AutoUpdate = 'off';                  % evita mudanças futuras
        hold(ax_for_legend, 'off');
    end

elseif experimento_a_rodar == 2
    % --- EXPERIMENTO 2: Variação de Kappa (Figura 3) ---
    disp('Rodando Experimento 2: Variando kappa para o choque fiscal (Figura 3)...');

    % >>>>> AJUSTE CRÍTICO PARA CORRIGIR A DINÂMICA <<<<<
    % 1. Salvar o valor original de psi_pi
    psi_pi_original = M_.params(find(strcmp(M_.param_names, 'psi_pi')));
    % 2. Definir o valor temporário para este experimento
    set_param_value('psi_pi', 1.5);

    kappa_values = [0.10, 0.05, 0.20];   % [Baseline, Flatter, Steeper]
    lambda_fixo  = 0.5;

    y_results   = zeros(irf_periods, length(kappa_values));
    pi_results  = zeros(irf_periods, length(kappa_values));
    R_results   = zeros(irf_periods, length(kappa_values));
    tau_results = zeros(irf_periods, length(kappa_values));
    r_results   = zeros(irf_periods, length(kappa_values));
    b_results   = zeros(irf_periods, length(kappa_values));

    set_param_value('lambda', lambda_fixo);

    for i = 1:length(kappa_values)
        set_param_value('kappa', kappa_values(i));
        stoch_simul(order=1, irf=40, nograph, noprint, irf_shocks=(eps_tau));

        y_results(:,i)   = oo_.irfs.y_eps_tau';
        pi_results(:,i)  = oo_.irfs.pi_eps_tau';
        R_results(:,i)   = oo_.irfs.R_eps_tau';
        tau_results(:,i) = oo_.irfs.tau_eps_tau';
        r_results(:,i)   = oo_.irfs.r_eps_tau';
        b_results(:,i)   = oo_.irfs.b_eps_tau';
    end

    % >>>>> AJUSTE CRÍTICO: 3. Restaurar o valor original de psi_pi <<<<<
    set_param_value('psi_pi', psi_pi_original);

    results_cell = {y_results, pi_results, R_results, tau_results, r_results, b_results};

    % -------- Paleta: rosa (κ baixo) -> vermelho escuro (κ alto) --------
    nK = length(kappa_values);
    startRGB = [1.00, 0.85, 0.85];    % claro
    stopRGB  = [0.55, 0.00, 0.00];    % escuro
    cmap_sorted = [linspace(startRGB(1), stopRGB(1), nK)', ...
                   linspace(startRGB(2), stopRGB(2), nK)', ...
                   linspace(startRGB(3), stopRGB(3), nK)'];   % em ordem ascendente de κ
    % atribui cor por ranking de κ (menor=claro, maior=escuro) preservando a ordem original
    [~, ascIdx] = sort(kappa_values, 'ascend');   % índices de κ do menor p/ maior
    colors_in_order = zeros(nK,3);
    colors_in_order(ascIdx, :) = cmap_sorted;

    line_w = 1.5;

    figure('Name', 'Figure 3: Transfer shock with different Phillips curve slopes', 'Color', 'w');

    ax_array = gobjects(1,6);
    for ii = 1:6
        ax_array(ii) = subplot(2, 3, ii);
        ax = ax_array(ii);

        set(ax, 'ColorOrder', colors_in_order, 'NextPlot', 'replacechildren');
        hold(ax, 'on');

        % plota nas cores já alinhadas à ordem original [Baseline, Flatter, Steeper]
        h_lines = plot(ax, 1:irf_periods, results_cell{ii}, 'LineWidth', line_w);
        set(h_lines, 'HandleVisibility', 'off');  % evita poluir a legenda

        yline(ax, 0, 'k-');
        xlim(ax, [0, irf_periods]);
        grid(ax, 'on'); box(ax, 'off');
        title(ax, title_list{ii}, 'Interpreter','latex');
        xlabel(ax, 'Periods', 'Interpreter','latex');
        ylabel(ax, 'PP dev from SS', 'Interpreter','latex');
        hold(ax, 'off');
    end

    % ---- Legenda dentro do painel "Real rate" (posição automática) ----
    ax_for_legend = ax_array(5);  % Real rate
    hold(ax_for_legend, 'on');
    lbls = {'Baseline ($\kappa=0.10$)', 'Flatter ($\kappa=0.05$)', 'Steeper ($\kappa=0.20$)'};
    h_dummy = gobjects(nK,1);
    for j = 1:nK
        h_dummy(j) = plot(ax_for_legend, NaN, NaN, '-', ...
                          'Color', colors_in_order(j,:), 'LineWidth', line_w);
    end
    lgd = legend(ax_for_legend, h_dummy, lbls, ...
                 'Location','best', 'Interpreter','latex', ...
                 'FontSize',8, 'Box','off');
    lgd.AutoUpdate = 'off';
    hold(ax_for_legend, 'off');

elseif experimento_a_rodar == 3
    % --- EXPERIMENTO 3: Variação de psi_piF (Figura 4) ---
    disp('Rodando Experimento 3: Variando psi_piF para o choque fiscal (Figura 4)...');

    psi_piF_values = [0.0, 0.4, 0.8];
    lambda_fixo    = 0.5;

    y_results   = zeros(irf_periods, length(psi_piF_values));
    pi_results  = zeros(irf_periods, length(psi_piF_values));
    R_results   = zeros(irf_periods, length(psi_piF_values));
    tau_results = zeros(irf_periods, length(psi_piF_values));
    r_results   = zeros(irf_periods, length(psi_piF_values));
    b_results   = zeros(irf_periods, length(psi_piF_values));

    set_param_value('lambda', lambda_fixo);

    for i = 1:length(psi_piF_values)
        set_param_value('psi_piF', psi_piF_values(i));
        stoch_simul(order=1, irf=40, nograph, noprint, irf_shocks=(eps_tau));

        y_results(:,i)   = oo_.irfs.y_eps_tau';
        pi_results(:,i)  = oo_.irfs.pi_eps_tau';
        R_results(:,i)   = oo_.irfs.R_eps_tau';
        tau_results(:,i) = oo_.irfs.tau_eps_tau';
        r_results(:,i)   = oo_.irfs.r_eps_tau';
        b_results(:,i)   = oo_.irfs.b_eps_tau';
    end

    results_cell = {y_results, pi_results, R_results, tau_results, r_results, b_results};

    % -------- Paleta: rosa (ψπ^F baixo) -> vermelho escuro (ψπ^F alto) --------
    nP = length(psi_piF_values);
    startRGB = [1.00, 0.85, 0.85];  % claro
    stopRGB  = [0.55, 0.00, 0.00];  % escuro
    cmap_sorted = [linspace(startRGB(1), stopRGB(1), nP)', ...
                   linspace(startRGB(2), stopRGB(2), nP)', ...
                   linspace(startRGB(3), stopRGB(3), nP)'];   % cores do menor p/ maior ψπ^F
    % mapeia essas cores para a ordem original de psi_piF_values
    [~, ascIdx] = sort(psi_piF_values, 'ascend');  % indices do menor p/ maior
    colors_in_order = zeros(nP,3);
    colors_in_order(ascIdx, :) = cmap_sorted;

    line_w = 1.5;

    figure('Name', 'Figure 4: Transfer shock with different Fisherian effects', 'Color', 'w');

    ax_array = gobjects(1,6);
    for ii = 1:6
        ax_array(ii) = subplot(2, 3, ii);
        ax = ax_array(ii);

        set(ax, 'ColorOrder', colors_in_order, 'NextPlot', 'replacechildren');
        hold(ax, 'on');

        % plota nas cores alinhadas à ordem original de psi_piF_values
        h_lines = plot(ax, 1:irf_periods, results_cell{ii}, 'LineWidth', line_w);
        set(h_lines, 'HandleVisibility', 'off');  % evita "data1..." na legenda

        yline(ax, 0, 'k-');
        xlim(ax, [0, irf_periods]);
        grid(ax, 'on'); box(ax, 'off');
        title(ax, title_list{ii}, 'Interpreter','latex');
        xlabel(ax, 'Periods', 'Interpreter','latex');
        ylabel(ax, 'PP dev from SS', 'Interpreter','latex');
        hold(ax, 'off');
    end

    % ---- Legenda dentro do painel "Real rate" (posição automática) ----
    ax_for_legend = ax_array(5);  % Real rate
    hold(ax_for_legend, 'on');
    h_dummy = gobjects(nP,1);
    for j = 1:nP
        h_dummy(j) = plot(ax_for_legend, NaN, NaN, '-', ...
                          'Color', colors_in_order(j,:), 'LineWidth', line_w);
    end
    lbls = arrayfun(@(v) sprintf('$\\psi_{\\pi}^{F}=%.1f$', v), psi_piF_values, 'UniformOutput', false);
    lgd = legend(ax_for_legend, h_dummy, lbls, ...
                 'Location','best', 'Interpreter','latex', ...
                 'FontSize',8, 'Box','off');
    lgd.AutoUpdate = 'off';
    hold(ax_for_legend, 'off');
end

disp('Processo concluído.');