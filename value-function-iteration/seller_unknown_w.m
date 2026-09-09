% seller_unknown_w.m
% Seller with Unknown Demand — homogeneity reduction: w(a) = v(a,1)
% State: a in [0,1), Control: p in [a,1]
% Bellman: w(a) = max_p ((1-p)/(1-a))*p
%                   + beta*[ ((p-a)/(1-a))*w(a/p) + ((1-p)/(1-a))*w(p) ]

clear; clc;

%% ---- params ----
beta = 0.95;
Na   = 2000;                 % discretization for a \in [0,1)
Np   = 2000;                 % discretization for p \in [a,1]
eps1 = 1e-10;                % to avoid division by zero at a ~ 1

%% ---- grids ----
a      = linspace(0, 1-eps1, Na)';  % column vector, avoid a=1 exactly
p_full = linspace(0, 1, Np);        % control grid, later sliced to [a,1]

%% ---- init ----
w0  = a/(1-beta);     % good starting guess from slide
w1  = w0;
pol = a;              % placeholder for p*(a)

% anchor at a = 1 (theoretical): w(1) = 1/(1-beta)
w_at_1 = 1/(1-beta);

tol   = 1e-12;        % tighter tolerance
maxit = 5000;
dist  = Inf; it = 0;

%% ---- value function iteration (with PCHIP interpolation) ----
while dist > tol && it < maxit
    w_prev = w0;                        % keep a clean copy

    % extended grid including a=1 to avoid extrapolation artifacts
    a_ext = [a; 1];
    w_ext = [w_prev; w_at_1];

    for i = 1:Na
        ai = a(i);

        % feasible controls p in [ai, 1]
        pj = p_full(p_full >= ai);
        if isempty(pj), pj = ai; end

        den   = max(1 - ai, eps1);      % (1 - a)
        psell = (1 - pj) ./ den;        % Prob(sell)
        pnos  = 1 - psell;              % Prob(no sell) = (pj - ai)/(1-ai)

        % continuation values via PCHIP interpolation on current w_prev
        ap_over_p = ai ./ pj;                 % a' = a/p in [0,1]
        ap_over_p = min(ap_over_p, 1);        % cap at 1 to stay within a_ext
        w_ap = interp1(a_ext, w_ext, ap_over_p, 'pchip'); % w(a/p)
        w_p  = interp1(a_ext, w_ext, min(pj,1),  'pchip'); % w(p)

        curr = psell .* pj;                           % current profit if sells
        cont = beta .* ( pnos .* w_ap + psell .* w_p );
        val  = curr + cont;                           % objective at each p

        [w1(i), idx] = max(val);
        pol(i) = pj(idx);                             % argmax p*(a)
    end

    dist = max(abs(w1 - w0));
    w0   = w1;
    it   = it + 1;
end
fprintf('Converged in %d iterations; dist = %.3e\n', it, dist);

%% ---- analytics ----
a_hat = 1/(2 - beta);   % threshold in slides where g(a)=a is optimal
share_id = mean(abs(pol - a) < 1e-8);
fprintf('Share with p*(a)=a (identity on the grid): %.1f%% | a_hat = %.4f\n', ...
        100*share_id, a_hat);

%% ---- plots (com âncora em a=1 no gráfico) ----
% usa os vetores estendidos já criados acima
a_ext = [a; 1];
w_ext = [w1; w_at_1];

% Value function (mostra w(1) corretamente)
figure;
plot(a_ext, w_ext, 'LineWidth', 1.6); hold on;
plot(1, w_at_1, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4); % marca a=1
grid on; xlim([0 1]); ylim([min(w1)*0.98, 1.02*w_at_1]);
xlabel('a'); ylabel('w(a)'); title('Value Function  w(a)=v(a,1)');

% Policy vs identity
figure; plot(a, pol, 'LineWidth', 1.6); hold on;
plot(a, a, '--', 'LineWidth', 1.2);
xline(a_hat, ':', 'LineWidth', 1.2);
grid on; xlim([0 1]); ylim([0 1]);
xlabel('a'); ylabel('p^*(a)'); title('Optimal Price Policy');
legend('p^*(a)','identity p=a','a\_hat','Location','NorthWest');

% Distance to identity
figure; plot(a, pol - a, 'LineWidth', 1.6); grid on;
xlim([0 1]); xlabel('a'); ylabel('p^*(a)-a'); title('Distance to Identity');

%% ===== Heatmap of objective Q(a,p) with p*(a) overlaid =====
% What is this heatmap?
% - Each cell (a,p) shows the objective Q(a,p) = E[profit today + beta*value next]
%   computed *before* choosing p, given current belief a about the buyer’s R.
% - Colors: higher values (brighter) indicate more attractive price choices p
%   at that belief a.
% - The black line overlays the maximizing price p*(a) for each a (the policy).
%
% Mathematically:
%   Q(a,p) = ((1-p)/(1-a))*p + beta*[ ((p-a)/(1-a)) * w(a/p) + ((1-p)/(1-a)) * w(p) ]
% where:
%   (1-p)/(1-a) = Prob(sell) when R ~ U[a,1], and the posteriors are:
%     - if sell:   R ~ U[p,1] -> next state is p
%     - if not:    R ~ U[a,p] -> next state is a/p

NAh = 300;  NPh = 300;                     % heatmap resolution
a_h = linspace(0, 1-eps1, NAh);
p_h = linspace(0, 1, NPh);
Q   = nan(NAh, NPh);

% use the *final* value function with the same anchor at a=1
a_ext = [a; 1];
w_ext = [w1; w_at_1];

for ii = 1:NAh
    ai  = a_h(ii);
    den = max(1-ai, eps1);
    for jj = 1:NPh
        pj = p_h(jj);
        if pj < ai || pj > 1
            Q(ii,jj) = NaN;                          % infeasible p
        else
            psell = (1 - pj)/den;                   % Prob(sell)
            pnos  = 1 - psell;                      % Prob(no sell)
            w_ap  = interp1(a_ext, w_ext, min(ai/max(pj,eps1),1), 'pchip');
            w_p   = interp1(a_ext, w_ext, min(pj,1),                'pchip');
            curr  = psell * pj;
            cont  = beta * ( pnos * w_ap + psell * w_p );
            Q(ii,jj) = curr + cont;
        end
    end
end

% p*(a) on the same a_h grid (for overlay)
pstar_h = nan(size(a_h));
for ii = 1:NAh
    ai = a_h(ii);
    pj = p_h(p_h >= ai); if isempty(pj), pj = ai; end
    den = max(1-ai, eps1);
    psell = (1 - pj)/den;  pnos = 1 - psell;
    w_ap  = interp1(a_ext, w_ext, min(ai./pj,1), 'pchip');
    w_p   = interp1(a_ext, w_ext, min(pj,1),     'pchip');
    curr  = psell .* pj;
    cont  = beta .* ( pnos .* w_ap + psell .* w_p );
    [~, idx] = max(curr + cont);
    pstar_h(ii) = pj(idx);
end

% Heatmap; x = a, y = p
figure;
imagesc(a_h, p_h, Q');                    % transpose so rows map to y = p
set(gca,'YDir','normal'); colorbar;
hold on; plot(a_h, pstar_h, 'k-', 'LineWidth', 1.6);  % p*(a) on y-axis
xlabel('a'); ylabel('p');
title('Objective Q(a,p) and optimal p^*(a)');
legend('p^*(a)','Location','southeast');
axis tight

% =====================================================================
% HEATMAP: Objective Q(a,p) and optimal p*(a)
%
% O que está sendo plotado?
%   • Cada ponto (a,p) mostra o valor da função-objetivo
%       Q(a,p) = ((1-p)/(1-a))*p
%                + beta*[ ((p-a)/(1-a))*w(a/p) + ((1-p)/(1-a))*w(p) ].
%     Q(a,p) é o “payoff esperado” (lucro de hoje + valor futuro) dado o
%     estado a (crença sobre R ~ U[a,1]) e um preço candidato p.
%
% Eixos:
%   • Eixo x (horizontal): a ∈ [0,1) — crença/limite inferior da reserva do comprador.
%   • Eixo y (vertical):   p ∈ [a,1] — preço escolhido pelo vendedor.
%   • A região factível é o triângulo acima da diagonal p=a (abaixo dessa
%     diagonal é inviável, por isso aparece vazia/NaN).
%
% Cores:
%   • Cores mais claras (valores maiores) indicam preços p que deixam
%     Q(a,p) mais alto para aquele a — ou seja, opções mais atraentes.
%   • Cores mais escuras indicam preços com pior trade-off (vender pouco,
%     atualizar mal a crença, etc.).
%
% Linha preta (sobreposta):
%   • Curva p*(a) — o preço ótimo obtido pelo VFI para cada a.
%   • Em cada a, p*(a) coincide com o argmax da coluna correspondente de Q.
%
% Observação numérica:
%   • Para estabilizar a borda a≈1 ou p≈1, a interpolação usa a âncora
%     w(1)=1/(1-beta) e o método PCHIP (mais suave/monótono que linear).
% =====================================================================
