% human_capital_vfi.m
% Growth with Human Capital Accumulation (slides)
% u(c) = c^sigma / sigma,  f(L) = L^alpha
% FE: v(k) = max_{y in Γ(k)} { u( f( k * phi(y/k) ) ) + beta * v(y) }
% Γ(k) = [ (1-delta)k , (1+lambda)k ]

clear; clc;

%% ===== Parâmetros =====
alpha = 0.30;          % produção f(L) = L^alpha
sigma = 0.50;          % u(c) = c^sigma / sigma   (0<sigma<1)
beta  = 0.9;
delta = 0.08;          % depreciação
lambda = 0.05;         % limite de crescimento

p = 3; rL = 1 - delta; rU = 1 + lambda;
phi = @(r) max(1e-12, (rU.^p - r.^p) ./ (rU.^p - rL.^p));


%% ===== Grade de k =====
Ng = 1001;                                  % # pontos da grade de estado
kmin = 1e-6;                               % evitar zero estrito
kmax = 200;                                 % limite superior numérico
k    = linspace(kmin, kmax, Ng);           % 1xNg
N    = Ng;

%% ===== Matriz de payoff F(k_i, k'_j) com Γ(k) = [(1-δ)k, (1+λ)k] =====
F = -inf(N,N);
for i = 1:N
    ki   = k(i);
    ymin = (1 - delta) * ki;
    ymax = (1 + lambda) * ki;

    % projeta os limites contínuos para índices da grade
    jmin = find(k >= ymin, 1, 'first');   if isempty(jmin), jmin = N; end
    jmax = find(k <= ymax, 1, 'last');    if isempty(jmax), jmax = 1; end
    if jmin > jmax
        % se o intervalo é mais estreito que o passo, pega o ponto mais próximo
        [~, jstar] = min(abs(k - 0.5*(ymin + ymax)));
        jmin = jstar; jmax = jstar;
    end

    for j = jmin:jmax
        yj = k(j);
        r  = yj / ki;                 % razão y/k
        L  = ki * phi(r);             % trabalho efetivo
        c  = L^alpha;                 % toda produção é consumo

        if c <= 0 || ~isfinite(c)
            F(i,j) = -Inf;
        else
            if abs(sigma) < 1e-12
                F(i,j) = log(c);
            else
                F(i,j) = (c^sigma)/sigma;
            end
        end
    end
end


%% ===== Iteração do valor (parada por tolerância) =====
v0    = zeros(N,1);
tol   = 1e-10;
maxit = 5e4;
dist  = inf; it = 0;

while dist > tol && it < maxit
    Tv = F + beta .* v0';                    % Bellman
    [v1, pol] = max(Tv, [], 2);              % valor e argmax
    dist = max(abs(v1 - v0));
    v0   = v1;
    it   = it + 1;
end
fprintf('Convergiu em %d iterações; dist = %.3e\n', it, dist);

%% ===== Política ótima, consumo e gráficos =====
g = k(pol)';                                 % g(k)=k'
% Consumo (= produção): c(k) = [k * phi(g(k)/k)]^alpha
c = ( k'.*phi(g./k') ).^alpha;

% V(k)
figure; plot(k, v1, 'LineWidth',1.6); grid on;
title('Value Function'); xlabel('k'); ylabel('V(k)');

% g(k) com 45°
figure; plot(k, g, 'LineWidth',1.6); hold on; plot(k, k, '--');
grid on; title('Optimal Policy'); xlabel('k'); ylabel("g(k)=k'");
legend('g(k)','45°','Location','NorthWest');

% c(k) com 45°
figure; plot(k, c, 'LineWidth',1.6); hold on; plot(k, k, '--');
grid on; title('Optimal Consumption'); xlabel('k'); ylabel('c(k)');
legend('c(k)','45°','Location','NorthWest');

% --- Optimal policy com banda viável ---
figure;
plot(k, g, 'LineWidth',1.6); hold on;
plot(k, (1-delta)*k, ':', 'LineWidth',1.2);
plot(k, (1+lambda)*k, ':', 'LineWidth',1.2);
grid on;
title('Optimal Policy vs. Feasible Set');
xlabel('k'); ylabel("g(k)=k'");
legend('g(k)','(1-\delta)k','(1+\lambda)k','Location','NorthWest');

% --- Ratio r*(k) ---
rstar = g ./ k';
figure;
plot(k, rstar, 'LineWidth',1.6); hold on;
yline(1-delta, ':'); yline(1+lambda, ':');
grid on; title('r^*(k) = g(k)/k');
xlabel('k'); ylabel('r^*(k)');
legend('r^*(k)','1-\delta','1+\lambda','Location','SouthEast');

% diagnóstico “canto”
share_upper = mean(abs(g - (1+lambda)*k') < 1e-8);
share_lower = mean(abs(g - (1-delta)*k') < 1e-8);
fprintf('Pct no limite superior: %.1f%% | inferior: %.1f%%\n', ...
        100*share_upper, 100*share_lower);

% -------------------------------------------------------------------------
% INTERPRETAÇÃO DOS GRÁFICOS (r*(k) e Política vs. Conjunto Viável)
%
% 1) Gráfico r*(k) = g(k)/k
%    - A curva azul mostra o ratio ótimo r*(k) = k' / k escolhido pelo VFI.
%    - Os "degraus" (pequenas quebras) vêm da discretização: g(k) é escolhido
%      entre pontos de grade para k'. Quando o índice ótimo j muda de um
%      ponto para outro, r*(k)=k'(k)/k sofre saltos discretos.
%
% 2) Gráfico "Optimal Policy vs. Feasible Set"
%    - Azul: g(k) = k' ótimo. Pontilhado vermelho: (1-delta)k (limite inferior).
%      Pontilhado amarelo: (1+lambda)k (limite superior do conjunto viável).
%    - Observação: g(k) fica quase colado em (1+lambda)k, confirmando que a
%      solução ótima está próxima do limite superior da banda viável para
%      a maior parte da grade (coerente com o r*(k) acima).
%    - A comparação correta para esta política é com as linhas da banda
%      viável, não com a linha de 45°. A 45° não é a restrição relevante
%      aqui: o modelo impõe k' ∈ [(1-delta)k, (1+lambda)k].