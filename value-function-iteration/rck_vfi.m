% rck_vfi.m

clear; clc;

% ==== parâmetros (do slide) ====
gamma = 0.2;   A = 40;   alpha = 0.3;   beta = 0.9;

% ==== grade ====
Ng = 1001;                               % # pontos
xm = 0.0;
xM = A^(1/(1-alpha));                   % domínio estável
d  = (xM - xm)/(Ng-1);
x  = xm:d:xM;                           % 1xNg
N  = Ng;

% ==== matriz de payoff F(i,j) no estilo dos slides ====
F = -inf(N,N);
for i = 1:N
    for j = 1:N
        % viabilidade: c = A*x(i)^alpha - x(j) >= 0
        if A*x(i)^alpha >= x(j)
            c = A*x(i)^alpha - x(j);
            if abs(gamma - 1) < 1e-12
                F(i,j) = log(c);                       % caso limite 
            else
                % >>> forma exatamente como no slide:
                F(i,j) = (1-gamma)^(-1) * ( c^(1-gamma) );
                % (sem o "-1" constante)
            end
        else
            F(i,j) = -Inf;                              % inviável
        end
    end
end

% ==== iteração do valor (parada por tolerância) ====
v0    = zeros(N,1);
tol   = 1e-10;
maxit = 5e4;
dist  = inf; it = 0;

while dist > tol && it < maxit
    Tv = F + beta.*v0';                % aplica o operador de Bellman
    [v1, pol] = max(Tv, [], 2);        % valor e argmax
    dist = max(abs(v1 - v0));
    v0   = v1;  it = it + 1;
end
fprintf('Convergiu em %d iterações; dist = %.3e; Ng = %d\n', it, dist, Ng);

% ==== política e gráficos ====
g = x(pol)';                              % k'(k)
c = A*(x'.^alpha) - g;                    % consumo ótimo

% --- Value function (sem 45°) ---
figure;
plot(x, v1, 'LineWidth',1.6); grid on;
title('Value Function'); xlabel('k'); ylabel('V(k)');
legend('V(k)','Location','NorthWest');

% --- Optimal policy + 45° ---
figure;
plot(x, g, 'LineWidth',1.6); hold on; plot(x, x, '--');
grid on;
title('Optimal Policy'); xlabel('k'); ylabel("g(k)=k'");
legend('g(k)','45°','Location','NorthWest');

% --- Consumo + 45° ---
figure;
plot(x, c, 'LineWidth',1.6); hold on; plot(x, x, '--');
grid on;
title('Optimal Consumption'); xlabel('k'); ylabel('c(k)');
legend('c(k)','45°','Location','NorthWest');


