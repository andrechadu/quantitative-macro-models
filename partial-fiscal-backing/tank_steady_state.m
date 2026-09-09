function [bar_cR, bar_cH, bar_nR, bar_nH, bar_y, diagnostics] = ...
    tank_steady_state(beta, h, sigma, varphi, omega, ...
    bar_w, bar_YN, bar_IY, bar_g, bar_z, bar_tau)
%TANK_STEADY_STATE Compute the positive steady-state TANK allocation.
%
%   [bar_cR, bar_cH, bar_nR, bar_nH, bar_y] = ...
%       tank_steady_state(beta, h, sigma, varphi, omega, ...
%       bar_w, bar_YN, bar_IY, bar_g, bar_z, bar_tau)
%
%   solves the household-allocation and aggregate-resource block of the
%   deterministic steady state. The function is designed to be called from
%   a Dynare steady_state_model block after prices, returns, and aggregate
%   capital and investment ratios have been computed analytically.
%
%   INPUTS
%   ------
%   beta      Household discount factor.
%   h         Internal consumption-habit parameter.
%   sigma     Curvature of utility from habit-adjusted consumption.
%   varphi    Curvature of disutility from hours worked.
%   omega     Population share of hand-to-mouth households.
%   bar_w     Steady-state real wage.
%   bar_YN    Steady-state output per aggregate hour, bar_y/bar_n.
%   bar_IY    Steady-state aggregate investment-to-output ratio.
%   bar_g     Steady-state government-purchases-to-output ratio.
%   bar_z     Steady-state transfers-to-output ratio.
%   bar_tau   Steady-state revenue-to-output ratio.
%
%   OUTPUTS
%   -------
%   bar_cR    Consumption per Ricardian household.
%   bar_cH    Consumption per hand-to-mouth household.
%   bar_nR    Hours worked per Ricardian household.
%   bar_nH    Hours worked per hand-to-mouth household.
%   bar_y     Aggregate output.
%   diagnostics
%              Optional structure containing expenditure shares, aggregate
%              hours, root information, and equation residuals.
%
%   ECONOMIC SYSTEM
%   ---------------
%   The function solves the following five steady-state equations:
%
%     bar_nR^varphi =
%         (1-beta*h)*bar_w*((1-h)*bar_cR)^(-sigma)
%
%     bar_nH^varphi =
%         (1-beta*h)*bar_w*((1-h)*bar_cH)^(-sigma)
%
%     bar_cH =
%         bar_w*bar_nH + (bar_z-bar_tau)*bar_y
%
%     bar_y =
%         bar_YN*((1-omega)*bar_nR + omega*bar_nH)
%
%     bar_y =
%         (1-omega)*bar_cR + omega*bar_cH
%         + (bar_IY+bar_g)*bar_y
%
%   NUMERICAL METHOD
%   ----------------
%   For 0 < omega < 1, the five-equation system is reduced analytically
%   to one equation in xH = bar_cH/bar_y. The admissible interval is
%
%       0 < xH < (1-bar_IY-bar_g)/omega,
%
%   which guarantees positive consumption for both household types. The
%   scalar equation is solved with fzero. Logarithmic formulas are used
%   when reconstructing the allocation to improve numerical stability.
%   The representative-agent limit omega = 0 is handled separately. In
%   that case, Ricardian consumption is fixed by the aggregate resource
%   constraint, while the zero-mass HtM allocation is recovered from its
%   labor-supply and within-period budget equations.
%
%   No Optimization Toolbox is required. The function uses only standard
%   MATLAB/Octave functions.

    narginchk(11, 11);

    input_names = { ...
        'beta', 'h', 'sigma', 'varphi', 'omega', ...
        'bar_w', 'bar_YN', 'bar_IY', 'bar_g', 'bar_z', 'bar_tau'};
    input_values = { ...
        beta, h, sigma, varphi, omega, ...
        bar_w, bar_YN, bar_IY, bar_g, bar_z, bar_tau};

    for input_index = 1:numel(input_values)
        assert_finite_real_scalar(input_values{input_index}, ...
            input_names{input_index});
    end

    % Dynare's initial_estimation_checks routine multiplies every estimated
    % parameter by 1.01 when testing whether the steady state changes.  A
    % mode close to the upper boundary (for example beta=0.9907) therefore
    % reaches this function as beta=1.000607 during that diagnostic, even
    % though the actual estimated value remains inside (0,1).  Keep rejecting
    % genuinely invalid values, but project this narrow diagnostic overshoot
    % to the admissible boundary for the steady-state calculation.
    beta_for_steady_state = beta;
    dynare_check_upper_bound = 1.01 + 10*eps(1.01);
    if beta <= 0 || beta > dynare_check_upper_bound
        error('tank_steady_state:InvalidBeta', ...
            'beta must satisfy 0 < beta < 1.');
    end
    if beta >= 1
        beta_for_steady_state = 1 - 1e-10;
    end

    if h < 0 || h >= 1
        error('tank_steady_state:InvalidHabit', ...
            'h must satisfy 0 <= h < 1.');
    end

    if sigma <= 0
        error('tank_steady_state:InvalidSigma', ...
            'sigma must be strictly positive.');
    end

    if varphi <= 0
        error('tank_steady_state:InvalidVarphi', ...
            'varphi must be strictly positive.');
    end

    if omega < 0 || omega >= 1
        error('tank_steady_state:InvalidOmega', ...
            ['omega must satisfy 0 <= omega < 1. The all-HtM boundary ', ...
             'omega = 1 is not supported because Ricardian asset-pricing ', ...
             'conditions remain part of the model.']);
    end

    if bar_w <= 0
        error('tank_steady_state:InvalidWage', ...
            'bar_w must be strictly positive.');
    end

    if bar_YN <= 0
        error('tank_steady_state:InvalidProductivity', ...
            'bar_YN must be strictly positive.');
    end

    if bar_IY < 0
        error('tank_steady_state:InvalidInvestmentShare', ...
            'bar_IY must be nonnegative.');
    end

    if bar_g < 0
        error('tank_steady_state:InvalidGovernmentShare', ...
            'bar_g must be nonnegative.');
    end

    private_consumption_share = 1 - bar_IY - bar_g;

    if private_consumption_share <= 0
        error('tank_steady_state:NoPrivateConsumption', ...
            ['The resource constraint requires ', ...
             '1-bar_IY-bar_g > 0.']);
    end

    consumption_hours_elasticity = sigma / varphi;
    net_transfer_share = bar_z - bar_tau;

    % From labor supply:
    % bar_nJ = labor_scale*bar_cJ^(-sigma/varphi), J in {R,H}.
    log_labor_scale = ...
        (log(1-beta_for_steady_state*h) + log(bar_w) ...
         - sigma*log(1-h)) / varphi;
    labor_scale = exp(log_labor_scale);

    if ~isfinite(labor_scale) || labor_scale <= 0
        error('tank_steady_state:InvalidLaborScale', ...
            ['The implied labor-supply scale is nonpositive or ', ...
             'nonfinite. Check beta, h, sigma, varphi, and bar_w.']);
    end

    solver_options = optimset( ...
        'Display', 'off', ...
        'TolX', 1e-12);

    if omega == 0
        % Representative-agent limit. The HtM household has zero measure,
        % but its individual allocation remains well defined by its two
        % intratemporal conditions.
        ricardian_consumption_share = private_consumption_share;
        residual_function = @(log_share_h) ...
            representative_agent_htm_residual( ...
                log_share_h, ...
                ricardian_consumption_share, ...
                consumption_hours_elasticity, ...
                bar_w, ...
                bar_YN, ...
                net_transfer_share);

        [transformed_root, root_residual, exit_flag] = ...
            fzero(residual_function, [-50, 50], solver_options);

        if exit_flag <= 0
            error('tank_steady_state:RootSolverFailure', ...
                ['fzero did not converge for the zero-mass HtM ', ...
                 'allocation in the representative-agent limit. ', ...
                 'Exit flag: %d.'], exit_flag);
        end

        hand_to_mouth_consumption_share = exp(transformed_root);
        candidate_residuals = root_residual;
        candidate_exit_flags = exit_flag;
    else
        % The upper bound follows from positive Ricardian consumption:
        % xR = (private_consumption_share-omega*xH)/(1-omega) > 0.
        hand_to_mouth_share_upper_bound = ...
            private_consumption_share / omega;

        residual_function = @(transformed_share) ...
            hand_to_mouth_share_residual( ...
                transformed_share, ...
                hand_to_mouth_share_upper_bound, ...
                private_consumption_share, ...
                omega, ...
                consumption_hours_elasticity, ...
                bar_w, ...
                bar_YN, ...
                net_transfer_share);

        % Search for sign-changing brackets in a logistic transformation of
        % the full admissible interval. This avoids evaluating the residual
        % at zero consumption or at zero Ricardian consumption.
        transformed_grid = linspace(-50, 50, 2001);
        residual_grid = arrayfun(residual_function, transformed_grid);

        valid_grid = isfinite(residual_grid);
        bracket_indices = find( ...
            valid_grid(1:end-1) & valid_grid(2:end) & ...
            sign(residual_grid(1:end-1)) .* ...
            sign(residual_grid(2:end)) <= 0);

        if isempty(bracket_indices)
            finite_residuals = residual_grid(valid_grid);

            if isempty(finite_residuals)
                residual_range = [NaN, NaN];
            else
                residual_range = ...
                    [min(finite_residuals), max(finite_residuals)];
            end

            error('tank_steady_state:RootNotBracketed', ...
                ['No root was bracketed inside the economically admissible ', ...
                 'consumption-share interval. The scanned residual range ', ...
                 'was [%g, %g]. Check the fiscal targets and expenditure ', ...
                 'shares.'], residual_range(1), residual_range(2));
        end

        candidate_shares = zeros(numel(bracket_indices), 1);
        candidate_residuals = zeros(numel(bracket_indices), 1);
        candidate_exit_flags = zeros(numel(bracket_indices), 1);

        for root_index = 1:numel(bracket_indices)
            grid_index = bracket_indices(root_index);
            bracket = transformed_grid(grid_index:grid_index+1);

            [transformed_root, root_residual, exit_flag] = ...
                fzero(residual_function, bracket, solver_options);

            if exit_flag <= 0
                error('tank_steady_state:RootSolverFailure', ...
                    ['fzero did not converge inside a valid sign-changing ', ...
                     'bracket. Exit flag: %d.'], exit_flag);
            end

            candidate_shares(root_index) = logistic_share( ...
                transformed_root, hand_to_mouth_share_upper_bound);
            candidate_residuals(root_index) = root_residual;
            candidate_exit_flags(root_index) = exit_flag;
        end

        % Adjacent grid intervals can occasionally identify the same root.
        candidate_shares = sort(candidate_shares);
        root_separation_tolerance = ...
            1e-9 * max(1, hand_to_mouth_share_upper_bound);
        distinct_root = [true; ...
            abs(diff(candidate_shares)) > root_separation_tolerance];
        candidate_shares = candidate_shares(distinct_root);

        if numel(candidate_shares) > 1
            error('tank_steady_state:MultipleSteadyStates', ...
                ['The scalar allocation equation has %d distinct positive ', ...
                 'roots. The calibration does not select a unique steady ', ...
                 'state.'], numel(candidate_shares));
        end

        hand_to_mouth_consumption_share = candidate_shares(1);
        ricardian_consumption_share = ...
            (private_consumption_share ...
            - omega*hand_to_mouth_consumption_share) ...
            / (1-omega);
    end

    if hand_to_mouth_consumption_share <= 0 || ...
            ricardian_consumption_share <= 0
        error('tank_steady_state:NonpositiveConsumptionShare', ...
            'The computed consumption share is not strictly positive.');
    end

    log_labor_mix = log_weighted_labor_mix( ...
        ricardian_consumption_share, ...
        hand_to_mouth_consumption_share, ...
        omega, ...
        consumption_hours_elasticity);

    log_bar_y = ...
        (log(bar_YN) + log_labor_scale + log_labor_mix) ...
        / (1 + consumption_hours_elasticity);
    bar_y = exp(log_bar_y);

    bar_cR = ricardian_consumption_share * bar_y;
    bar_cH = hand_to_mouth_consumption_share * bar_y;

    bar_nR = exp( ...
        log_labor_scale ...
        - consumption_hours_elasticity*log(bar_cR));
    bar_nH = exp( ...
        log_labor_scale ...
        - consumption_hours_elasticity*log(bar_cH));

    solution_vector = [bar_cR, bar_cH, bar_nR, bar_nH, bar_y];

    if any(~isfinite(solution_vector)) || any(solution_vector <= 0)
        error('tank_steady_state:InvalidSolution', ...
            ['The reconstructed steady state contains a nonpositive or ', ...
             'nonfinite allocation.']);
    end

    aggregate_hours = ...
        (1-omega)*bar_nR + omega*bar_nH;

    labor_rhs_R = ...
        (1-beta_for_steady_state*h)*bar_w*((1-h)*bar_cR)^(-sigma);
    labor_rhs_H = ...
        (1-beta_for_steady_state*h)*bar_w*((1-h)*bar_cH)^(-sigma);

    htm_budget_rhs = ...
        bar_w*bar_nH + net_transfer_share*bar_y;
    production_rhs = ...
        bar_YN*aggregate_hours;
    resource_rhs = ...
        (1-omega)*bar_cR ...
        + omega*bar_cH ...
        + (bar_IY+bar_g)*bar_y;

    equation_residuals = [ ...
        scaled_difference(bar_nR^varphi, labor_rhs_R); ...
        scaled_difference(bar_nH^varphi, labor_rhs_H); ...
        scaled_difference(bar_cH, htm_budget_rhs); ...
        scaled_difference(bar_y, production_rhs); ...
        scaled_difference(bar_y, resource_rhs)];

    residual_tolerance = 1e-9;
    maximum_scaled_residual = max(abs(equation_residuals));

    if maximum_scaled_residual > residual_tolerance
        error('tank_steady_state:ResidualCheckFailed', ...
            ['The candidate steady state fails the equation check. ', ...
             'Maximum scaled residual: %.3e.'], ...
            maximum_scaled_residual);
    end

    diagnostics = struct();
    diagnostics.private_consumption_share = private_consumption_share;
    diagnostics.ricardian_consumption_share = ...
        ricardian_consumption_share;
    diagnostics.hand_to_mouth_consumption_share = ...
        hand_to_mouth_consumption_share;
    diagnostics.aggregate_hours = aggregate_hours;
    diagnostics.net_transfer_share = net_transfer_share;
    diagnostics.equation_residuals = equation_residuals;
    diagnostics.maximum_scaled_residual = maximum_scaled_residual;
    diagnostics.candidate_root_residuals = candidate_residuals;
    diagnostics.candidate_root_exit_flags = candidate_exit_flags;
end


function residual = hand_to_mouth_share_residual( ...
    transformed_share, share_upper_bound, private_consumption_share, ...
    omega, consumption_hours_elasticity, bar_w, bar_YN, ...
    net_transfer_share)
%HAND_TO_MOUTH_SHARE_RESIDUAL Evaluate the scalar steady-state equation.

    share_h = logistic_share(transformed_share, share_upper_bound);
    share_r = ...
        (private_consumption_share-omega*share_h)/(1-omega);

    if share_h <= 0 || share_r <= 0
        residual = NaN;
        return;
    end

    log_labor_mix = log_weighted_labor_mix( ...
        share_r, share_h, omega, consumption_hours_elasticity);

    log_hand_to_mouth_labor_term = ...
        -consumption_hours_elasticity*log(share_h);

    normalized_hand_to_mouth_labor = exp( ...
        log_hand_to_mouth_labor_term-log_labor_mix);

    residual = ...
        share_h ...
        - net_transfer_share ...
        - (bar_w/bar_YN)*normalized_hand_to_mouth_labor;
end


function residual = representative_agent_htm_residual( ...
    log_share_h, share_r, consumption_hours_elasticity, ...
    bar_w, bar_YN, net_transfer_share)
%REPRESENTATIVE_AGENT_HTM_RESIDUAL Recover the zero-mass HtM allocation.

    share_h = exp(log_share_h);
    relative_labor = ...
        exp(consumption_hours_elasticity*(log(share_r)-log(share_h)));

    residual = ...
        share_h ...
        - net_transfer_share ...
        - (bar_w/bar_YN)*relative_labor;
end


function share = logistic_share(transformed_share, share_upper_bound)
%LOGISTIC_SHARE Map the real line to the open admissible share interval.

    if transformed_share >= 0
        share = share_upper_bound / (1+exp(-transformed_share));
    else
        exponential_term = exp(transformed_share);
        share = ...
            share_upper_bound*exponential_term/(1+exponential_term);
    end
end


function log_mix = log_weighted_labor_mix( ...
    share_r, share_h, omega, consumption_hours_elasticity)
%LOG_WEIGHTED_LABOR_MIX Compute the labor aggregator in logarithms.

    log_ricardian_term = ...
        log(1-omega) ...
        - consumption_hours_elasticity*log(share_r);
    log_hand_to_mouth_term = ...
        log(omega) ...
        - consumption_hours_elasticity*log(share_h);

    largest_log_term = max( ...
        log_ricardian_term, log_hand_to_mouth_term);

    log_mix = ...
        largest_log_term ...
        + log( ...
            exp(log_ricardian_term-largest_log_term) ...
            + exp(log_hand_to_mouth_term-largest_log_term));
end


function residual = scaled_difference(left_hand_side, right_hand_side)
%SCALED_DIFFERENCE Return a dimensionless equation residual.

    residual_scale = max( ...
        [1, abs(left_hand_side), abs(right_hand_side)]);
    residual = ...
        (left_hand_side-right_hand_side)/residual_scale;
end


function assert_finite_real_scalar(value, value_name)
%ASSERT_FINITE_REAL_SCALAR Validate a scalar function argument.

    if ~isnumeric(value) || ~isscalar(value) || ...
            ~isreal(value) || ~isfinite(value)
        error('tank_steady_state:InvalidInput', ...
            '%s must be a finite real numeric scalar.', value_name);
    end
end
