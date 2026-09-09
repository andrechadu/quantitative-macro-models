// ============================================================================
// TANK NEW KEYNESIAN MODEL WITH PARTIAL FISCAL BACKING
// ============================================================================
//
// Main features:
// - Ricardian and hand-to-mouth households;
// - internal consumption habits;
// - capital owned only by Ricardian households;
// - investment adjustment costs;
// - monopolistic competition and Rotemberg price adjustment;
// - long-duration nominal government debt;
// - separate rules for revenue, transfers, and government purchases;
// - a fiscal-led shadow economy implementing partial fiscal backing.
//
// Mechanism version 08: ten-observable specification with three iid measurement errors:
// - default Dynare timing for capital, without a special timing declaration;
// - SAMBA-compatible parameters calibrated or assigned SAMBA priors;
// - one marginal-efficiency-of-investment state in each economy, driven by
//   the same structural innovation e_i;
// - seven economic innovations and exactly three measurement errors;
// - no semistructural term-premium state or innovation;
// - b_target enters both debt laws and the debt observable directly.
//
// All nonfiscal endogenous variables in model(linear) are log deviations
// from the deterministic steady state. Fiscal ratios and debt ratios are
// level deviations from their steady-state values.
//
// Capital uses Dynare's default timing. Investment at t determines end-of-
// period capital k, while production at t uses installed capital k(-1).
// ============================================================================




// ============================================================================
// CALIBRATED MECHANISM FILE
// ============================================================================
//
// This file reproduces the economic structure of tank_pfb_mechanism_v08.mod and
// generates calibrated IRF grids only. It contains no estimation or MCMC.
// Composite grids are saved as PNG files in the three mechanism folders.
//

// ============================================================================
// ENDOGENOUS VARIABLES: MAIN ECONOMY
// ============================================================================

var
    cR              // Ricardian consumption
    cH              // Hand-to-mouth consumption
    nR              // Ricardian hours
    nH              // Hand-to-mouth hours
    lamR            // Ricardian marginal value of wealth
    lamH            // Hand-to-mouth marginal value of wealth
    n               // Aggregate hours
    y               // Aggregate output
    inv             // Aggregate investment
    k               // End-of-period capital stock
    q               // Value of installed capital
    muI             // Main marginal-efficiency-of-investment state
    rK              // Real rental rate of capital
    w               // Real wage
    mc              // Real marginal cost
    pi              // Gross inflation
    R               // Gross short-term nominal interest rate
    pB              // Price of the long-duration government bond
    Rb              // Realized long-bond return from t-1 to t
    b               // Debt-to-output ratio: level deviation
    tau             // Revenue-to-output ratio: level deviation
    z               // Transfer-to-output ratio: level deviation
    g               // Government-purchases-to-output ratio: level deviation
    a               // Total factor productivity
    xi              // Consumption-preference shifter
    ystar           // Technology-based output target
    r_real          // Ex-ante real short-term interest rate
    primary_surplus // Revenue net of transfers and purchases
    mp_innov        // Monetary-policy MA(1) innovation process. // (SAMBA, 2023)

// ============================================================================
// EMPIRICAL AUXILIARIES AND OBSERVABLES
// ============================================================================

    c_model         // Aggregate-consumption log deviation
    y1_model        // Log deviation of the geometric-bond yield
    g_ratio_model   // Government purchases in p.p. of quarterly GDP
    z_ratio_model   // Transfers in p.p. of quarterly GDP
    pb_ratio_model  // Primary surplus in p.p. of quarterly GDP
    me_y            // Iid measurement error for output growth
    me_y1           // Iid measurement error for the one-year-yield change
    me_b            // Iid measurement error for gross debt

    y_obs           // Real GDP growth, quarterly percent
    c_obs           // Real household-consumption growth, quarterly percent
    inv_obs         // Real investment growth, quarterly percent
    pi_obs          // IPCA inflation, quarterly percent
    dR_obs          // Change in effective quarterly Selic, percentage points
    dy1_obs         // Change in effective quarterly one-year yield, p.p.
    dg_obs          // Change in purchases/GDP, percentage points
    dz_obs          // Change in transfers/GDP, percentage points
    dpb_obs         // Change in primary-balance/GDP, percentage points
    b_obs           // Gross debt, percent of annual GDP

// ============================================================================
// ENDOGENOUS VARIABLES: FISCAL-LED SHADOW ECONOMY
// ============================================================================

    cRF             // Shadow Ricardian consumption
    cHF             // Shadow hand-to-mouth consumption
    nRF             // Shadow Ricardian hours
    nHF             // Shadow hand-to-mouth hours
    lamRF           // Shadow Ricardian marginal value of wealth
    lamHF           // Shadow hand-to-mouth marginal value of wealth
    nF              // Shadow aggregate hours
    yF              // Shadow aggregate output
    invF            // Shadow aggregate investment
    kF              // Shadow end-of-period capital stock
    qF              // Shadow value of installed capital
    muIF            // Shadow MEI state, with the same innovation as muI
    rKF             // Shadow real rental rate of capital
    wF              // Shadow real wage
    mcF             // Shadow real marginal cost
    piF             // Fiscal-led shadow inflation
    RF              // Shadow short-term nominal interest rate
    pBF             // Shadow long-duration bond price
    RbF             // Shadow realized long-bond return
    bF              // Shadow unfunded-debt ratio: level deviation
    tauF            // Shadow revenue ratio: level deviation
    zF              // Shadow transfer ratio: level deviation
    gF              // Shadow purchases ratio: level deviation
    aF              // Shadow total factor productivity
    xiF             // Shadow preference shifter
    yFstar          // Shadow technology-based output target
;

// ============================================================================
// EXOGENOUS INNOVATIONS
// ============================================================================

varexo
    e_a             // Technology innovation
    e_d             // Preference or demand innovation
    e_i             // Marginal-efficiency-of-investment innovation
    e_tau           // Government-revenue innovation
    e_z             // Government-transfer innovation
    e_g             // Government-purchases innovation
    e_mp            // Monetary-policy innovation
    e_me_y          // Measurement innovation for output growth
    e_me_y1         // Measurement innovation for the one-year-yield change
    e_me_b          // Measurement innovation for gross debt
;

// ============================================================================
// PARAMETERS
// ============================================================================

parameters
    beta            // Discount factor
    h               // Internal consumption-habit parameter
    sigma           // Consumption curvature and inverse IES
    varphi          // Labor curvature and inverse Frisch elasticity
    omega           // Hand-to-mouth population share

    delta           // Capital depreciation rate
    phi_i           // Investment-adjustment-cost curvature
    alpha           // Capital elasticity in production
    eps_p           // Elasticity of substitution across goods
    kappa_p         // Rotemberg price-adjustment parameter

    rho_m           // Long-bond coupon-decay parameter
    lambda          // Fiscal-backing parameter

    A_ss            // Steady-state technology normalization
    pi_target       // Gross inflation target
    b_target        // Debt-to-quarterly-output target
    g_target        // Government-purchases-to-output target
    z_target        // Transfer-to-output target
    gamma_obs       // Common quarterly trend growth in percent
    pi_ss_pct       // Steady-state quarterly inflation in percent

    bar_pi          // Steady-state gross inflation
    bar_R           // Steady-state gross nominal interest rate
    bar_b           // Steady-state debt-to-output ratio
    bar_g           // Steady-state purchases-to-output ratio
    bar_z           // Steady-state transfer-to-output ratio
    bar_tau         // Steady-state revenue-to-output ratio
    bar_mc          // Steady-state real marginal cost
    bar_rK          // Steady-state real rental rate of capital
    bar_KY          // Capital-to-output ratio
    bar_IY          // Investment-to-output ratio
    bar_YN          // Output per aggregate hour
    bar_w           // Steady-state real wage
    bar_cR          // Ricardian consumption per household
    bar_cH          // Hand-to-mouth consumption per household
    bar_nR          // Ricardian hours per household
    bar_nH          // Hand-to-mouth hours per household
    bar_y           // Steady-state aggregate output
    bar_iR          // Investment per Ricardian household

    rho_a           // Persistence of technology
    rho_d           // Persistence of preference shifter
    rho_i           // Persistence of marginal investment efficiency

    rho_tau         // Persistence of revenue
    rho_z           // Persistence of transfers
    rho_g           // Persistence of government purchases
    db_tau          // Revenue response to backed debt
    db_z            // Signed transfer response to backed debt
    db_g            // Signed purchases response to backed debt
    dy              // Revenue response to lagged output
    dDy             // Revenue response to output growth
    dDn             // Transfer-response magnitude to employment growth

    rho_R           // Interest-rate smoothing
    psi_pi          // Response to monetary-led inflation #Psi_pi = 1 + psi_pi
    psi_piF         // Response to fiscal-led inflation
    psi_y           // Response to the relative output gap
    psi_Dy          // Response to relative output-gap growth
    mu_mp           // MA(1) coefficient for monetary-policy innovation.        // (SAMBA, 2023)

    sigma_me_y      // Fraction of the output-growth measurement-error cap
    sigma_me_y1     // Fraction of the yield-change measurement-error cap
    sigma_me_b      // Debt measurement-error standard deviation, annual-GDP p.p.
    me_cap_y        // Output-growth measurement-error volatility cap
    me_cap_y1       // Yield-change measurement-error volatility cap
;

// ============================================================================
// BASELINE QUARTERLY CALIBRATION
// ============================================================================

A_ss = 1.0;

// Five-percent annual neutral real interest rate.
beta = (1+0.05)^(-1/4);

// SAMBA annual inflation target; the estimable object is the
// nonannualized quarterly log rate pi_ss_pct.
pi_ss_pct = 100*log(1.045)/4; // (SAMBA, 2023)
pi_target = exp(pi_ss_pct/100);

// Common balanced-growth intercept used by y_obs, c_obs, and inv_obs.
// Partial analogue: SAMBA's productivity trend centers the broader TANK-PFB observable trend.
gamma_obs = 100*log(1.02)/4; // (SAMBA, 2023)

// The cap-fraction construction follows the measurement-error framework in
// SAMBA (2023). The caps are the sample standard deviations in the observable
// units used below.
me_cap_y = 0.578773517817218; // (SAMBA, 2023)
me_cap_y1 = 0.370026954058522; // Adapted from the SAMBA (2023) cap-fraction framework.

mu_mp = 0.421; // (SAMBA, 2023)
sigma_me_y = 0.50; // Adapted from the SAMBA (2023) cap-fraction framework.
sigma_me_y1 = 0.50; // Adapted from the SAMBA (2023) cap-fraction framework.
sigma_me_b = 1.00; // Smets and Wouters (2024), debt measurement-error prior mean.

// Gross debt equal to 80 percent of annual GDP. Multiplication by four
// converts the ratio to units of quarterly GDP.
b_target = 4*0.80;

g_target = 0.171;
z_target = 1755035/10856112;

// Average federal-debt maturity of 4.08 years.
rho_m = 1-1/(4*4.08);

// Household preferences.
h       = 0.784; // (SAMBA, 2023)
sigma   = 1.0001; // (SAMBA, 2023)
varphi  = 2.0;
omega   = 0.20; // (SAMBA, 2023)

// Capital, investment, and production.
delta   = 0.027; // (SAMBA, 2023)
phi_i   = 11.88; // (SAMBA, 2023)
eps_p   = 11.0; // (SAMBA, 2023)

// Capital elasticity chosen to match I/Y = 0.183 while accounting for
// the steady-state monopolistic-competition markup.
bar_IY = 0.183; // (SAMBA, 2023)
alpha = bar_IY*(1/beta-(1-delta))/(delta*((eps_p-1)/eps_p)); // (SAMBA, 2023)

// The target reduced-form NKPC slope is 0.0075.
kappa_p = (eps_p-1)/(0.0075*pi_target^2); // (SAMBA, 2023)

// Monetary policy.
rho_R   = 0.804; // (SAMBA, 2023)
psi_pi  = 1.211; // (SAMBA, 2023) #Psi_pi = 1 + psi_pi
psi_y   = 0.0;
psi_Dy  = 0.120; // (SAMBA, 2023)
psi_piF = 0.0;

// Partial fiscal backing.
lambda  = 0.75; //0.642;

// Fiscal-flow persistence.
// Adaptation of SAMBA lump-sum-tax persistence to aggregate revenue.
rho_tau = 0.405; // (SAMBA, 2023)
rho_z   = 0.99;
rho_g   = 0.262; // (SAMBA, 2023)

// Signed fiscal coefficients: revenue rises with backed debt, whereas
// transfers and purchases fall.
db_tau  = 0.01;
db_z    = -0.07;
db_g    = -0.03;
dy      = 0.49;
dDy     = 0.39;
dDn     = 0.18;

// Structural-shock persistence.
rho_a   = 0.428; // (SAMBA, 2023)
rho_d   = 0.581; // (SAMBA, 2023)
rho_i   = 0.274; // (SAMBA, 2023)

// ============================================================================
// LINEAR EQUILIBRIUM SYSTEM
// ============================================================================

model(linear);

    // Endogenous steady-state expenditure and employment shares.
    # s_cR = (1-omega)*bar_cR/bar_y;
    # s_cH = omega*bar_cH/bar_y;
    # s_i  = (1-omega)*bar_iR/bar_y;
    # s_g  = bar_g;
    # s_c  = s_cR+s_cH;
    # s_nR = (1-omega)*bar_nR
             /((1-omega)*bar_nR+omega*bar_nH);
    # s_nH = omega*bar_nH
             /((1-omega)*bar_nR+omega*bar_nH);
    # Psi_pi = 1 + psi_pi;

    // ========================================================================
    // MAIN ECONOMY
    // ========================================================================

    // Marginal values of real wealth with internal consumption habits.
    lamR =
        (xi-beta*h*xi(+1))/(1-beta*h)
        -sigma/((1-h)*(1-beta*h))
        *((1+beta*h^2)*cR-h*cR(-1)-beta*h*cR(+1));

    lamH =
        (xi-beta*h*xi(+1))/(1-beta*h)
        -sigma/((1-h)*(1-beta*h))
        *((1+beta*h^2)*cH-h*cH(-1)-beta*h*cH(+1));

    // Labor supply and the Ricardian Euler equation.
    varphi*nR = lamR+w;
    varphi*nH = lamH+w;
    lamR-lamR(+1) = R-pi(+1);

    // Hand-to-mouth budget constraint, aggregation, and market clearing.
    bar_cH*cH =
        bar_w*bar_nH*(w+nH)
        +bar_y*(z-tau+(bar_z-bar_tau)*y);

    n = s_nR*nR+s_nH*nH;

    (1-s_g)*y =
        s_cR*cR+s_cH*cH+s_i*inv+g;

    // Default Dynare timing: investment at t determines end-of-period capital
    // k, and k(-1) is the installed stock used in production at t.
    k = (1-delta)*k(-1)+delta*(inv+muI);

    // Investment optimality and the value of installed capital.
    q+muI = phi_i*((1+beta)*inv-inv(-1)-beta*inv(+1));

    q =
        beta*(1-delta)*q(+1)
        +beta*bar_rK*rK(+1)
        -(R-pi(+1));

    // Production and factor demand.
    y  = a+alpha*k(-1)+(1-alpha)*n;
    rK = mc+y-k(-1);
    w  = mc+y-n;

    // Rotemberg New Keynesian Phillips curve.
    pi = beta*pi(+1)+(eps_p-1)/(kappa_p*bar_pi^2)*mc;

    // Long-duration government bonds.
    Rb = (rho_m/bar_R)*pB-pB(-1);
    R  = Rb(+1);

    // Government budget constraint.
    b =
        (1/beta)*b(-1)
        +(b_target/beta)*(Rb-pi+y(-1)-y)
        -tau+z+g;

    // Reporting identity for the primary surplus.
    primary_surplus = tau-z-g;

    // Fiscal reaction functions. The debt coefficients use the signed
    // convention adopted in the paper.
    tau =
        rho_tau*tau(-1)
        +(1-rho_tau)*(db_tau*(b(-1)-bF(-1))+dy*y(-1))
        +dDy*(y-y(-1))
        +e_tau;

    z =
        rho_z*z(-1)
        +(1-rho_z)*db_z*(b(-1)-bF(-1))
        -dDn*(n-n(-1))
        +e_z;

    g =
        rho_g*g(-1)
        +(1-rho_g)*db_g*(b(-1)-bF(-1))
        +e_g;

    // Ex-ante real rate and monetary-policy rule.
    r_real = R-pi(+1);

    // The paper establishes an MA(1) monetary shock but does not print its sign.
    // The negative sign implements the short-run mean-reverting convention and
    // is recorded as an inference rather than an exact quoted SAMBA equation.
    mp_innov = e_mp-mu_mp*e_mp(-1); // (SAMBA, 2023)

    R =
        rho_R*R(-1)
        +(1-rho_R)
         *(Psi_pi*(pi-piF)
           +psi_piF*piF
           +psi_y*((y-ystar)-(yF-yFstar)))
        +psi_Dy
         *(((y-ystar)-(y(-1)-ystar(-1)))
           -((yF-yFstar)-(yF(-1)-yFstar(-1))))
        +mp_innov; // (SAMBA, 2023)

    // Structural shock processes.
    a  = rho_a*a(-1)+e_a; // (SAMBA, 2023)
    xi = rho_d*xi(-1)+e_d; // (SAMBA, 2023)
    muI = rho_i*muI(-1)+e_i; // (SAMBA, 2023)

    // Technology-based output target.
    ystar = a;

    // ========================================================================
    // FISCAL-LED SHADOW ECONOMY
    // ========================================================================

    // Shadow marginal values of real wealth.
    lamRF =
        (xiF-beta*h*xiF(+1))/(1-beta*h)
        -sigma/((1-h)*(1-beta*h))
        *((1+beta*h^2)*cRF-h*cRF(-1)-beta*h*cRF(+1));

    lamHF =
        (xiF-beta*h*xiF(+1))/(1-beta*h)
        -sigma/((1-h)*(1-beta*h))
        *((1+beta*h^2)*cHF-h*cHF(-1)-beta*h*cHF(+1));

    // Shadow labor supply and the Ricardian Euler equation.
    varphi*nRF = lamRF+wF;
    varphi*nHF = lamHF+wF;
    lamRF-lamRF(+1) = RF-piF(+1);

    // Shadow hand-to-mouth budget, aggregation, and market clearing.
    bar_cH*cHF =
        bar_w*bar_nH*(wF+nHF)
        +bar_y*(zF-tauF+(bar_z-bar_tau)*yF);

    nF = s_nR*nRF+s_nH*nHF;

    (1-s_g)*yF =
        s_cR*cRF+s_cH*cHF+s_i*invF+gF;

    // Shadow capital with the same default Dynare timing.
    kF = (1-delta)*kF(-1)+delta*(invF+muIF);

    qF+muIF = phi_i*((1+beta)*invF-invF(-1)-beta*invF(+1));

    qF =
        beta*(1-delta)*qF(+1)
        +beta*bar_rK*rKF(+1)
        -(RF-piF(+1));

    // Shadow production, factor demand, and price setting.
    yF   = aF+alpha*kF(-1)+(1-alpha)*nF;
    rKF  = mcF+yF-kF(-1);
    wF   = mcF+yF-nF;
    piF  = beta*piF(+1)+(eps_p-1)/(kappa_p*bar_pi^2)*mcF;

    // Shadow long-duration government bonds.
    RbF = (rho_m/bar_R)*pBF-pBF(-1);
    RF  = RbF(+1);

    // Shadow government budget constraint.
    bF =
        (1/beta)*bF(-1)
        +(b_target/beta)*(RbF-piF+yF(-1)-yF)
        -tauF+zF+gF;

    // Shadow fiscal rules retain persistence and cyclical feedback but have
    // no response to debt. Every innovation is scaled by 1-lambda.
    tauF =
        rho_tau*tauF(-1)
        +(1-rho_tau)*dy*yF(-1)
        +dDy*(yF-yF(-1))
        +(1-lambda)*e_tau;

    zF =
        rho_z*zF(-1)
        -dDn*(nF-nF(-1))
        +(1-lambda)*e_z;

    gF =
        rho_g*gF(-1)
        +(1-lambda)*e_g;

    // Passive shadow monetary policy.
    RF =
        rho_R*RF(-1)
        +(1-rho_R)*psi_piF*piF
        +(1-lambda)*mp_innov; // (SAMBA, 2023)

    // Shadow structural shocks.
    aF  = rho_a*aF(-1)+(1-lambda)*e_a; // (SAMBA, 2023)
    xiF = rho_d*xiF(-1)+(1-lambda)*e_d; // (SAMBA, 2023)
    muIF = rho_i*muIF(-1)+(1-lambda)*e_i; // (SAMBA, 2023)

    // Shadow technology-based output target.
    yFstar = aF;

    // ========================================================================
    // EMPIRICAL AUXILIARIES
    // ========================================================================

    // Aggregate consumption is a log deviation. The weights are the
    // steady-state expenditure shares of Ricardian and HtM consumption.
    c_model = (s_cR/s_c)*cR+(s_cH/s_c)*cH;

    // The geometric coupon bond has price P^B=1/(Y^B-rho_m). Therefore the
    // first-order log deviation of its yield is
    // -(bar_R-rho_m)/bar_R times the log deviation of its price. This is a
    // yield object, not the realized holding-period return Rb.
    // The one-year yield is the expectations component implied by the
    // geometric long-bond price. Its observation equation has an iid error.
    y1_model = -(bar_R-rho_m)/bar_R*pB;

    // Fiscal variables g, z, and primary_surplus are already level deviations
    // of ratios to quarterly output. Multiplication by 100 converts the ratios
    // into percentage points; no extra subtraction of y is appropriate.
    g_ratio_model  = 100*(bar_g+g);
    z_ratio_model  = 100*(bar_z+z);
    pb_ratio_model = 100*(bar_tau-bar_z-bar_g+primary_surplus);

    // ========================================================================
    // MEASUREMENT EQUATIONS: DATA_TANK.XLS, 2000Q1--2019Q4
    // ========================================================================

    // Iid measurement errors. Output and yield use SAMBA-style volatility caps;
    // debt uses the observable's own annual-GDP percentage-point units.
    me_y  = sigma_me_y*me_cap_y*e_me_y; // Adapted from SAMBA (2023).
    me_y1 = sigma_me_y1*me_cap_y1*e_me_y1; // Adapted from SAMBA (2023).
    me_b  = sigma_me_b*e_me_b; // Smets and Wouters (2024).

    // Quarter-on-quarter log growth of seasonally adjusted real GDP,
    // multiplied by 100 and not annualized. Source: IBGE SCNT Table 6612.
    y_obs = gamma_obs+100*(y-y(-1))+me_y; // (SAMBA, 2023)

    // Quarter-on-quarter log growth of seasonally adjusted real household
    // consumption, multiplied by 100 and not annualized. Source: IBGE SCNT.
    c_obs = gamma_obs+100*(c_model-c_model(-1));

    // Quarter-on-quarter log growth of seasonally adjusted real FBCF,
    // multiplied by 100 and not annualized. Source: IBGE SCNT.
    inv_obs = gamma_obs+100*(inv-inv(-1));

    // IPCA end-of-quarter log inflation, multiplied by 100 and not annualized.
    // pi is the log deviation of quarterly gross inflation from pi_target.
    pi_obs = pi_ss_pct+100*pi;

    // First difference of the effective quarterly net Selic in percentage
    // points. R is the log deviation of the quarterly gross nominal rate, so
    // d[100*(bar_R*exp(R)-1)] = 100*bar_R*dR at first order.
    dR_obs = 100*bar_R*(R-R(-1));

    // First difference of the effective quarterly one-year yield in percentage
    // points, augmented by its iid measurement error.
    dy1_obs = 100*bar_R*(y1_model-y1_model(-1))+me_y1;

    // First differences of the seasonally adjusted fiscal ratios. Sources:
    // STN/RTN Tables 4.1 and 1.2, with BCB SGS 4380 as common denominator.
    dg_obs  = g_ratio_model-g_ratio_model(-1);
    dz_obs  = z_ratio_model-z_ratio_model(-1);
    dpb_obs = pb_ratio_model-pb_ratio_model(-1);

    // Gross debt is kept in percent of annual GDP. Since b and b_target are
    // ratios to quarterly GDP, 100/4=25 appears exactly once in this mapping.
    // Source: BCB SGS 4537.
    b_obs = 25*(b_target+b)+me_b;

end;




// ============================================================================
// EPILOGUE VARIABLES: FUNDED AND UNFUNDED COMPONENTS
// ============================================================================
//
// These variables implement the accounting decomposition used in the
// Smets-Wouters funded-versus-unfunded historical decomposition.
//
// The shadow economy measures the fiscal-led component. The residual between
// the actual and shadow economies measures the funded monetary-led component.
//
// These definitions do not add equations to the structural model and do not
// enlarge the state space.
//
// ============================================================================

epilogue;
    pi_funded   = pi-piF;
    pi_unfunded = piF;

    pb_unfunded = tauF-zF-gF;
    pb_funded   = primary_surplus-pb_unfunded;
end;





// ============================================================================
// DETERMINISTIC STEADY STATE
// ============================================================================

steady_state_model;

    // Long-run policy targets.
    // pi_target is derived from the estimated quarterly-percent intercept.
    pi_target = exp(pi_ss_pct/100);
    bar_pi = pi_target;
    bar_b  = b_target;
    bar_g  = g_target;
    bar_z  = z_target;

    // Interest rates and the fiscal identity.
    bar_R   = bar_pi/beta;
    bar_tau = bar_z+bar_g+(1/beta-1)*bar_b;

    // Firms, capital, and factor prices.
    bar_mc  = (eps_p-1)/eps_p;
    bar_rK  = 1/beta-1+delta;
    bar_KY  = alpha*bar_mc/bar_rK;
    bar_IY  = delta*bar_KY;
    bar_YN  = A_ss^(1/(1-alpha))
              *bar_KY^(alpha/(1-alpha));
    bar_w   = (1-alpha)*bar_mc*bar_YN;

    // Positive TANK household allocation.
    [bar_cR,bar_cH,bar_nR,bar_nH,bar_y] =
        tank_steady_state(
            beta,h,sigma,varphi,omega,
            bar_w,bar_YN,bar_IY,
            bar_g,bar_z,bar_tau);

    // Capital and investment are held only by Ricardian households.
    bar_iR = bar_IY*bar_y/(1-omega);

    // All variables in model(linear) are deviations from this reference.
    cR = 0;
    cH = 0;
    nR = 0;
    nH = 0;
    lamR = 0;
    lamH = 0;
    n = 0;
    y = 0;
    inv = 0;
    k = 0;
    q = 0;
    muI = 0;
    rK = 0;
    w = 0;
    mc = 0;
    pi = 0;
    R = 0;
    pB = 0;
    Rb = 0;
    b = 0;
    tau = 0;
    z = 0;
    g = 0;
    a = 0;
    xi = 0;
    ystar = 0;
    r_real = 0;
    primary_surplus = 0;
    mp_innov = 0; // (SAMBA, 2023)

    c_model = 0;
    y1_model = 0;
    g_ratio_model = 100*bar_g;
    z_ratio_model = 100*bar_z;
    pb_ratio_model = 100*(bar_tau-bar_z-bar_g);
    me_y = 0;
    me_y1 = 0;
    me_b = 0;

    y_obs = gamma_obs;
    c_obs = gamma_obs;
    inv_obs = gamma_obs;
    pi_obs = pi_ss_pct;
    dR_obs = 0;
    dy1_obs = 0;
    dg_obs = 0;
    dz_obs = 0;
    dpb_obs = 0;
    b_obs = 25*b_target;

    cRF = 0;
    cHF = 0;
    nRF = 0;
    nHF = 0;
    lamRF = 0;
    lamHF = 0;
    nF = 0;
    yF = 0;
    invF = 0;
    kF = 0;
    qF = 0;
    muIF = 0;
    rKF = 0;
    wF = 0;
    mcF = 0;
    piF = 0;
    RF = 0;
    pBF = 0;
    RbF = 0;
    bF = 0;
    tauF = 0;
    zF = 0;
    gF = 0;
    aF = 0;
    xiF = 0;
    yFstar = 0;

end;

// ============================================================================
// OBSERVED VARIABLES
// ============================================================================

varobs
    y_obs
    c_obs
    inv_obs
    pi_obs
    dR_obs
    dy1_obs
    dg_obs
    dz_obs
    dpb_obs
    b_obs
;


// ============================================================================
// MODEL CHECKS AND SHOCK CALIBRATION
// ============================================================================

steady;
resid;
model_diagnostics;
check;

// ============================================================================
// OUTPUT TREE AND CALIBRATED STRUCTURAL DIAGNOSTICS
// ============================================================================
//
// All paths are derived from M_.fname so that robustness specifications never
// overwrite one another. The static-Jacobian diagnostics are evaluated at the
// deterministic steady state before any posterior parameter vector is loaded.
//
// ============================================================================

figures_root_v08 = [M_.fname '_figures'];

results_root_v08 = [M_.fname '_results'];
calibrated_irf_results_dir_v08 = fullfile(results_root_v08,'calibrated_irfs');
posterior_irf_results_dir_v08 = fullfile(results_root_v08,'posterior_irfs');
posterior_summary_results_dir_v08 = fullfile(results_root_v08,'posterior_summary');
historical_results_dir_v08 = fullfile(results_root_v08,'historical_decomposition');
diagnostics_results_dir_v08 = fullfile(results_root_v08,'diagnostics');
rank_results_dir_v08 = fullfile(results_root_v08,'rank_tests');
identification_results_dir_v08 = fullfile(results_root_v08,'identification');
logs_results_dir_v08 = fullfile(results_root_v08,'logs');

required_output_directories_v08 = {calibrated_irf_results_dir_v08;posterior_irf_results_dir_v08;posterior_summary_results_dir_v08;historical_results_dir_v08;diagnostics_results_dir_v08;rank_results_dir_v08;identification_results_dir_v08;logs_results_dir_v08};

for directory_index_v08 = 1:length(required_output_directories_v08)
    if exist(required_output_directories_v08{directory_index_v08},'dir') ~= 7
        mkdir(required_output_directories_v08{directory_index_v08});
    end
end

M_.dname = results_root_v08;

static_model_handle_v08 = str2func([M_.fname '.static']);
[static_residuals_v08,static_jacobian_v08] = static_model_handle_v08(oo_.steady_state,oo_.exo_steady_state,M_.params);
static_singular_values_v08 = svd(full(static_jacobian_v08));
static_largest_singular_value_v08 = static_singular_values_v08(1);
static_smallest_singular_value_v08 = static_singular_values_v08(end);
static_rank_tolerance_v08 = max(size(static_jacobian_v08))*eps(static_largest_singular_value_v08);
static_jacobian_rank_v08 = sum(static_singular_values_v08>static_rank_tolerance_v08);
static_condition_number_v08 = static_largest_singular_value_v08/static_smallest_singular_value_v08;
maximum_static_residual_v08 = max(abs(static_residuals_v08));
if isempty(options_.qz_criterium)
    bk_eigenvalue_threshold_v08 = 1+1e-6;
elseif isscalar(options_.qz_criterium)
    bk_eigenvalue_threshold_v08 = options_.qz_criterium;
else
    error('options_.qz_criterium must be empty or scalar.');
end

unstable_root_count_v08 = sum(abs(oo_.dr.eigval(:))>bk_eigenvalue_threshold_v08);
forward_looking_count_v08 = M_.nsfwrd;
bk_rank_count_matches_v08 = unstable_root_count_v08==forward_looking_count_v08;

if ~isfinite(maximum_static_residual_v08) || ~isfinite(static_largest_singular_value_v08) || static_smallest_singular_value_v08 <= 0
    error('The calibrated structural diagnostics contain a nonfinite residual or a singular static Jacobian.');
end

if static_jacobian_rank_v08 ~= min(size(static_jacobian_v08))
    error('The calibrated static Jacobian is rank deficient.');
end

if ~bk_rank_count_matches_v08
    error('The calibrated number of unstable roots does not match the number of forward-looking variables.');
end

structural_diagnostics_file_v08 = fullfile(diagnostics_results_dir_v08,'structural_diagnostics.mat');
save(structural_diagnostics_file_v08,'maximum_static_residual_v08','static_jacobian_rank_v08','static_rank_tolerance_v08','static_largest_singular_value_v08','static_smallest_singular_value_v08','static_condition_number_v08','unstable_root_count_v08','forward_looking_count_v08','bk_rank_count_matches_v08','static_singular_values_v08');

structural_diagnostics_csv_v08 = fullfile(diagnostics_results_dir_v08,'structural_diagnostics.csv');
structural_diagnostics_fid_v08 = fopen(structural_diagnostics_csv_v08,'w');

if structural_diagnostics_fid_v08 < 0
    error('Could not open the structural-diagnostics CSV file.');
end

fprintf(structural_diagnostics_fid_v08,'specification,original_endogenous,generated_endogenous,equations,max_static_residual,static_rank,static_rank_tolerance,smallest_singular_value,largest_singular_value,condition_number,unstable_roots,forward_looking_variables,bk_count_match\n');
fprintf(structural_diagnostics_fid_v08,'%s,%d,%d,%d,%.15g,%d,%.15g,%.15g,%.15g,%.15g,%d,%d,%d\n',M_.fname,M_.orig_endo_nbr,M_.endo_nbr,M_.eq_nbr,maximum_static_residual_v08,static_jacobian_rank_v08,static_rank_tolerance_v08,static_smallest_singular_value_v08,static_largest_singular_value_v08,static_condition_number_v08,unstable_root_count_v08,forward_looking_count_v08,bk_rank_count_matches_v08);
fclose(structural_diagnostics_fid_v08);

fprintf('\nCALIBRATED STRUCTURAL DIAGNOSTICS: %s\n',M_.fname);
fprintf('Original endogenous variables: %d; generated endogenous variables: %d; equations: %d\n',M_.orig_endo_nbr,M_.endo_nbr,M_.eq_nbr);
fprintf('Maximum absolute static residual: %.12g\n',maximum_static_residual_v08);
fprintf('Static Jacobian rank: %d; smallest singular value: %.12g; condition number: %.12g\n',static_jacobian_rank_v08,static_smallest_singular_value_v08,static_condition_number_v08);
fprintf('Unstable roots: %d; forward-looking variables: %d; BK count match: %d\n',unstable_root_count_v08,forward_looking_count_v08,bk_rank_count_matches_v08);

// Optional structural report. It is left commented out because Dynare's
// block-normalization diagnostic can issue a nonfatal matching warning for
// the long-bond return variable even when the full model is well posed.
// model_info;

shocks;

    // Nonfiscal innovations are in decimal log-deviation units.
    var e_a;
    stderr 0.039712;

    var e_d;
    stderr 0.010895;

    // Provisional MEI calibration used only before estimation.
    var e_i;
    stderr 0.010000;

    // Converted from quarterly net-rate percentage points to a log deviation
    // of the gross nominal interest rate.
    var e_mp;
    stderr 0.009232;

    // Fiscal innovations are decimal level deviations of fiscal-flow ratios.
    var e_tau;
    stderr 0.0002816;

    var e_z;
    stderr 0.0001516;

    var e_g;
    stderr 0.00008123;

    // Unit-variance innovations; their observable-unit scales are parameters
    // in the three iid measurement equations.
    var e_me_y;
    stderr 1;

    var e_me_y1;
    stderr 1;

    var e_me_b;
    stderr 1;

end;






// ============================================================================


// ============================================================================
// FIGURE AND GRID CONTROL
// ============================================================================

save_mechanism_grids = 1;
keep_mechanism_figures_open = 0;
irf_periods = 40;
period_axis = 0:(irf_periods-1);
lambda_original = get_param_by_name('lambda');
omega_original = get_param_by_name('omega');
lambda_values = [0.00 0.25 0.50 0.75 1.00];
number_lambda_values = length(lambda_values);

// Editable baseline for Exercise 3: lambda is held fixed while omega varies.
lambda_for_omega_exercise = 0.75;
omega_values = [0.00 0.05 0.10 0.15 0.20 0.25];
number_omega_values = length(omega_values);
// Safety regularization retained in case the all-HtM boundary is reintroduced;
// tank_steady_state.m excludes the exact omega=1 endpoint.
omega_boundary_epsilon = 1e-6;

grid_font_size = 20;
grid_dpi = 200;
grid_position = [40 40 2600 1800];

figures_root_directory = fullfile(pwd,[M_.fname '_figures']);
mechanism_results_directory = fullfile(pwd,[M_.fname '_results']);

// Clean only the mechanism-figure root. Existing subfolders and files
// are intentionally removed so every run starts with exactly three folders.
if exist(figures_root_directory,'dir') ~= 7
    mkdir(figures_root_directory);
else
    existing_figure_entries = dir(figures_root_directory);
    for figure_entry_index = 1:length(existing_figure_entries)
        figure_entry = existing_figure_entries(figure_entry_index);
        if strcmp(figure_entry.name,'.') || strcmp(figure_entry.name,'..')
            continue;
        end
        figure_entry_path = fullfile(figures_root_directory,figure_entry.name);
        if figure_entry.isdir
            [removed_entry_v08,removal_message_v08] = rmdir(figure_entry_path,'s');
        else
            delete(figure_entry_path);
            removed_entry_v08 = exist(figure_entry_path,'file') ~= 2;
            removal_message_v08 = 'delete failed';
        end
        if ~removed_entry_v08
            error('Could not remove old mechanism-figure entry %s: %s',figure_entry.name,removal_message_v08);
        end
    end
end

lambda_grid_directory = fullfile(figures_root_directory,'lambda_grid');
omega_grid_directory = fullfile(figures_root_directory,'omega_grid');
mechanism_grid_directory = fullfile(figures_root_directory,'mechanism_grid');

required_mechanism_directories = {lambda_grid_directory;omega_grid_directory;mechanism_grid_directory;mechanism_results_directory};
for directory_index = 1:length(required_mechanism_directories)
    if exist(required_mechanism_directories{directory_index},'dir') ~= 7
        mkdir(required_mechanism_directories{directory_index});
    end
end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

light_red = [1.00 0.87 0.87];
dark_red = [0.48 0.00 0.00];
red_colormap = [linspace(light_red(1),dark_red(1),number_lambda_values)' linspace(light_red(2),dark_red(2),number_lambda_values)' linspace(light_red(3),dark_red(3),number_lambda_values)'];
light_blue = [0.86 0.92 1.00];
dark_blue = [0.00 0.18 0.55];
blue_colormap = [linspace(light_blue(1),dark_blue(1),number_omega_values)' linspace(light_blue(2),dark_blue(2),number_omega_values)' linspace(light_blue(3),dark_blue(3),number_omega_values)'];

// ============================================================================
// LAMBDA-SENSITIVITY GRIDS FOR THE SEVEN ECONOMIC INNOVATIONS
// ============================================================================

shocks_to_plot = {'e_a','e_d','e_i','e_mp','e_tau','e_z','e_g'};
shock_display_names = {'Technology shock','Preference shock','Investment-efficiency shock','Monetary-policy shock','Government-revenue shock','Government-transfer shock','Government-purchases shock'};
shock_file_names = {'technology','preference','investment_efficiency','monetary_policy','government_revenue','government_transfer','government_purchases'};

var_list_to_plot = {'pi','y','R','r_real','primary_surplus','b'};
title_list = {'Inflation','Output','Nominal rate','Real rate','Primary surplus','Debt ratio'};
ylabel_list = {'Annualized p.p.','Percent from steady state','Annualized p.p.','Annualized p.p.','p.p. of GDP','p.p. of annual GDP'};
scale_factors = [400 100 400 400 100 25];

number_shocks = length(shocks_to_plot);
number_variables = length(var_list_to_plot);
if number_shocks ~= 7 || number_variables ~= 6
    error('The v08 mechanism grid requires seven economic shocks and six plotted variables.');
end

lambda_irf_results = struct();
for shock_index = 1:number_shocks
    for variable_index = 1:number_variables
        lambda_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index}) = zeros(irf_periods,number_lambda_values);
    end
end

for lambda_index = 1:number_lambda_values
    current_lambda = lambda_values(lambda_index);
    fprintf('Computing mechanism grids for lambda = %.2f\n',current_lambda);
    set_param_value('lambda',current_lambda);
    stoch_simul(order=1,irf=40,nograph,noprint,irf_shocks=(e_a,e_d,e_i,e_mp,e_tau,e_z,e_g)) pi y R r_real primary_surplus b;

    if info(1) ~= 0
        set_param_value('lambda',lambda_original);
        error('IRF computation failed for lambda = %.2f with Dynare code %d.',current_lambda,info(1));
    end

    for shock_index = 1:number_shocks
        for variable_index = 1:number_variables
            irf_field_name = [var_list_to_plot{variable_index} '_' shocks_to_plot{shock_index}];
            if isfield(oo_.irfs,irf_field_name)
                current_irf = oo_.irfs.(irf_field_name);
                available_periods = min(irf_periods,length(current_irf));
                lambda_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index})(1:available_periods,lambda_index) = current_irf(1:available_periods);
            end
        end
    end
end

set_param_value('lambda',lambda_original);
stoch_simul(order=1,irf=0,nograph,noprint);

generated_png_files = cell(0,1);
for shock_index = 1:number_shocks
    mechanism_figure = figure('Visible','off','Name',shock_display_names{shock_index}, ...
        'NumberTitle','off','Color','w','Units','pixels','Position',grid_position);
    axes_array = gobjects(1,number_variables);

    for variable_index = 1:number_variables
        axes_array(variable_index) = subplot(2,3,variable_index);
        current_axes = axes_array(variable_index);
        hold(current_axes,'on');
        yline(current_axes,0,'-','Color',[0.30 0.30 0.30], ...
            'LineWidth',0.60,'HandleVisibility','off');
        scaled_irfs = scale_factors(variable_index)*lambda_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index});

        for lambda_index = 1:number_lambda_values
            line_width = 0.95;
            if lambda_index == number_lambda_values
                line_width = 1.80;
            end
            plot(current_axes,period_axis,scaled_irfs(:,lambda_index), ...
                'Color',red_colormap(lambda_index,:),'LineWidth',line_width, ...
                'HandleVisibility','off');
        end

        title(current_axes,title_list{variable_index},'FontSize',grid_font_size);
        xlabel(current_axes,'Quarters','FontSize',grid_font_size);
        ylabel(current_axes,ylabel_list{variable_index},'FontSize',grid_font_size);
        xlim(current_axes,[0 irf_periods-1]);
        grid(current_axes,'on');
        box(current_axes,'off');
        set(current_axes,'FontSize',grid_font_size,'Layer','top');
        hold(current_axes,'off');
    end

    legend_axes = axes_array(4);
    hold(legend_axes,'on');
    dummy_lines = gobjects(number_lambda_values,1);
    legend_labels = cell(number_lambda_values,1);
    for legend_index = 1:number_lambda_values
        lambda_index = number_lambda_values-legend_index+1;
        dummy_lines(legend_index) = plot(legend_axes,NaN,NaN,'-', ...
            'Color',red_colormap(lambda_index,:),'LineWidth',1.4);
        legend_labels{legend_index} = sprintf('$\\lambda=%.2f$',lambda_values(lambda_index));
    end
    legend(legend_axes,dummy_lines,legend_labels,'Location','best', ...
        'FontSize',grid_font_size,'Box','off');
    hold(legend_axes,'off');

    drawnow;
    png_file_name = fullfile(lambda_grid_directory, ...
        ['irf_grid_' shock_file_names{shock_index} '.png']);
    print(mechanism_figure,png_file_name,'-dpng',['-r' num2str(grid_dpi)]);
    if exist(png_file_name,'file') ~= 2 || dir(png_file_name).bytes == 0
        error('The required lambda grid PNG was not created: %s',png_file_name);
    end
    imfinfo(png_file_name);
    generated_png_files{end+1,1} = png_file_name;

    if keep_mechanism_figures_open == 0
        close(mechanism_figure);
    else
        set(mechanism_figure,'Visible','on');
    end
end

// ============================================================================
// TECHNOLOGY-SHOCK MECHANISM GRID (EXERCISE 2)
// ============================================================================

diagnostic_variables = {'cR','cH','c_model','k','mc','w'};
diagnostic_titles = {'Ricardian consumption','HtM consumption','Aggregate consumption','Capital','Marginal cost','Real wage'};
diagnostic_scales = [100 100 100 100 100 100];
diagnostic_ylabels = {'Percent','Percent','Percent','Percent','Percent','Percent'};

e_a_results = zeros(irf_periods,6,number_lambda_values);

for lambda_index = 1:number_lambda_values
    set_param_value('lambda',lambda_values(lambda_index));
    stoch_simul(order=1,irf=40,nograph,noprint,irf_shocks=(e_a)) cR cH c_model k mc w;

    if info(1) ~= 0
        set_param_value('lambda',lambda_original);
        error('Technology-shock mechanism grid failed for lambda = %.2f.',lambda_values(lambda_index));
    end

    for variable_index = 1:6
        diagnostic_field = [diagnostic_variables{variable_index} '_e_a'];
        if isfield(oo_.irfs,diagnostic_field)
            current_irf = oo_.irfs.(diagnostic_field);
            e_a_results(1:min(irf_periods,length(current_irf)),variable_index,lambda_index) = current_irf(1:min(irf_periods,length(current_irf)));
        end
    end
end

set_param_value('lambda',lambda_original);
diagnostic_figure = figure('Visible','off','Name','Technology-shock mechanism', ...
    'NumberTitle','off','Color','w','Units','pixels','Position',grid_position);
diagnostic_axes = gobjects(1,6);
for variable_index = 1:6
    diagnostic_axes(variable_index) = subplot(2,3,variable_index);
    current_axes = diagnostic_axes(variable_index);
    hold(current_axes,'on');
    yline(current_axes,0,'-','Color',[0.30 0.30 0.30], ...
        'LineWidth',0.60,'HandleVisibility','off');
    for lambda_index = 1:number_lambda_values
        plot(current_axes,period_axis,diagnostic_scales(variable_index)* ...
            e_a_results(:,variable_index,lambda_index), ...
            'Color',red_colormap(lambda_index,:),'LineWidth',1.40, ...
            'HandleVisibility','off');
    end
    title(current_axes,diagnostic_titles{variable_index},'FontSize',grid_font_size);
    xlabel(current_axes,'Quarters','FontSize',grid_font_size);
    ylabel(current_axes,diagnostic_ylabels{variable_index},'FontSize',grid_font_size);
    xlim(current_axes,[0 irf_periods-1]);
    grid(current_axes,'on');
    box(current_axes,'off');
    set(current_axes,'FontSize',grid_font_size,'Layer','top');
    hold(current_axes,'off');
end

legend_axes = diagnostic_axes(1);
hold(legend_axes,'on');
dummy_lines = gobjects(number_lambda_values,1);
legend_labels = cell(number_lambda_values,1);
for legend_index = 1:number_lambda_values
    lambda_index = number_lambda_values-legend_index+1;
    dummy_lines(legend_index) = plot(legend_axes,NaN,NaN,'-', ...
        'Color',red_colormap(lambda_index,:),'LineWidth',1.4);
    legend_labels{legend_index} = sprintf('$\\lambda=%.2f$',lambda_values(lambda_index));
end
legend(legend_axes,dummy_lines,legend_labels,'Location','best', ...
    'FontSize',grid_font_size,'Box','off');
hold(legend_axes,'off');

drawnow;
png_file_name = fullfile(mechanism_grid_directory,'irf_grid_technology_mechanism.png');
print(diagnostic_figure,png_file_name,'-dpng',['-r' num2str(grid_dpi)]);
if exist(png_file_name,'file') ~= 2 || dir(png_file_name).bytes == 0
    error('The technology-shock mechanism PNG was not created: %s',png_file_name);
end
imfinfo(png_file_name);
generated_png_files{end+1,1} = png_file_name;

if keep_mechanism_figures_open == 0
    close(diagnostic_figure);
else
    set(diagnostic_figure,'Visible','on');
end

// ============================================================================
// OMEGA-SENSITIVITY GRIDS (EXERCISE 3)
// ============================================================================

omega_irf_results = struct();
for shock_index = 1:number_shocks
    for variable_index = 1:number_variables
        omega_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index}) = zeros(irf_periods,number_omega_values);
    end
end

for omega_index = 1:number_omega_values
    current_omega = omega_values(omega_index);
    omega_for_solver = current_omega;
    if current_omega >= 1
        omega_for_solver = 1-omega_boundary_epsilon;
    end

    set_param_value('lambda',lambda_for_omega_exercise);
    set_param_value('omega',omega_for_solver);
    steady;
    if info(1) ~= 0
        set_param_value('lambda',lambda_original);
        set_param_value('omega',omega_original);
        error('Steady state failed for omega = %.2f (solver value %.8f).',current_omega,omega_for_solver);
    end

    fprintf('Computing mechanism grids for omega = %.2f at lambda = %.2f\n',current_omega,lambda_for_omega_exercise);
    stoch_simul(order=1,irf=40,nograph,noprint,irf_shocks=(e_a,e_d,e_i,e_mp,e_tau,e_z,e_g)) pi y R r_real primary_surplus b;
    if info(1) ~= 0
        set_param_value('lambda',lambda_original);
        set_param_value('omega',omega_original);
        error('IRF computation failed for omega = %.2f at lambda = %.2f.',current_omega,lambda_for_omega_exercise);
    end

    for shock_index = 1:number_shocks
        for variable_index = 1:number_variables
            irf_field_name = [var_list_to_plot{variable_index} '_' shocks_to_plot{shock_index}];
            if isfield(oo_.irfs,irf_field_name)
                current_irf = oo_.irfs.(irf_field_name);
                available_periods = min(irf_periods,length(current_irf));
                omega_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index})(1:available_periods,omega_index) = current_irf(1:available_periods);
            end
        end
    end
end

set_param_value('lambda',lambda_original);
set_param_value('omega',omega_original);
steady;

for shock_index = 1:number_shocks
    omega_figure = figure('Visible','off','Name', ...
        ['Omega sensitivity - ' shock_display_names{shock_index}], ...
        'NumberTitle','off','Color','w','Units','pixels','Position',grid_position);
    omega_axes = gobjects(1,number_variables);

    for variable_index = 1:number_variables
        omega_axes(variable_index) = subplot(2,3,variable_index);
        current_axes = omega_axes(variable_index);
        hold(current_axes,'on');
        yline(current_axes,0,'-','Color',[0.30 0.30 0.30], ...
            'LineWidth',0.60,'HandleVisibility','off');
        scaled_irfs = scale_factors(variable_index)*omega_irf_results.(shocks_to_plot{shock_index}).(var_list_to_plot{variable_index});

        for omega_index = 1:number_omega_values
            line_width = 1.40;
            if omega_index == number_omega_values
                line_width = 1.90;
            end
            plot(current_axes,period_axis,scaled_irfs(:,omega_index), ...
                'Color',blue_colormap(omega_index,:),'LineWidth',line_width, ...
                'HandleVisibility','off');
        end

        title(current_axes,title_list{variable_index},'FontSize',grid_font_size);
        xlabel(current_axes,'Quarters','FontSize',grid_font_size);
        ylabel(current_axes,ylabel_list{variable_index},'FontSize',grid_font_size);
        xlim(current_axes,[0 irf_periods-1]);
        grid(current_axes,'on');
        box(current_axes,'off');
        set(current_axes,'FontSize',grid_font_size,'Layer','top');
        hold(current_axes,'off');
    end

    legend_axes = omega_axes(1);
    hold(legend_axes,'on');
    dummy_lines = gobjects(number_omega_values,1);
    legend_labels = cell(number_omega_values,1);
    for legend_index = 1:number_omega_values
        omega_index = number_omega_values-legend_index+1;
        dummy_lines(legend_index) = plot(legend_axes,NaN,NaN,'-', ...
            'Color',blue_colormap(omega_index,:),'LineWidth',1.4);
        legend_labels{legend_index} = sprintf('$\\omega=%.2f$',omega_values(omega_index));
    end
    legend(legend_axes,dummy_lines,legend_labels,'Location','best', ...
        'FontSize',grid_font_size,'Box','off');
    hold(legend_axes,'off');

    drawnow;
    png_file_name = fullfile(omega_grid_directory, ...
        ['irf_grid_omega_' shock_file_names{shock_index} '.png']);
    print(omega_figure,png_file_name,'-dpng',['-r' num2str(grid_dpi)]);
    if exist(png_file_name,'file') ~= 2 || dir(png_file_name).bytes == 0
        error('The omega-sensitivity PNG was not created: %s',png_file_name);
    end
    imfinfo(png_file_name);
    generated_png_files{end+1,1} = png_file_name;

    if keep_mechanism_figures_open == 0
        close(omega_figure);
    else
        set(omega_figure,'Visible','on');
    end
end

// ============================================================================
// OUTPUT VALIDATION
// ============================================================================

expected_grid_count = 2*number_shocks+1;
if length(generated_png_files) ~= expected_grid_count
    error('The v08 mechanism file expected exactly %d composite PNG grids.',expected_grid_count);
end

for output_index = 1:length(generated_png_files)
    output_info = dir(generated_png_files{output_index});
    if isempty(output_info) || output_info.bytes == 0
        error('A mechanism PNG is missing or empty: %s',generated_png_files{output_index});
    end
end

save(fullfile(mechanism_results_directory,'mechanism_grid_results.mat'),'lambda_values','lambda_irf_results','e_a_results','omega_values','lambda_for_omega_exercise','omega_irf_results','generated_png_files');
fprintf('Verified %d nonempty composite PNG grids under %s\n',length(generated_png_files),figures_root_directory);
fprintf('Mechanism model counts: %d original endogenous variables, %d generated endogenous variables, %d equations, %d innovations, %d observables.\n',M_.orig_endo_nbr,M_.endo_nbr,M_.eq_nbr,M_.exo_nbr,length(options_.varobs));
fprintf('No estimation or MCMC was executed.\n');
try
    if ispc
        winopen(figures_root_directory);
    elseif ismac
        system(['open "' figures_root_directory '"']);
    elseif isunix
        system(['xdg-open "' figures_root_directory '"']);
    end
catch
    fprintf('Could not open the figures folder automatically: %s\n',figures_root_directory);
end
