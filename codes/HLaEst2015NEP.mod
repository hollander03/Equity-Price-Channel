// FROM WORKING VERSION: HLaEst2015BEP.mod // fix delta_B=0.4
// PAPER 1 2015 10 AUGUST: RE-ESTIMATION with updated sample
// SAVER HOUSEHOLD FLOW OF FUNDS CONSTRAINT; entrep flow of funds constraint
// Note: SAMPLE PERIOD 1982:01 - 2015:01 USED.

var c_b i_d pi i_h lambda l_h w nu_h h_s h_b
    d c_s 
    x
    y k_e h l_e q_k i_e nu_e lambda_e v 
    k_B pi_B i_l l k_BL 
    c  
    xi_z xi_d xi_i varepsilon_e varepsilon_h xi_p 
    pi_obs y_obs l_h_obs l_e_obs d_obs i_d_obs i_h_obs i_e_obs;  // div_e

varexo epsilon_p epsilon_z epsilon_i epsilon_d epsilon_h epsilon_e epsilon_nu_h epsilon_nu_e ; // epsilon_psi

parameters beta_b phi R_d R_h Nu_h phi_w eta gamma gamma_b
           beta_s
           theta_R varepsilon_p gamma_p
           alpha Nu_e phi_k R_e beta_e delta_e
           kappa_v
           phi_s 
           delta_B kappa_k tau kappa_e varepsilon_ess beta_B kappa_h varepsilon_hss L_hL L_eL
           kappa_i kappa_pi kappa_y  
           CY K_BY LY   
           rho_z rho_d rho_i rho_e rho_h rho_nuh rho_nue rho_p 
           i_h_ss i_e_ss i_d_ss pi_ss; //  //  K_eY  PsiY //  theta_w varepsilon_w gamma_w

//Calibrated Parameters

//households

beta_b = 0.96;     // 0.95 0.96 discount factor for borrower must be less than (1/R_h)
phi = 0.5;         // habit formation
R_d = 1.01;     // 1.01187 gross real return to deposit
R_h = 1.033;     // gross return to HH loans #Lambda = 1/R_h - beta_b must be bigger than zero in steady state. Therefore R_h corresponds to beta_b
Nu_h = 0.75;        // loan-to-value for HH
phi_w = 1;       // weight on wage income in borrowing constraint
eta = 1;           // inverse labour elasticity
gamma_b = 2;      // borrower RRA

beta_s = 0.99;    // 0.988 related to steady state equity return discount factor for savers = 1/R_d = opportunity cost for holding funds for savers
gamma = 2;        // constant relative risk aversion coefficient

//entrepreneurs

beta_e = 0.95;     // 0.95 0.96 entrepreneurs discount factor must be less than 1/R_e
R_e = 1.039;      // 1.0247 gross return to Entrepreneurs loans see R_h above
Nu_e = 0.65;       // loan-to-value ratio for entrepreneurs
phi_k = 1;       // weight on the value of capital in borrowing constraint
alpha = 0.33;      // 0.25 share of capital in firm production
delta_e = 0.025;    // depreciation rate of capital
varepsilon_p = 11;   // 7.667 11 // 21 (Iacoviello 2005) or 6 (Christiano et al 2010 for intermediate goods) price elasticity of demand(substitution) across differentiated retail goods
kappa_v = 2;         // physical capital adjustment cost

//Retailers and unions

theta_R = 0.65;     // probability of retail prices remaining unchanged in NKPC, theta_R = 0 for fully flexible prices
gamma_p = 0.25;     // degree of price indexation

//theta_w = 0.75;    // probability of wages remaining unchanged, theta_w = 0 for fully flexible wages
//gamma_w = 0.5;     // degree of wage indexation, gamma_w = 0 no indexation and gamma_w = 1 full indexation
//varepsilon_w = 5;  // 5 wage elasticity of substitution across different types of HH
phi_s = 0.53;       // weight on (or share of) saver households in unions/economy

//banking sector

beta_B = 0.99;            // 0.988 banking discount factor is equal to the patient HHs' discount factor in Gerali et al. 2010
delta_B = 0.4;         // Gerali 0.1049 USmodel 0.1044 Cost for managing the banks capital position
kappa_k = 4;              // parameter governing adjustment costs in banking
kappa_e = 8;            // parameter governing firm loan rate adjustment costs
kappa_h = 8;            // parameter governing HH loan rate adjustment costs
tau = 0.11;               // target capital-to-loans ratio or capital requirement (leverage ratio) this is inconsistent with the data setup!! should be 0.5-0.7
varepsilon_ess = 1.341;   // USdata 3.843 steady state markup, varepsilon_ess/(varepsilon_ess - 1) is the markup on the loan rate to entrepreneurs
varepsilon_hss = 1.427;   // USdata 3.295 steady state markup, varepsilon_hss/(varepsilon_hss - 1) is the markup on the loan rate to HHs
L_hL = 0.45;              // 0.475 HH loans to total loans ratio
L_eL = 0.55;              // 0.525 Entrepreneurial loans to total loans ratio
//DL = 0.89;                 // 0.3547 deposit loan ratio, originally 1 - tau but too small in relation to data, but 0.89 in new bank data!

//monetary policy and aggregates

kappa_i = 0.65;          // Taylor rule coefficient on i
kappa_pi = 1.5;            // Taylor rule coefficient on pi
kappa_y = 0.25;          // Taylor rule coefficient on y

CY = 0.653;             // 0.6734 0.7 parameter ratios in aggregate resource eqn
K_BY = 0.165;            // 0.1804 0.165  // 0.115 gives 0.07 K_B/L
LY = 1.5;                // 0.9 1.5 
//K_eY = 10.7;             // 8.8 // 10.7 Christiano et al. 2010 USdata

//shock parameters

rho_z = 0.75;             // AR(1) parameter for productivity shock
rho_d = 0.5;             // AR(1) parameter for deposit shock
rho_i = 0.5;             // AR(1) parameter for MP shock
rho_e = 0.5;             // AR(1) parameter for entrep shock
rho_h = 0.5;             // AR(1) parameter for HH interest elasticity  shock
rho_nuh = 0.75;         //AR(1) parameter for LTV shock
rho_nue = 0.75;         //AR(1) parameter for LTV shock
rho_p = 0.5;            //AR(1) parameter for Cost-push shock
//rho_w = 0.5;            //AR(1) parameter for wage markup shock

pi_ss = 0.005887;           // inflation steady-state
//i_ss = 0.05;
i_h_ss = 0.078998;         // steady state interest rate level
i_e_ss = 0.08487;
i_d_ss = 0.04558;         // fed funds i_ss replaces i_d_ss 0.0596 in this setup

model(linear);

#Gamma_w = ((1/R_h) - beta_b)*Nu_h*phi_w;                                    // 
#gamma_ds = 0.789;     // 0.7762 calibrated

#Upsilon_k = ((1/R_e) - beta_e)*Nu_e*phi_k;            // parameter in FOC K_e

#X = varepsilon_p/(varepsilon_p - 1);                             //steady state markup
#beta_R = beta_s;                                            // discount factor for retailers

#HY = (1-alpha)/X;
#K_eY = (alpha/X)*(1/(1 - beta_e*(1 - delta_e) - Upsilon_k)); //  model steady state approx. 5

#DL = 1 - tau;                                             // bank deposit to loan ratio 0.89
#Pi_BL = (R_h - 1)*L_hL + (R_e - 1)*L_eL - (R_d - 1)*DL;   // bank profit to loans ratio approx. 0.0272

//EQUATIONS

//BORROWERS
     
l_h = (1/R_h)*(w(+1) + h_b) - i_h + (1/R_h)*pi(+1) + (1/R_h)*nu_h;                    

((1/R_h) - beta_b)*lambda = beta_b*((gamma_b/(1 - phi))*(c_b(+1) - phi*c_b) + pi(+1)) - (1/R_h)*((gamma_b/(1 - phi))*(c_b - phi*c_b(-1)) + i_h);                                                   

w = (gamma_b/(1 - phi))*(c_b - phi*c_b(-1)) + (1 + Gamma_w)*eta*h_b - Gamma_w*(lambda + nu_h + w(+1) + pi(+1));      // FOC H with w see p. 20 dynare manual

//SAVERS

d = (gamma_ds/(1 - phi))*(gamma/(1 - phi))*(c_s - phi*c_s(-1)) - (beta_s*R_d*(gamma_ds/(1 - phi)))*((gamma/(1 - phi))*(c_s(+1) - phi*c_s) - i_d + pi(+1)); // FOC for D_s

w = eta*h_s + (gamma/(1 - phi))*(c_s - phi*c_s(-1)); // 

   //Saver household flow of funds constraint
c_s + gamma_ds*(d) = ((1-alpha)/X)*(1/(CY))*(w + h_s) + R_d*gamma_ds*(i_d(-1) + d(-1) - pi);

//RETAILERS

pi = (beta_R/(1 + beta_R*gamma_p))*pi(+1) + (gamma_p/(1 + beta_R*gamma_p))*pi(-1) - (1/(1 + beta_R*gamma_p))*(((1 - theta_R)*(1 - theta_R*beta_R))/theta_R)*x + xi_p; //with indexation

// ENTREPRENEURS

h = y - x - w;                                             // FOC H
((1/R_e) - beta_e)*lambda_e = beta_e*(pi(+1)) - (1/R_e)*(i_e); // FOC L_e

(v - k_e(-1)) = (beta_e)*(v(+1) - k_e) + (1/kappa_v)*(1 - beta_e*(1 - delta_e) - Upsilon_k)*(y(+1) - x(+1) - k_e) + (1/kappa_v)*(Upsilon_k)*(lambda_e + nu_e + q_k(+1)); // FOC K_e determines the investment schedule

q_k = kappa_v*(v - k_e(-1));   // Shadow price of capital

l_e = (1/R_e)*(q_k(+1) + k_e) + (1/R_e)*pi(+1) + (1/R_e)*nu_e; // Borrowing constraint

y = alpha*k_e(-1) + (1 - alpha)*h + xi_z;        // production function
k_e = (1 - delta_e)*k_e(-1) + delta_e*v;       // Capital accumulation equation

//Unions and Wage setting equations
   // indexation
//w = (theta_w/(1 + (theta_w^2)*beta))*w(-1) + (theta_w*beta/(1 + (theta_w^2)*beta))*w(+1) + (theta_w*beta/(1 + (theta_w^2)*beta))*pi(+1) - (theta_w/(1 + (theta_w^2)*beta))*pi
//     - (((theta_w^2)*beta*gamma_w)/(1 + (theta_w^2)*beta))*pi + ((theta_w*gamma_w)/(1 + (theta_w^2)*beta))*pi(-1) 
//     + (((1 - theta_w)*(1 - theta_w*beta))/((1 + (theta_w^2)*beta)*(1 + varepsilon_w*eta)))*(varepsilon_w*eta*w + chi_s*((gamma/(1 - phi))*(c_s - phi*c_s(-1)) + eta*h) + chi_b*((gamma_b/(1 - phi))*(c_b - phi*c_b(-1)) + eta*h));  //    + xi_w
   // no indexation
//w = (theta_w/(1 + (theta_w^2)*beta))*w(-1) + (theta_w*beta/(1 + (theta_w^2)*beta))*w(+1) + (theta_w*beta/(1 + (theta_w^2)*beta))*pi(+1) - (theta_w/(1 + (theta_w^2)*beta))*pi 
//     + (((1 - theta_w)*(1 - theta_w*beta))/((1 + (theta_w^2)*beta)*(1 + varepsilon_w*eta)))*(varepsilon_w*eta*w + chi_s*((gamma/(1 - phi))*(c_s - phi*c_s(-1)) + eta*h) + (chi_b*((gamma_b/(1 - phi))*(c_b - phi*c_b(-1)) + eta*h)))
//     + epsilon_w;

//Banking Sector
 //Wholesale branch

k_B = (1 - delta_B)*k_B(-1) + delta_B*pi_B(-1) - pi; //

i_l = i_d - ((1/(R_d - 1))*kappa_k*((tau)^3))*(k_B - l);  //  - xi_t  spread between the wholesale loan and wholesale deposit rate

 //Retail: loan branch
i_e = (kappa_e/(varepsilon_ess - 1 + kappa_e + kappa_e*beta_B))*i_e(-1) + (beta_B*kappa_e/(varepsilon_ess - 1 + kappa_e + kappa_e*beta_B))*i_e(+1)
      + ((varepsilon_ess - 1)/(varepsilon_ess - 1 + kappa_e + kappa_e*beta_B))*i_l 
      - (varepsilon_e/(varepsilon_ess - 1 + kappa_e + kappa_e*beta_B));              // loan rate setting charged to entrepreneurs; the simplified case where varepsilon is a constant, drop last term. 
// not log-linearized// i_e = varepsilon_e/(varepsilon_e-1)*i_l;                       //perfectly flexible interest rates gives the markup over the wholesale loan rate, see derivation for i_e - i spread

i_h = (kappa_h/(varepsilon_hss - 1 + kappa_h + kappa_h*beta_B))*i_h(-1) + (beta_B*kappa_h/(varepsilon_hss - 1 + kappa_h + kappa_h*beta_B))*i_h(+1)
      + ((varepsilon_hss - 1)/(varepsilon_hss - 1 + kappa_h + kappa_h*beta_B))*i_l 
      - (varepsilon_h/(varepsilon_hss - 1 + kappa_h + kappa_h*beta_B));              // loan rate setting charged to HHs; the simplified case where varepsilon is a constant, drop last term.
// not log-linearized// i_h = varepsilon_h/(varepsilon_h-1)*i_l;                       // perfectly flexible interest rates gives the markup over the wholesale loan rate, see derivation for i_h - i spread

 //Retail: deposit branch
//i_d = (kappa_d/(1 - varepsilon_dss + kappa_d + kappa_d*beta_B))*i_d(-1) + (beta_B*kappa_d/(1 - varepsilon_dss + kappa_d + kappa_d*beta_B))*i_d(+1)
//      + ((1 - varepsilon_dss)/(1 - varepsilon_dss + kappa_d + kappa_d*beta_B))*i 
//      + (varepsilon_d/(1 - varepsilon_dss + kappa_d + kappa_d*beta_B));              // deposit rate setting to HHs; the simplified case where varepsilon is a constant, drop last term.
// Note: varepsilon_dss is negative, Gerali et al. 2010 working paper they write denominator with a positive sign, implying they convert varepsilon_dss to positive, perhaps for AR(1) shock.
// not log-linearized// i_d = varepsilon_d/(varepsilon_d-1)*i;                       // note for markdown varepsilon_h is negative, may need to convert to positive by subst. minus into above equation.

//Retained earnings net dividend payments
Pi_BL*pi_B = (R_h - 1)*L_hL*(i_h + l_h) + (R_e - 1)*L_eL*(i_e + l_e) - (R_d - 1)*DL*(i_d + d);  //

//MONETARY POLICY

i_d = kappa_i*i_d(-1) + kappa_pi*(1 - kappa_i)*pi + kappa_y*(1 - kappa_i)*(y - y(-1)) + xi_i; // conventional nominal interest rate (Taylor-type) rule

//MARKET CLEARING - aggregate resource constraints

y = phi_s*CY*c_s + (1-phi_s)*CY*c_b + delta_e*K_eY*v + delta_B*K_BY*k_B(-1) - (delta_B*K_BY)*pi;

c = phi_s*c_s + (1 - phi_s)*c_b;                                      // aggregate consumption
//d = phi_s*d_s + (1 - phi_s)*d_b;                                    // aggregate deposits
h = phi_s*h_s + (1 - phi_s)*h_b;
//h: aggregate hours worked

l = L_hL*l_h + L_eL*l_e;                                              // aggregate loans

l = tau*k_B + DL*(d  - xi_d);  //  - xi_d bank balance sheet identity for when deposits not in hh utility

//Shocks
xi_z = rho_z*xi_z(-1) + epsilon_z;                      // productivity shock
xi_d = rho_d*xi_d(-1) + epsilon_d;                      // deposit shock
xi_i = rho_i*xi_i(-1) + epsilon_i;                      // MP shock
varepsilon_e = rho_e*varepsilon_e(-1) + epsilon_e;      // interest elasticity shocks
varepsilon_h = rho_h*varepsilon_h(-1) + epsilon_h;      //
nu_h = rho_nuh*nu_h(-1) + epsilon_nu_h;
nu_e = rho_nue*nu_e(-1) + epsilon_nu_e;                 // LTV shocks
xi_p = rho_p*xi_p(-1) + epsilon_p;              // Cost-push shock
//xi_w = rho_w*xi_w(-1) + epsilon_w;              // wage shock

//measurement equation for all observed variables!!

k_BL = k_B - l;   // bank capital-asset ratio

pi_obs = pi + pi_ss;
y_obs = y - y(-1) + 0.001685;
l_h_obs = l_h - l_h(-1) + 0.003758;
l_e_obs = l_e - l_e(-1) + 0.0024686;
d_obs = d - d(-1) + 0.001761;
i_h_obs = i_h + i_h_ss;
i_e_obs = i_e + i_e_ss;
i_d_obs = i_d + i_d_ss;

end;

steady;
check;

estimated_params;

gamma, inv_gamma_pdf, 3, 0.25;                 // constant relative risk aversion coefficient
gamma_b, inv_gamma_pdf, 2.5, 0.25;
phi, beta_pdf, 0.5, 0.05;                // habit formation or 0.7
theta_R, beta_pdf, 0.7, 0.05;            // price stickiness or 0.75
gamma_p, beta_pdf, 0.5, 0.05;            // degree of price indexation
Nu_h, beta_pdf, 0.75, 0.05;               // LTV ratio

//parameters to be estimated only

   //entrepreneurs
Nu_e, beta_pdf, 0.55, 0.05;               // LTV ratio

   //adjustment costs
kappa_k, gamma_pdf, 4, 2;               // leverage deviation cost
kappa_e, gamma_pdf, 4, 2;             // Entrepreneur loan rate adjustment costs
kappa_h, gamma_pdf, 4, 2;             // HH loan rate adjustment costs

//coefficients in Monetary Policy function

kappa_i, beta_pdf, 0.5, 0.05;         // Taylor rule coefficient on i
kappa_pi, gamma_pdf, 2, 0.05;            // Taylor rule coefficient on pi
kappa_y, beta_pdf, 0.25, 0.05;         // Taylor rule coefficient on y

//coefficients of AR(1) shocks

rho_z, beta_pdf, 0.9, 0.02;              // AR(1) parameter for productivity shock
rho_d, beta_pdf, 0.9, 0.02;              // AR(1) parameter for deposit shock
rho_i, beta_pdf, 0.5, 0.05;              // AR(1) parameter for MP shock
rho_e, beta_pdf, 0.5, 0.05;              // AR(1) parameter for entrep shock
rho_h, beta_pdf, 0.5, 0.05;              // AR(1) parameter for HH interest elasticity  shock
rho_nuh, beta_pdf, 0.9, 0.02;
rho_nue, beta_pdf, 0.9, 0.02;            // AR(1) parameter for LTV shock
rho_p, beta_pdf, 0.5, 0.05;              // AR(1) parameter for Cost-push shock
//rho_w, beta_pdf, 0.5, 0.05;            // AR(1) parameter for wage markup shock

//innovations of the 11 shocks - 10 AR(1) shocks and 1 i.i.d shocks

stderr epsilon_z, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_i, inv_gamma_pdf, 0.01, inf;
stderr epsilon_d, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_e, inv_gamma_pdf, 0.015, inf; 
stderr epsilon_h, inv_gamma_pdf, 0.015, inf;   
stderr epsilon_nu_h, inv_gamma_pdf, 0.01, inf;
stderr epsilon_nu_e, inv_gamma_pdf, 0.01, inf;
stderr epsilon_p, inv_gamma_pdf, 0.005, inf;
//stderr epsilon_w, inv_gamma_pdf, 0.002, inf;

end;

varobs pi_obs y_obs l_h_obs l_e_obs d_obs i_d_obs i_h_obs i_e_obs;  // 8 shock NEP   

estimation(datafile=HL2015EstimationNEPshock, mh_nblocks=3, mh_replic=10000, mh_drop=0.5, mh_jscale=0.3, 
           mode_compute=0, mode_file=HLaEst2013NEP4bNEW_mode, mode_check, bayesian_irf, irf=20, moments_varendo, conditional_variance_decomposition=[1:20],
           graph_format=fig) y pi i_d i_e i_h d l_e l_h k_BL l i_l pi_B v q_k k_e c_s c_b w h; // w_s w_b h_s h_b 
                             
//shock_decomposition pi_obs y_obs i_d_obs;
