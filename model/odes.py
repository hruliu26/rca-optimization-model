import numpy as np
from scipy.integrate import solve_ivp

# -------------------------
# conversion
# -------------------------
def convert_enzyme_UuL_to_uM(enzyme_u_per_uL):
    """
    Convert enzyme activity units (U/uL) to an approximate enzyme concentration in µM.
    enzyme_u_per_uL: numeric (U/uL)

    For WT phi29 DNA polymerase, molecular weight is 67 kDa with a specific activity of 100 U/ug.

    Returns approximate enzyme concentration in µM.
    """
    specific_activity_U_per_ug = 100.0  # U/ug for phi29
    specific_activity_U_per_g = specific_activity_U_per_ug * 1e6  # U/g

    mw_kDa = 67.0  # kDa
    mw_g_per_mol = mw_kDa * 1e3  # g/mol

    enzyme_activity_g_per_uL = enzyme_u_per_uL / specific_activity_U_per_g  # g/µL
    enzyme_mol_per_uL = enzyme_activity_g_per_uL / mw_g_per_mol  # mol/µL
    enzyme_activity_uM = enzyme_mol_per_uL * 1e12  # µM

    return enzyme_activity_uM

# -------------------------
# core ODEs
# -------------------------
def rca_ode(t, y, params, k):
    """
    ODE system (internal µM units).
    D: DNA template product (ng/uL)
    params: dictionary with fields already converted to µM
      'T' : total template ng/ul (free + bound)
      'P' : free primer µM
      'E' : enzyme "concentration" in µM 
      'N' : dNTP concentration in µM
      'Mg': Mg concentration in µM
    k: kinetics dict
      'Km_N', 'Km_Mg' : µM
      'k*' : µM^-2 s^-1
    Returns dD_dt
    """
    # state vector: y = [D, N, E]
    D = y[0]
    N = y[1]
    E = y[2]

    T = params['T']   # ng/uL (template)
    P = params['P']   # µM (primers)
    Mg = params['Mg'] # µM

    # kinetics / helpers
    Km_N = k.get('Km_N', 25.0)
    Km_Mg = k.get('Km_Mg', 2000.0)
    k_syn = k.get('k*', 1.0)

    # additional effects
    s_N = k.get('s_N', 1.0)       # stoichiometric factor: µM dNTP consumed per unit of D formed
    k_deg = k.get('k_deg', 1e-5)  # enzyme first-order decay (s^-1)

    # ensure non-negative concentrations in rate expressions
    N_pos = max(N, 0.0)
    E_pos = max(E, 0.0)

    # base synthesis rate (same structure as before)
    r_syn = k_syn * E_pos * T * P * (N_pos / (Km_N + N_pos)) * (Mg / (Km_Mg + Mg))

    # ODEs: product formation, dNTP depletion, enzyme decay
    dD_dt = r_syn
    dN_dt = -s_N * r_syn
    dE_dt = -k_deg * E_pos

    return [dD_dt, dN_dt, dE_dt]

def convert_params_for_solver(raw_params):
    """
    raw_params: dict with param_ranges as read from parameters.py
    Returns dict with all concentrations in µM, except 'T' which is in ng/µL
    """
    return {
        'T': raw_params['template'],
        'P': raw_params['primers'],
        'E': convert_enzyme_UuL_to_uM(raw_params['polymerase']),
        'N': raw_params['dNTPs'] * 1e3,  # mM to µM
        'Mg': raw_params['Mg2'] * 1e3,   # mM to µM
    }

def solve_rca(raw_params, t_max=16*3600, k=None, D0=0.0):
    """
    Solve single-ODE RCA from t=0 to t=t_max seconds.
    Returns:
        sol: ODE solution (t, D in µM)
    """
    params = convert_params_for_solver(raw_params)

    # default kinetics
    if k is None:
        k = {
            'Km_N': 25,      # µM
            'Km_Mg': 2000,   # µM
            'k*': 1,       # µM^-2 s^-1
        }

    # initial state: D, N, E
    y0 = [D0, params['N'], params['E']]

    t_span = (0.0, float(t_max))
    sol = solve_ivp(rca_ode, t_span, y0, args=(params, k),
                    dense_output=True, atol=1e-9, rtol=1e-6)

    # convert D from ng/µL to µg/µL for backward compatibility with plotting
    sol.y[0] = sol.y[0] / 1000.0

    return sol
