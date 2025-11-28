# parameter ranges in units to be displayed
param_ranges = {
    'template': (0, 10, 0.01),    # circular DNA template in ng/µL
    'polymerase': (0, 1, 0.01),   # phi29 DNA polymerase in U/µL
    'primers': (0, 5, 0.1),       # µM
    'dNTPs': (0, 20, 0.1),        # mM
    'Mg2': (0, 10, 0.1),          # mM
    'time': (0, 4*60, 60),       # range of 0 to 4 hours
}

# default parameters 
default_params = {
    'template': 1,
    'polymerase': 0.5,
    'primers': 0.5,
    'dNTPs': 10,
    'Mg2': 5,
    'time': 120,  # in minutes
}
