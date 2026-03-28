# parameter ranges in units to be displayed
param_ranges = {
    'template': (0, 20, 1),    # circular DNA template in ng/µL
    'polymerase': (0, 2, 0.1),   # phi29 DNA polymerase in U/µL
    'primers': (0, 100, 1),       # µM
    'dNTPs': (0, 20, 1),        # mM
    'Mg2': (0, 20, 1),          # mM
    'time': (0, 4*60, 60),       # range of 0 to 4 hours
}

# default parameters 
default_params = {
    'template': 10,
    'polymerase': 0.5,
    'primers': 25,
    'dNTPs': 5,
    'Mg2': 15,
    'time': 180,  # in minutes
}
