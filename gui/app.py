import dash
from dash import dcc, html, no_update
from dash.dependencies import Input, Output, State
from model.odes import solve_rca
from model.parameters import default_params, param_ranges
from gui.plotting import plot_yield_vs_time, plot_yield_vs_param
import numpy as np

# initialize Dash app
app = dash.Dash(__name__)
server = app.server  # for deployment

# ---------------------------
# Layout
# ---------------------------
app.layout = html.Div([
    html.Div(
        html.H1("RCA Optimization Dashboard", style={'textAlign': 'center', 'marginBottom': '0', 'marginTop': '4px'}),
        style={
            'backgroundColor': "#edeff0",
            'padding': '10px 10px',
            'borderRadius': '8px',
            'marginBottom': '12px'
        }
    ),

    html.Div([
        html.H2("Parameters"),

        html.Label('Template (ng/µL)'),
        dcc.Slider(
            id='template-slider',
            min=param_ranges['template'][0],
            max=param_ranges['template'][1],
            step=param_ranges['template'][2],
            value=default_params['template'],
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

        html.Label('Polymerase (U/µL)'),
        dcc.Slider(
            id='polymerase-slider',
            min=param_ranges['polymerase'][0],
            max=param_ranges['polymerase'][1],
            step=param_ranges['polymerase'][2],
            value=default_params['polymerase'],
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

        html.Label('Primers (µM)'),
        dcc.Slider(
            id='primers-slider',
            min=param_ranges['primers'][0],
            max=param_ranges['primers'][1],
            step=param_ranges['primers'][2],
            value=default_params['primers'],
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

        html.Label('dNTPs (mM)'),
        dcc.Slider(
            id='dntps-slider',
            min=param_ranges['dNTPs'][0],
            max=param_ranges['dNTPs'][1],
            step=param_ranges['dNTPs'][2],
            value=default_params['dNTPs'],
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

        html.Label('Mg²⁺ (mM)'),
        dcc.Slider(
            id='mg-slider',
            min=param_ranges['Mg2'][0],
            max=param_ranges['Mg2'][1],
            step=param_ranges['Mg2'][2],
            value=default_params['Mg2'],
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

        html.Label('Incubation time (h)'),
        dcc.Slider(
            id='time-slider',
            min=0,
            max=4,
            step=0.1,
            value=4,
            tooltip={"placement": "bottom", "always_visible": True},
            marks=None
        ),

    ], style={'width': '30%', 'float': 'left', 'padding': '10px'}),

    html.Div([
        html.Div([
            html.H3('RCA Product yield (µg/µL) vs.', style={'textAlign': 'center', 'gridColumn': '1 / -1', 'marginTop': '0'}),

            # 3 columns, 2 rows grid of plots
            html.Div([
                html.Div(dcc.Graph(id='yield-time-plot', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'}),
                html.Div(dcc.Graph(id='yield-vs-template', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'}),
                html.Div(dcc.Graph(id='yield-vs-polymerase', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'}),
                html.Div(dcc.Graph(id='yield-vs-primers', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'}),
                html.Div(dcc.Graph(id='yield-vs-dntps', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'}),
                html.Div(dcc.Graph(id='yield-vs-mg', style={'height': '240px'}),
                         style={'border': '1px solid #e6e6e6', 'borderRadius': '6px', 'padding': '6px', 'background': 'transparent'})
            ], style={
                'display': 'grid',
                'gridTemplateColumns': 'repeat(3, 1fr)',
                'gridGap': '12px'
            })

        ], style={
            'width': '100%',
            'boxSizing': 'border-box',
        })

    ], style={
        'width': '65%',
        'float': 'right',
        'padding': '12px',
        'display': 'flex',
        'alignItems': 'flex-start',
        'boxSizing': 'border-box',
        'backgroundColor': "#edeff0",
        'borderRadius': '8px'
    })],

    style={'backgroundColor': '#f2f2f2', 'minHeight': '100vh', 'padding': '12px'}
)

# ---------------------------
# Callback
# ---------------------------
@app.callback(
    [
        Output('yield-time-plot', 'figure'),
        Output('yield-vs-template', 'figure'),
        Output('yield-vs-polymerase', 'figure'),
        Output('yield-vs-primers', 'figure'),
        Output('yield-vs-dntps', 'figure'),
        Output('yield-vs-mg', 'figure')
    ],
    [
        Input('template-slider', 'value'),
        Input('polymerase-slider', 'value'),
        Input('primers-slider', 'value'),
        Input('dntps-slider', 'value'),
        Input('mg-slider', 'value'),
        Input('time-slider', 'value')
    ]
)
def update_all_figures(template, polymerase, primers, dntps, mg, time_h):
    # base params dictionary
    params = {
        'template': template,
        'polymerase': polymerase,
        'primers': primers,
        'dNTPs': dntps,
        'Mg2': mg,
    }

    t_max = time_h * 3600  # seconds

    # solve RCA ODE
    sol = solve_rca(params, t_max=t_max, D0=0.0)

    # main time-course figure
    fig_time = plot_yield_vs_time(sol)

    # Sweep function for parameter sensitivity
    n_points = 20
    def sweep_param(param_name, param_range):
        values = np.linspace(param_range[0], param_range[1], n_points)
        yields = []
        for val in values:
            sweep_params = params.copy()
            sweep_params[param_name] = val
            sol_sweep = solve_rca(sweep_params, t_max=t_max, D0=0.0)
            yields.append(sol_sweep.y[0, -1])  # final yield in µg/µL
        return values, yields

    template_vals, template_yields = sweep_param('template', param_ranges['template'])
    polymerase_vals, polymerase_yields = sweep_param('polymerase', param_ranges['polymerase'])
    primers_vals, primers_yields = sweep_param('primers', param_ranges['primers'])
    dntps_vals, dntps_yields = sweep_param('dNTPs', param_ranges['dNTPs'])
    mg_vals, mg_yields = sweep_param('Mg2', param_ranges['Mg2'])

    fig_template = plot_yield_vs_param('Template (ng/µL)', template_vals, template_yields)
    fig_polymerase = plot_yield_vs_param('Polymerase (U/µL)', polymerase_vals, polymerase_yields)
    fig_primers = plot_yield_vs_param('Primers (µM)', primers_vals, primers_yields)
    fig_dntps = plot_yield_vs_param('dNTPs (mM)', dntps_vals, dntps_yields)
    fig_mg = plot_yield_vs_param('Mg²⁺ (mM)', mg_vals, mg_yields)

    return fig_time, fig_template, fig_polymerase, fig_primers, fig_dntps, fig_mg
