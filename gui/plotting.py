import plotly.graph_objects as go


def _base_figure(x, y):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', line={'width': 2}))
    # compact margins to maximize plotting area; enable light grid lines
    fig.update_layout(
        margin=dict(l=12, r=8, t=28, b=30),
        xaxis=dict(showgrid=True, gridcolor='#e6e6e6', zeroline=False),
        yaxis=dict(showgrid=True, gridcolor='#e6e6e6', zeroline=False),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    return fig


def _annotate_param(fig, label):
    # Small title above each subplot showing the parameter name
    fig.update_layout(annotations=[dict(
        text=label,
        x=0.5,
        y=1.05,
        xref='paper',
        yref='paper',
        showarrow=False,
        font=dict(size=12, color='#111')
    )])
    return fig


def plot_yield_vs_time(sol):
    t = sol.t / 3600.0  # hours
    y = sol.y[0]
    fig = _base_figure(t, y)
    # annotate parameter above plot
    fig = _annotate_param(fig, 'Time (h)')
    # hide axis titles to allow shared left label in layout
    fig.update_xaxes(title_text=None)
    fig.update_yaxes(title_text='RCA Product Yield (µg/µL)')
    return fig


def plot_yield_vs_param(param_name, param_values, yields):
    fig = _base_figure(param_values, yields)
    fig = _annotate_param(fig, param_name)
    fig.update_xaxes(title_text=None)

    if param_name == 'Primers (µM)':
        fig.update_yaxes(title_text='RCA Product Yield (µg/µL)')
    else:   
        fig.update_yaxes(title_text=None)
    return fig
