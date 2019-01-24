"""GEM Viewer.


python viewer.py "/media/tylerbiggs/genomic/data/gem_explore_test/" \
    "/media/tylerbiggs/genomic/data/sample_annots.txt"

"""


# ----------------------------------------------------------------------------
# General imports.
# ----------------------------------------------------------------------------
import os
# import sys
import argparse
import glob
import logging

import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc

import plotly.graph_objs as go

import pandas as pd


# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("source", action="store")
parser.add_argument("labels", action="store")


def build_options(options: list) -> list:
    return [dict(label=item, value=item) for item in options]


def slider_constructor(key):
    pass
    # values = UMAP_ARGS[key]
    # marks = {i: i for i in values}
    # return dcc.Slider(id=f"slider-{key}", marks=marks, min=np.min(values),
    # max=np.max(values), value=np.min(values))


def build_options(options: list) -> list:
    return [dict(label=item, value=item) for item in options]


def slider_constructor(name, values):
    marks = {i: str(i) for i in values}
    return html.Div(style={'margin': '25px 5px 30px 0px'},
                    children=[dcc.Slider(id=f"slider-{name}",
                                         marks=marks,
                                         min=min(values),
                                         max=max(values),
                                         step=None,
                                         value=values[0])])


class Prism:
    gem_globs = dict(
        gem_names="*",
        normalizations="*/*",
        reductions="*/*/*",
        tsne_pca_dims="*/*/tsne/*",
        tsne_perplexity="*/*/tsne/*/*/*",
        tsne_learning_rate="*/*/tsne/*/*/*/*",
        tsne_iterations="*/*/tsne/*/*/*/*/*",
        umap_n_neighbors="*/*/umap/*",
        umap_min_dist="*/*/umap/*/*",
        umap_metric="*/*/umap/*/*/*/*",
    )

    def __init__(self, gem_dir, label_dir):
        self.gem = gem_dir
        self.labels = label_dir
        self.label_df = pd.read_table(self.labels, index_col=0)

    def glob_path_names(self, path):
        return list(set(os.path.basename(x)
                        for x in glob.glob(os.path.join(self.gem, path))))

    @property
    def gem_dict(self):
        return {k: self.glob_path_names(v) for k, v in Prism.gem_globs.items()}

    @property
    def columns(self):
        output = {}
        output['all'] = sorted(self.label_df.columns)
        output['discrete'] = [col for col in output['all']
                              if self.label_df[col].dtype == object]
        output['continuous'] = [col for col in output['all']
                                if col not in output['discrete']]
        output['quantileable'] = [col for col in output['continuous']
                                  if len(self.label_df[col].unique()) > 20]
        return output


# Parse the given arguments.
args = parser.parse_args()

# Create a new Prism object.
prism = Prism(gem_dir=args.source, label_dir=args.labels)

# Create the server application.
app = dash.Dash('gem-viewer')
server = app.server
app.config['suppress_callback_exceptions'] = True
external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "//fonts.googleapis.com/css?family=Raleway:400,300,600",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
    "https://cdn.rawgit.com/plotly/dash-tsne/master/custom_styles.css",
    "https://cdn.rawgit.com/plotly/dash-app-stylesheets/2cc54b8c03f4126569a3440aae611bbef1d7a5dd/stylesheet.css"
]

for css in external_css:
    app.css.append_css({"external_url": css})

umap_controls = [
    "# of Neighbors",
    slider_constructor('umap-n_neighbors',
                       prism.gem_dict['umap_n_neighbors']),
    "Min Distance",
    slider_constructor('umap-min_distance',
                       prism.gem_dict['umap_min_dist']),
    "Distance Metric",
    dcc.Dropdown(id='dropdown-umap-metric',
                 options=prism.gem_dict['umap_metric'],
                 value=prism.gem_dict['umap_metric'][0])]


tsne_controls = [
    "PCA Dimensions Input",
    slider_constructor('tsne_pca_dims',
                       prism.gem_dict['tsne_pca_dims']),
    "Perplexity",
    slider_constructor('tsne_perplexity',
                       prism.gem_dict['tsne_perplexity']),
    "Learning Rate",
    slider_constructor('tsne_learning_rate',
                       prism.gem_dict['tsne_learning_rate']),
    "Iterations",
    slider_constructor('tsne_iterations',
                       prism.gem_dict['tsne_iterations'])]

pca_controls = [
    "Select X PCA Dimension",
    slider_constructor('tsne_pca_dims', range(20)),
    "Select Y PCA Dimension",
    slider_constructor('tsne_pca_dims', range(20)),
    "Select Z PCA Dimension",
    slider_constructor('tsne_pca_dims', range(20)),
]


def base_layout(plot_id):
    return [
        # Create the title div.
        html.Div(className="row", children=[html.H1("GEM Explorer")]),
        # Create the div to hold the plot.
        html.Div(className="row", children=[dcc.Graph(id=plot_id)]),
        # Create the control div.
        html.Div(children=[
            "Projection Dimensions",
            dcc.RadioItems(id='radio-projection',
                           labelStyle={'display': 'inline-block'},
                           options=[{"label": "2D", "value": "2"},
                                    {"label": "3D", "value": "3"}],
                           value="3"),
            "Normalization Method",
            dcc.Dropdown(id='dropdown-normalizations',
                         options=build_options(prism.gem_dict['normalizations']),
                         value=prism.gem_dict['normalizations'][0]),
            "Projection Method",
            dcc.Dropdown(id='dropdown-reductions',
                         options=build_options(prism.gem_dict['reductions']),
                         value=prism.gem_dict['reductions'][0]),
            "Color By",
            dcc.Dropdown(id='dropdown-color_by',
                         options=prism.columns['discrete'] + prism.columns['quantileable'],
                         value=prism.columns['discrete'][0]),
            "Size By",
            dcc.Dropdown(id='dropdown-size_by',
                         options=prism.columns['continuous']),
        ]),
    ]


# Define the layout of the application.
app.layout = html.Div([
    html.H1('GEM PRISM Viewer'),
    dcc.Tabs(id="tabs-main", value='tab-umap', children=[
        dcc.Tab(label='Tab One', value='tab-umap'),
        dcc.Tab(label='Tab Two', value='tab-tsne'),
        dcc.Tab(label='Tab Two', value='tab-pca'),
    ]),
    html.Div(id='tabs-content')
])


@app.callback(Output('tabs-content', 'children'),
              [Input('tabs-main', 'value')])
def render_content(tab):
    options = {'tab-umap': umap_controls,
               'tab-tsne': tsne_controls,
               'tab-pca': pca_controls}
    return html.Div(className="container",
                    style={'width': '90%', 'max-width': 'none',
                           'font-size': '1.5rem', 'padding': '10px 30px'},
                    children=base_layout(f"figure-{tab}") + options[tab])


def generate_3d_plot(groups, fig_layout, marker_dict):
    plot_data = []
    for name, values in groups:
        scatter = go.Scatter3d(
            name=name, x=values['x'], y=values['y'], z=values['z'],
            mode='markers', marker=marker_dict,
            text=list(values[[col for col in values.columns
                              if col not in ['x', 'y', 'z']]].to_records()))
        plot_data.append(scatter)
    return go.Figure(data=plot_data, layout=fig_layout)


def generate_2d_plot(groups, fig_layout, marker_dict):
    plot_data = []
    for name, values in groups:
        scatter = go.Scattergl(
            name=name, x=values['x'], y=values['y'], mode='markers',
            marker=marker_dict,
            text=list(values[[col for col in values.columns
                              if col not in ['x', 'y', 'z']]].to_records()))
        plot_data.append(scatter)
    return go.Figure(data=plot_data, layout=fig_layout)


@app.callback(
        Output('figure-tab-umap', 'figure'),
        [Input('slider-umap-n_neighbors', 'value'),
            Input('slider-umap-min_distance', 'value'),
            Input('dropdown-umap-metric', 'value'),
            Input('radio-projection', 'value'),
            Input('dropdown-normalizations', 'value'),
            Input('dropdown-reductions', 'value'),
            Input('dropdown-color_by', 'value'),
            Input('dropdown-size_by', 'value')])
def update_main_plot(n_neighbors, min_distance, metric, projection,
                     normalizations, reductions, color_by, size_by):

    logging.debug("UMAP callback activated.")
    logging.debug("arguemnts", n_neighbors, min_distance, metric, projection,
                         normalizations, reductions, color_by, size_by)
    # Ensure string conversion for path generation works. Ad hoc.
    if min_distance == 0:
        min_distance = "0.0"

    # Projection dictionary to avoid `if` statements.
    projections = {"2": {"func": generate_2d_plot, "columns": ["x", "y"]},
                   "3": {"func": generate_3d_plot, "columns": ["x", "y", "z"]}}

    # Plot and figure layout.
    axes = dict(title='', showgrid=True,
                zeroline=False, showticklabels=False)
    fig_layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0), autosize=True,
                           # plot_bgcolor='#252525', paper_bgcolor='#202020',
                           height=800,
                           font=dict(size=18, color='#CCCCCC'),
                           legend=dict(
                               bgcolor='rgba(255, 255, 255, 0.00)'),
                           scene=dict(xaxis=axes, yaxis=axes, zaxis=axes))

    # Try to load the selected embedding.
    # Create the path.
    path = os.path.join(
        prism.gem, normalizations, 'umap', n_neighbors, min_distance,
        projection, metric,
        f"umap_{n_neighbors}_{min_distance}_{projection}_{metric}.csv"),

    logging.debug(f"Attempting to load file at:\n\t{path}")

    try:
        xdf = pd.read_csv(path, names=projections[projection]["columns"])
        df = pd.concat([xdf, prism.label_df], axis=1)
        logging.debug(f"Data frame loaded for with columns:\n\t{df.columns}")
    except FileNotFoundError as error:
        logging.warning(
            f'Unable to load file:\n\t{path}.\nDue to\n\t{error}')
        return go.Figure(layout=fig_layout)

    # Create the size and marker dictionaries.
    sizes = df[size_by].astype('float64').values if size_by else 6
    marker_dict = dict(size=sizes, symbol='circle',
                       sizeref=(2.0 * max(sizes) / (3. ** 2)
                                ) if size_by else None,
                       sizemin=1, line=dict(width=0), opacity=0.8)
    # color=embedding_df[color_factor] if color_factor else None,
    # colorscale='Viridis' if color_factor else None)

    groups = df.groupby(color_by) if color_by else [("All", df), ]
    return projections[projection]["func"](groups, fig_layout, marker_dict)

#
# @app.callback(
#         Output('figure-tab-umap', 'figure'),
#         [Input('slider-umap-n_neighbors', 'value'),
#             Input('slider-umap-min_distance', 'value'),
#             Input('dropdown-umap-metric', 'value'),
#             Input('radio-projection', 'value'),
#             Input('dropdown-normalizations', 'value'),
#             Input('dropdown-reductions', 'value'),
#             Input('dropdown-color_by', 'value'),
#             Input('dropdown-size_by', 'value')])
# def update_main_plot(n_neighbors, min_distance, metric,
#                      projection, normalizations,
#                      reductions, color_by, size_by):
#
#     logging.debug("Main callback activated.")
#     # Ensure string conversion for path generation works. Ad hoc.
#     if min_distance == 0:
#         min_distance = "0.0"
#
#     # Projection dictionary to avoid `if` statements.
#     projections = {"2": {"func": generate_2d_plot, "columns": ["x", "y"]},
#                    "3": {"func": generate_3d_plot, "columns": ["x", "y", "z"]}}
#
#     # Plot and figure layout.
#     axes = dict(title='', showgrid=True,
#                 zeroline=False, showticklabels=False)
#     fig_layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0), autosize=True,
#                            # plot_bgcolor='#252525', paper_bgcolor='#202020',
#                            height=800,
#                            font=dict(size=18, color='#CCCCCC'),
#                            legend=dict(
#                                bgcolor='rgba(255, 255, 255, 0.00)'),
#                            scene=dict(xaxis=axes, yaxis=axes, zaxis=axes))
#
#     # Try to load the selected embedding.
#     # Create the path.
#     paths = {'tsne': os.path.join(prism.gem, normalizations, 'tsne',
#                                   tsne_pca_dims, projection,
#                                   tsne_perplexity, tsne_learning_rate,
#                                   tsne_iterations, "tsne.csv"),
#              'umap': os.path.join(prism.gem, normalizations, 'umap',
#                                   n_neighbors, min_distance, projection,
#                                   metric, f"umap_{n_neighbors}_{min_distance}_{projection}_{metric}.csv"),
#              'pca': os.path.join(prism.gem, normalizations, 'pca',
#                                  n_neighbors, min_distance, "pca.csv")}
#
#     path = paths[reductions]
#
#     logging.debug(f"Attempting to load file at:\n\t{path}")
#     try:
#         xdf = pd.read_csv(path, names=projections[projection]["columns"])
#         df = pd.concat([xdf, prism.label_df], axis=1)
#         logging.debug(f"Data frame loaded for with columns:\n\t{df.columns}")
#     except FileNotFoundError as error:
#         logging.warning(
#             f'Unable to load file:\n\t{path}.\nDue to\n\t{error}')
#         return go.Figure(layout=fig_layout)
#
#     # Create the size and marker dictionaries.
#     sizes = df[size_by].astype('float64').values if size_by else 6
#     marker_dict = dict(size=sizes, symbol='circle',
#                        sizeref=(2.0 * max(sizes) / (3. ** 2)
#                                 ) if size_by else None,
#                        sizemin=1, line=dict(width=0), opacity=0.8)
#     # color=embedding_df[color_factor] if color_factor else None,
#     # colorscale='Viridis' if color_factor else None)
#
#     groups = df.groupby(color_by) if color_by else [("All", df), ]
#     return projections[projection]["func"](groups, fig_layout, marker_dict)
#
#
# @app.callback(
#         Output('figure-tab-umap', 'figure'),
#         [Input('slider-umap-n_neighbors', 'value'),
#             Input('slider-umap-min_distance', 'value'),
#             Input('dropdown-umap-metric', 'value'),
#             # Input('slider-tsne_pca_dims', 'value'),
#             # Input('slider-tsne_perplexity', 'value'),
#             # Input('slider-tsne_learning_rate', 'value'),
#             # Input('slider-tsne_iterations', 'value'),
#             Input('radio-projection', 'value'),
#             Input('dropdown-normalizations', 'value'),
#             Input('dropdown-reductions', 'value'),
#             Input('dropdown-color_by', 'value'),
#             Input('dropdown-size_by', 'value')])
# def update_main_plot(n_neighbors, min_distance, metric, tsne_pca_dims,
#                      tsne_perplexity, tsne_learning_rate,
#                      tsne_iterations, projection, normalizations,
#                      reductions, color_by, size_by):
#
#     logging.debug("Main callback activated.")
#     # Ensure string conversion for path generation works. Ad hoc.
#     if min_distance == 0:
#         min_distance = "0.0"
#
#     # Projection dictionary to avoid `if` statements.
#     projections = {"2": {"func": generate_2d_plot, "columns": ["x", "y"]},
#                    "3": {"func": generate_3d_plot, "columns": ["x", "y", "z"]}}
#
#     # Plot and figure layout.
#     axes = dict(title='', showgrid=True,
#                 zeroline=False, showticklabels=False)
#     fig_layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0), autosize=True,
#                            # plot_bgcolor='#252525', paper_bgcolor='#202020',
#                            height=800,
#                            font=dict(size=18, color='#CCCCCC'),
#                            legend=dict(
#                                bgcolor='rgba(255, 255, 255, 0.00)'),
#                            scene=dict(xaxis=axes, yaxis=axes, zaxis=axes))
#
#     # Try to load the selected embedding.
#     # Create the path.
#     paths = {'tsne': os.path.join(prism.gem, normalizations, 'tsne',
#                                   tsne_pca_dims, projection,
#                                   tsne_perplexity, tsne_learning_rate,
#                                   tsne_iterations, "tsne.csv"),
#              'umap': os.path.join(prism.gem, normalizations, 'umap',
#                                   n_neighbors, min_distance, projection,
#                                   metric, f"umap_{n_neighbors}_{min_distance}_{projection}_{metric}.csv"),
#              'pca': os.path.join(prism.gem, normalizations, 'pca',
#                                  n_neighbors, min_distance, "pca.csv")}
#
#     path = paths[reductions]
#
#     logging.debug(f"Attempting to load file at:\n\t{path}")
#     try:
#         xdf = pd.read_csv(path, names=projections[projection]["columns"])
#         df = pd.concat([xdf, prism.label_df], axis=1)
#         logging.debug(f"Data frame loaded for with columns:\n\t{df.columns}")
#     except FileNotFoundError as error:
#         logging.warning(
#             f'Unable to load file:\n\t{path}.\nDue to\n\t{error}')
#         return go.Figure(layout=fig_layout)
#
#     # Create the size and marker dictionaries.
#     sizes = df[size_by].astype('float64').values if size_by else 6
#     marker_dict = dict(size=sizes, symbol='circle',
#                        sizeref=(2.0 * max(sizes) / (3. ** 2)
#                                 ) if size_by else None,
#                        sizemin=1, line=dict(width=0), opacity=0.8)
#     # color=embedding_df[color_factor] if color_factor else None,
#     # colorscale='Viridis' if color_factor else None)
#
#     groups = df.groupby(color_by) if color_by else [("All", df), ]
#     return projections[projection]["func"](groups, fig_layout, marker_dict)


if __name__ == '__main__':
    app.run_server(debug=True)
