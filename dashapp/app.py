import json
import dash 
from dash import html, dcc, Input, Output, State, ALL, MATCH
import dash_bootstrap_components as dbc
import os
import polars as pl
from polars import DataFrame
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

#TODO: Borrar sys cuando se arme el package completo y se publique
import sys
sys.path.append('/home/faivel/Documents/Projects/Tesis/kripr/')

from kripr import get_read_length_histogram, BioString, DNAString, DNAStringSet, GTFobject,BAMhandler,FastaHandler,GTFhandler
import time

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP])

# region auxiliar functions

class Timer:
    def __init__(self, name) -> None:
        self.name = name
        self.start = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self,exception_type, exception_value, exception_traceback):
        print(f'{self.name} has ended in {time.time() - self.start} seconds')



def process_bam_histogram(histogram_as_dict: dict[int,int]):
    max_key = max(histogram_as_dict.keys())

    keys_list = []
    values_list = []

    for i in range(max_key+1):
        keys_list.append(i)
        if i in histogram_as_dict:
            values_list.append(histogram_as_dict[i])
        else:
            values_list.append(0)
    return [keys_list, values_list]

def filter_by_rna_coverage(regions, rna_cov, miminum_cov):
    filtered_regions = []
    for region in regions:
        region_cov = rna_cov[region['start']:region['end']]
        if np.mean(np.array(region_cov)) >= miminum_cov:
            filtered_regions.append(region)

    return filtered_regions

def get_regions_from_outliers(l: list[list[int,int]], min_distance: int, min_len: int):
    if len(l) == 0:
        return []
    regions = [[l[0]]]
    last_outlier = l[0]
    for outlier in l[1:]:
        if outlier[0] - last_outlier[0] < min_distance:
            regions[-1].append(outlier)
        else:
            regions.append([outlier])
        last_outlier = outlier
    filtered_regions = [region for region in regions if len(region) >= min_len ]

    parsed_regions =[]
    region = {'start': float('inf'), 'end': float('-inf'), 'values': []}
    for reg in filtered_regions:
        for r in reg:
            if r[0] < region['start']:
                region['start'] = r[0]
            if r[0] > region['end']:
                region['end'] = r[0]
            region['values'].append(r[1])
        parsed_regions.append(region)
        region = {'start': float('inf'), 'end': float('-inf'), 'values': []}

    return parsed_regions
    

def detect_pauses(ts: list[int], min_distance, min_len):
    l = np.array(ts)
    [q1,q3] = np.quantile(l, [.25,.75])
    inter_quantile = q3-q1
    max = q3+inter_quantile
    return get_regions_from_outliers([[i, l[i]] for i in range(len(l)) if l[i] > max ],min_distance, min_len)

def detect_rolling_pauses(ts: list[int], ws: int, step: int, min_distance: int, min_len: int):
    outliers = []
    outliers_indexes = []
    for i in range(len(ts)-ws):
        l = np.array(ts[i:i+ws])
        [q1,q3] = np.quantile(l, [.25,.75])
        inter_quantile = q3-q1
        max = q3+inter_quantile
        for b in l:
            if b > max and i not in outliers_indexes:
                outliers_indexes.append(i)
                outliers.append([i, ts[i]])
    return get_regions_from_outliers(sorted(outliers, key= lambda x: x[0]),min_distance,min_len)

def detect_rolling_z_score_pauses(ts: list[int], ws: int, step: int, min_distance: int, min_len: int, z_score: float = 1.5):
    outliers = []
    outliers_indexes = []
    for i in range(0, len(ts)-ws, step):
        l = np.array(ts[i:i+ws])
        mean = np.mean(l)
        std = np.std(l)
        for index in range(i, i+ws):
            if (ts[index]-mean)/std > z_score and index not in outliers_indexes:
                outliers_indexes.append(index)
                outliers.append([index, ts[index]])
    return get_regions_from_outliers(sorted(outliers, key= lambda x: x[0]),min_distance,min_len)

def detect_z_score_pauses(ts: list[int], min_distance: int, min_len: int, z_score: float = 1.5):
    l = np.array(ts)
    mean = np.mean(l)
    std = np.std(l)
    return get_regions_from_outliers([[i, l[i]] for i in range(len(l)) if (l[i]-mean)/std > z_score ],min_distance, min_len)

def detect_rolling_percentile_pauses(ts: list[int], ws: int, step:int, min_distance: int, min_len: int, percentile: int = 85):
    outliers: list[list[int,int]] = []
    outliers_indexes = []
    for i in range(0, len(ts)-ws, step):
        l = np.array(ts[i:i+ws])
        quantile = np.quantile(l, percentile)
        for index in range(i, i+ws):
            if ts[index] > quantile and index not in outliers_indexes:
                outliers_indexes.append(index)
                outliers.append([index, ts[index]])
    return get_regions_from_outliers(sorted(outliers, key= lambda x: x[0]),min_distance,min_len)

def detect_percentile_pauses(ts: list[int], min_distance: int, min_len: int, percentile: float = .85):
    l = np.array(ts)
    max = np.quantile(l, percentile)
    return get_regions_from_outliers([[i, l[i]] for i in range(len(l)) if l[i] > max ],min_distance, min_len)
                
def get_pause_to_graph(pauses, xaxis_offset = 0):
    traces = []
    for p in pauses:
        trace = [[],[]]
        trace[0].append(p['start'])        
        trace[1].append(p['values'][0])
        if p['start'] != p['end']:
            trace[0].append(p['end'])
            trace[1].append(p['values'][-1])
        traces.append(trace)
    return [[[t+xaxis_offset for t in x[0]], x[1]] for x in traces]

#TODO: cambiar methods a una clase de contantes
def get_normalize_rpf_over_rna(cov_rpf: list[int], cov_rna: list[int], offset: int=1, window_size:int = 0, method: str= 'Offset'):
    if method == "Offset":
        new_rna = [c if c > 0 else offset for c in cov_rna]
    else:
        new_rna = [int(np.mean(cov_rna[i-window_size: i+window_size])) if cov_rna[i] == 0 else cov_rna[i]for i in range(window_size, len(cov_rna)-window_size)]
        new_rna = [c for c in cov_rna[:window_size]] + new_rna + [c for c in cov_rna[-window_size:]]
        new_rna = [c if c > 0 else offset for c in new_rna]
    return [cov_rpf[i]/new_rna[i] for i in range(len(cov_rpf))]

def get_second_elements(l: list[str,str])->list[int]:
    new_list = []
    for elem in l:
        new_list.append(int(elem[1]))

    return new_list


def get_feature_sequence_and_coverage(df: DataFrame, dna: DNAStringSet, bam: BAMhandler):
    iter_df = df.iterrows(named=True)
    feature_sequences = {}
    for row in iter_df:
        feature_seq = []
        feature_depth = []
        ref_seq = dna.get_sequence(row.seqname)
        for index in range(len(row.start)):
            seq = ref_seq.get_sub_sequence(row.start[index], row.end[index])
            if row.strand[index] == '-':
                seq = seq.reverse_complement()
            # print(row.seqname, row.start[index], row.end[index])
            cov = bam.region_coverage(row.seqname, row.start[index], row.end[index])
            feature_depth += get_second_elements(cov)
            feature_seq.append(seq.get_as_string())
        feature_sequences[row.transcript_id] = [''.join(feature_seq), feature_depth]
    return feature_sequences

def zip_lists_to_dict(keys: list[str], values: list[any]) -> dict:
    if len(keys) != len(values):
        raise IndexError("The lenght of the list must be the same")
    new_dict = {}
    for i in range(len(keys)):
        new_dict[keys[i]] = values[i]
    return new_dict

# endregion


# region Constants and Hardcoded
DATA_PATH = './data/'
RNA_PATH = './data/new_rna.parquet'
RPF_PATH = './data/new_rpf.parquet'
CDS_PATH = './data/cds_by_frame.json'
# endregion


# region variable preparation
selectable_files = os.listdir(DATA_PATH)

df_rpf: DataFrame = pl.read_parquet(RPF_PATH)
df_rna: DataFrame = pl.read_parquet(RNA_PATH)

transcript_ids = df_rpf.select('transcript_id').to_series().to_list()
transcript_to_genes = (df_rna
                .lazy()
                .groupby(by='transcript_id')
                .agg(
                    [
                        pl.col('gene_id').first()
                    ]
                )
            ).collect().to_dict(as_series=False)
        
transcript_to_genes = zip_lists_to_dict(transcript_to_genes['transcript_id'],transcript_to_genes['gene_id'])


bam_files_info = {}

cdss: DataFrame = pl.read_json(CDS_PATH)

protein_ids = cdss.select('protein_ids').unique().to_series().to_list()
# endregion

def create_file_card(file_name):
    return dbc.ListGroupItem(file_name)
            

home_tab_content = dbc.Card(
    dbc.CardBody(
        [
            html.P("This is tab 1!", className="card-text"),
            dbc.Button("Click here", color="success"),
        ]
    ),
    className="mt-3",
)

data_selection_tab_content = html.Div(
    [
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.ListGroup(
                        [
                            dbc.ListGroupItem(file_name, id={"index": i})
                            for i, file_name in enumerate(selectable_files)
                        ],
                        style={
                            "maxHeight": "550px",  # set the maximum height of the ListGroup
                            "overflowY": "scroll",  # set the y-axis overflow behavior to scroll
                        }
                    )
                ),
                dbc.Col(
                    dbc.Spinner
                        (dbc.Card(
                            [
                                dbc.CardHeader(children="Files to Select", id='selected_file'),
                                dbc.CardBody(
                                    [
                                        html.H4(id='bam_size'),
                                        html.H4(id='read_count'),
                                        html.H4(id='mapped_reads'),
                                        html.H4(id='unmapped_reads'),
                                        dcc.Graph(id="bam_histogram_graph"),
                                    ]
                                )
                            ]
                        )), md=8
                )
            ]
        )
    ]
)

transcript_selection_params = html.Div(
    [
        html.H4("Transcript Selection"),
        html.Div([
            dbc.Row(
                [
                dbc.Label("Transcript ID",width="auto"),
                dbc.Col(dbc.Select(transcript_ids, value=transcript_ids[0], id='feature-select'), md=8)
                ]
            )
        ]),
        html.Div([
            dbc.Label("Coverage % per base (filter)"),
            dcc.Slider(0, 80, 1, value=10, id='coverage-filter-slider', marks=None,
                tooltip={"placement": "bottom", "always_visible": False})
        ]),
    ]
)

cds_selection_params = html.Div(
    [
        html.H4("CDS Selection"),
        html.Div([
            dbc.Row(
                [
                dbc.Label("Protein ID",width="auto"),
                # TODO lista de protein id para seleccionar
                dbc.Col(dbc.Select(protein_ids, value=protein_ids[0], id='cds-select'), md=8)
                ]
            ),
            dbc.Row(
                [
                dbc.Label("Frame",width="auto"),
                dbc.RadioItems(['0', '1', '2'], value='0', id='cds-frame', inline=True)
                ]
            ),
        ]),
        # html.Div([
        #     dbc.Label("Coverage % per base (filter)"),
        #     dcc.Slider(0, 80, 1, value=10, id='coverage-filter-slider', marks=None,
        #         tooltip={"placement": "bottom", "always_visible": False})
        # ]),
    ]
)

normalization_params = html.Div(
    [
        html.H4("RPF/RNA Normalization"),
        html.Div([
            dbc.Label("DESeq2 Normalization"),
            # dcc.RadioItems(['Offset', 'Mean'], 'Offset', id='method', inline=True),
            dbc.RadioItems(['Yes', 'No'], value='No', id='use_deseq2', inline=True)
        ]),
        html.Div([
            dbc.Label("Division Method"),
            # dcc.RadioItems(['Offset', 'Mean'], 'Offset', id='method', inline=True),
            dbc.RadioItems(['Offset', 'Mean'], value='Offset', id='zero-handler-method', inline=True)
        ]),
        html.Div([
            dbc.Row(
                [
                dbc.Label("Window Size",width="auto"),
                dbc.Col(dbc.Input(type="number", min=1, max=10, step=1, value=3, id = 'mean-window-size')),
                dbc.Label("RNA Offset",width="auto"),
                dbc.Col(dbc.Input(type="number", min=1, max=5, step=1, value=3, id='RNA-offset')),
                ]
            )
        ]),
    ]
)

pause_detection_params = html.Div(
    [
        html.H4("Pause Detection Parameters"),
        html.Div([
            dbc.Row(
                [
                    dbc.Label("Method", width="auto"),
                    dbc.RadioItems(['Boxplot', 'Percentile', 'Z-Score'], value='Boxplot', id='pause-detection-method', inline=True),
                    dbc.RadioItems(['Fixed', 'Rolling'], value='Fixed', id='pause-rolling', inline=True),
                    dbc.Label("Step",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=100, step=1, value=10, id="rolling-step")),
                    dbc.Label("Window",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=100, step=1, value=20, id="rolling-window")),
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Boxplot",width="auto"),
                ]
            ),
            dbc.Row(
                [    
                    
                    dbc.Label("Max Distance",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=20, step=1, value=2, id="max-distance-between-outliers")),
                    dbc.Label("Min Length",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=20, step=1, value=3, id="min-region-length")),
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Percentile", width="auto"),
                    dcc.Slider(1, 99, 1, value=80, id='percentile-slider', marks=None,
                    tooltip={"placement": "bottom", "always_visible": False})
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Z-Score", width="auto"),
                    dcc.Slider(0.5, 3.5, 0.1, value=1.0, id='z-score-slider', marks=None,
                    tooltip={"placement": "bottom", "always_visible": False})
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Filter Pause Region by RNA coverage", width="auto"),
                    dcc.Slider(0, 100, 1, value=0, id='rna-coverage-region-filter', marks=None,
                    tooltip={"placement": "bottom", "always_visible": False})
                ]
            ),

        ]),
    ]
)

cds_pause_detection_params = html.Div(
    [
        html.H4("Pause Detection Parameters"),
        html.Div([
            dbc.Row(
                [
                    dbc.Label("Method", width="auto"),
                    dbc.RadioItems(['Boxplot', 'Percentile', 'Z-Score'], value='Boxplot', id='cds-pause-detection-method', inline=True),
                    dbc.RadioItems(['Fixed', 'Rolling'], value='Fixed', id='cds-pause-rolling', inline=True),
                    dbc.Label("Step",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=100, step=1, value=10, id="cds-rolling-step")),
                    dbc.Label("Window",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=100, step=1, value=20, id="cds-rolling-window")),
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Boxplot",width="auto"),
                ]
            ),
            dbc.Row(
                [    
                    
                    dbc.Label("Max Distance",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=20, step=1, value=2, id="cds-max-distance-between-outliers")),
                    dbc.Label("Min Length",width="auto"),
                    dbc.Col(dbc.Input(type="number", min=1, max=20, step=1, value=3, id="cds-min-region-length")),
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Percentile", width="auto"),
                    dcc.Slider(1, 99, 1, value=80, id='cds-percentile-slider', marks=None,
                    tooltip={"placement": "bottom", "always_visible": False})
                ]
            ),
            dbc.Row(
                [
                    dbc.Label("Z-Score", width="auto"),
                    dcc.Slider(0.5, 3.5, 0.1, value=1.0, id='cds-z-score-slider', marks=None,
                    tooltip={"placement": "bottom", "always_visible": False})
                ]
            ),
        ]),
    ]
)

controls_card = dbc.Card(
    [
        transcript_selection_params,
        normalization_params,
        pause_detection_params,
    ],
    body=True,
)

cds_controls_card = dbc.Card(
    [
        cds_selection_params,
        cds_pause_detection_params,
    ],
    body=True,
)

pause_detection_tab_content = html.Div(
    [
        html.Hr(),
        html.Div([
            dbc.Container(
                [
                    html.H1("Ribosome Profiling Pause Detection Thesis"),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(controls_card, md=3),
                            dbc.Col(dbc.Spinner([dbc.Card([
                                dcc.Graph(id='comparative-figure'),
                                dcc.RangeSlider(min=0, max=80, step=1, value=[0,20], id='active-region-slider', marks=None,
                                tooltip={"placement": "bottom", "always_visible": False}),
                                dcc.Graph(id="graph"),    
                            ])], size='lg', type='border'), md=9, align='stretch'),
                        ],
                        align="center",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(md=4),
                            # dbc.Col(dbc.Card([dcc.Graph(id="graph")]), md=8),
                        ]
                    )
                ],
                fluid=True,
            )
        ]),
    ]
)

cds_analisis_tab_content = html.Div(
    [
        html.Hr(),
        html.Div([
            dbc.Container(
                [
                    html.H1("Ribosome Profiling Pause Detection Thesis"),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(cds_controls_card, md=3),
                            dbc.Col(dbc.Spinner([dbc.Card([
                                dcc.Graph(id='cds-analisys-graph'),
                                dcc.RangeSlider(min=0, max=80, step=1, value=[0,20], id='active-cds-region-slider', marks=None,
                                tooltip={"placement": "bottom", "always_visible": False}),
                                dcc.Graph(id="cds_boxplot"),    
                            ])], size='lg', type='border'), md=9, align='stretch'),
                        ],
                        align="center",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(md=4),
                            # dbc.Col(dbc.Card([dcc.Graph(id="graph")]), md=8),
                        ]
                    )
                ],
                fluid=True,
            )
        ]),
    ]
)

app.layout = html.Div(
    [  
        html.Div("RIBO - VEXA", style={'fontSize': 50, 'textAlign': 'center'}),
        dbc.Tabs(
            [
                dbc.Tab(
                    home_tab_content,
                    label="Home",
                ),
                dbc.Tab(
                    data_selection_tab_content,
                    label="Data Selection"
                ),
                dbc.Tab(
                    pause_detection_tab_content,
                    label="Pause Detection"
                ),
                dbc.Tab(
                    cds_analisis_tab_content,
                    label="CDS Analysis"
                ),
            ]
        ),
        html.Hr(),

        dash.page_container
    ]
)
@app.callback(
    Output('cds-analisys-graph', 'figure'),
    Input('cds-select', 'value'),
    Input('cds-frame', 'value'),
    Input('cds-pause-detection-method', 'value'),
    Input('cds-pause-rolling', 'value'),
    Input('cds-rolling-step', 'value'),
    Input('cds-rolling-window', 'value'),
    Input('cds-max-distance-between-outliers', 'value'),
    Input('cds-min-region-length', 'value'),
    Input('cds-percentile-slider', 'value'),
    Input('cds-z-score-slider', 'value'),
)
def analyze_cds(p_id, frame, pause_method, pause_rolling, rolling_step, rolling_window, max_distance, min_len, percentile_value, z_score_limit):
    counts = cdss.filter((pl.col('protein_ids') == p_id) & (pl.col('frame') == int(frame))).select('read_count').to_series().to_list()[0]
    position = np.arange(0, len(counts))

    if pause_rolling == 'Rolling':
        if pause_method == 'Boxplot':
            outliers = detect_rolling_pauses(counts, rolling_window, rolling_step, max_distance, min_len)
        elif pause_method == 'Percentile':
            outliers = detect_rolling_percentile_pauses(counts, rolling_window, rolling_step, max_distance, min_len, percentile_value/100)
        elif pause_method == 'Z-Score':
            outliers = detect_rolling_z_score_pauses(counts, rolling_window, rolling_step, max_distance, min_len, z_score_limit)
    else:
        if pause_method == 'Boxplot':
            outliers = detect_pauses(counts, max_distance, min_len)
        elif pause_method == 'Percentile':
            outliers = detect_percentile_pauses(counts, max_distance, min_len, percentile_value/100)
        elif pause_method == 'Z-Score':
            outliers = detect_z_score_pauses(counts, max_distance, min_len, z_score_limit)

    pauses = get_pause_to_graph(outliers)

    fig = go.Figure()
    fig.update_yaxes(fixedrange=True)
    fig.add_trace(
        go.Bar(x=position, y=counts, showlegend=False))   
    print(pauses)
    for pause in pauses:
        fig.add_vrect(
            x0=pause[0][0], x1=pause[0][1],
            fillcolor="LightSalmon", opacity=0.2,
            layer="above", line_width=0
        )
    fig.update_layout(
        margin=dict(l=20, r=20, t=30, b=0),
    )
    return fig


@app.callback(
    Output({'index': ALL}, 'active'),
    Output('selected_file', 'children'),
    Output('bam_histogram_graph', "figure"),
    Output('bam_size', "children"),
    Output('mapped_reads', "children"),
    Output('unmapped_reads', "children"),
    Output('read_count', "children"),
    Input({'index': ALL}, 'n_clicks'),
    Input({'index': ALL}, 'children')

)
def select_bam_from_list(items, items_names):
    actives = [False] * len(items)
    item = dash.callback_context.triggered[0]["prop_id"].split('.')[0]
    if item:
        item_dict = json.loads(item)
    else:
        return actives, "No file yet selected", go.Figure(), "0 bytes", "0 reads", '0 mapped reads', '0 unmapped reads'

    actives = actives[:item_dict['index']] + [True] + actives[item_dict['index']+1:]  

    file_name = items_names[item_dict['index']]
    
    if file_name not in bam_files_info:
        with Timer("Dictionary"):
            histogram, mapped, unmapped = get_read_length_histogram(DATA_PATH+file_name)

        histogram_as_list = process_bam_histogram(histogram)

        file_size = f'{os.path.getsize(DATA_PATH+file_name)  / 1024 / 1024:.2f} MB'
        bam_files_info[file_name] = {
            'size': file_size,
            'histogram': histogram_as_list,
            'reads': f"{mapped+unmapped:,} reads".replace(',','.'),
            'mapped_reads': f"{mapped:,} mapped reads".replace(',','.'),
            'unmapped_reads': f"{unmapped:,} unmapped reads".replace(',','.'),
        }

    fig = go.Figure([go.Bar(x=bam_files_info[file_name]['histogram'][0], y=bam_files_info[file_name]['histogram'][1])])

    return actives, file_name, fig, bam_files_info[file_name]['size'], bam_files_info[file_name]['reads'], bam_files_info[file_name]['mapped_reads'],bam_files_info[file_name]['unmapped_reads']

@app.callback(
    Output('active-region-slider', 'min'),
    Output('active-region-slider', 'max'),
    Output('active-region-slider', 'value'),
    Input('feature-select', 'value'),
)
def update_range_slider(selected_feature):
    # feature_sequences_rpf = df_rpf.filter(pl.col('transcript_id') == selected_feature)
    feature_sequences_rna = df_rna.filter(pl.col('transcript_id') == selected_feature)

    # cov_rpf = feature_sequences_rpf[selected_feature][1]  
    # cov_rna = feature_sequences_rna[selected_feature][1]

    # cov_rpf = feature_sequences_rpf.select('coverage_per_base').to_series().to_list()[0]
    cov_rna = feature_sequences_rna.select('coverage_per_base').to_series().to_list()[0]

    min_slider_window = 0
    max_slider_window = len(cov_rna)

    return min_slider_window, max_slider_window, [min_slider_window, max_slider_window]

@app.callback(
    Output('comparative-figure', 'figure'),
    Output('graph', 'figure'),
    Input('feature-select', 'value'),
    Input('mean-window-size', 'value'),
    Input('zero-handler-method', 'value'),
    Input('max-distance-between-outliers', 'value'),
    Input('min-region-length', 'value'),
    Input('RNA-offset','value'),
    Input('active-region-slider', 'value'),
    Input('pause-detection-method','value'),
    Input('percentile-slider', 'value'),
    Input('z-score-slider', 'value'),
    Input('pause-rolling', 'value'),
    Input('rolling-step', 'value'),
    Input('rolling-window', 'value'),
    Input('rna-coverage-region-filter', 'value'),
    Input('use_deseq2', 'value'),
    # Input('pause-method', 'value'),
    # Input('pause-window', 'value')
)
def update_figure(  selected_feature, window_size, method, max_distance, 
                    min_len, rna_offset, region, pause_method, percentile_value, 
                    z_score_limit, pause_rolling, rolling_step, rolling_window,
                    rna_coverage_minimum, use_deseq2):
    
    # df = filtered_gtf.get_transcripts_data([selected_feature]) 

    # feature_sequences_rpf = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rpf)
    # feature_sequences_rna = get_feature_sequence_and_coverage(df, dna.get_reference(), bam_rna)

    feature_sequences_rpf = df_rpf.filter(pl.col('transcript_id') == selected_feature)
    feature_sequences_rna = df_rna.filter(pl.col('transcript_id') == selected_feature)

    # cov_rpf = feature_sequences_rpf[selected_feature][1]  
    # cov_rna = feature_sequences_rna[selected_feature][1]
    if use_deseq2 == 'Yes':
        cov_rpf = feature_sequences_rpf.select('normalized_count').to_series().to_list()[0]
        cov_rna = feature_sequences_rna.select('normalized_count').to_series().to_list()[0]
    else:
        cov_rpf = feature_sequences_rpf.select('coverage_per_base').to_series().to_list()[0]
        cov_rna = feature_sequences_rna.select('coverage_per_base').to_series().to_list()[0]

    cov_rpf_over_rna = get_normalize_rpf_over_rna(cov_rpf, cov_rna, offset=rna_offset, window_size=window_size,method=method)

    if pause_rolling == 'Rolling':
        if pause_method == 'Boxplot':
            outliers = detect_rolling_pauses(cov_rpf_over_rna[region[0]:region[1]], rolling_window, rolling_step, max_distance, min_len)
        elif pause_method == 'Percentile':
            outliers = detect_rolling_percentile_pauses(cov_rpf_over_rna[region[0]:region[1]], rolling_window, rolling_step, max_distance, min_len, percentile_value/100)
        elif pause_method == 'Z-Score':
            outliers = detect_rolling_z_score_pauses(cov_rpf_over_rna[region[0]:region[1]], rolling_window, rolling_step, max_distance, min_len, z_score_limit)
    else:
        if pause_method == 'Boxplot':
            outliers = detect_pauses(cov_rpf_over_rna[region[0]:region[1]], max_distance, min_len)
        elif pause_method == 'Percentile':
            outliers = detect_percentile_pauses(cov_rpf_over_rna[region[0]:region[1]], max_distance, min_len, percentile_value/100)
        elif pause_method == 'Z-Score':
            outliers = detect_z_score_pauses(cov_rpf_over_rna[region[0]:region[1]], max_distance, min_len, z_score_limit)

    if rna_coverage_minimum != 0:
        outliers = filter_by_rna_coverage(outliers, cov_rna, rna_coverage_minimum)

    pauses = get_pause_to_graph(outliers, region[0])

    # print(pauses)

    position = np.arange(0, len(cov_rpf))
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,subplot_titles=(f"RPF for {selected_feature}", f"RNA for {selected_feature}", f"RPF/RNA for {selected_feature}"))
    fig.update_yaxes(fixedrange=True)
    fig.add_trace(
        go.Bar(x=position, y=cov_rpf,showlegend=False), row=1, col=1)   
    fig.add_trace(
        go.Bar(x=position, y=cov_rna,showlegend=False), row=2, col=1)
    fig.add_trace(
        go.Bar(x=position, y=cov_rpf_over_rna,showlegend=False), row=3, col=1)
    for pause in pauses:
        fig.add_vrect(
            x0=pause[0][0], x1=pause[0][1],
            fillcolor="LightSalmon", opacity=0.2,
            layer="above", line_width=0, row=3, col=1
        )
    



# region Show IGV

# """No se ve como se debe. El IGV considera todos las partes de gen, mientras que el transctipro muestra el arn msj (es decir, los exones)"""
    
    # gene_traces = gtf.format_gene_igv(transcript_to_genes[selected_feature])[transcript_to_genes[selected_feature]]
    # print(gene_traces)
    # offset = 1
    
    # x_axis_offset = gene_traces[selected_feature]['exon'][0][0]
    # for e in gene_traces[selected_feature]['exon']:
    #     if e[0] < x_axis_offset:
    #         x_axis_offset = e[0]
    # for transcript in gene_traces.keys():
    #     fig.add_trace(
    #         go.Scatter(
    #             x = np.linspace(gene_traces[transcript]['exon'][0][0]-x_axis_offset, gene_traces[transcript]['exon'][-1][1]-x_axis_offset, 20), 
    #             y = [offset] * 20, 
    #             mode='lines+markers', 
    #             marker= dict(
    #                 symbol="arrow",
    #                 size=10,
    #                 angleref="previous",
    #             ),
    #             line=dict(
    #                 color='royalblue', 
    #                 width=1
    #             ),
    #             showlegend=False
    #         ), row=4, col=1
    #     )
    #     for exon in gene_traces[transcript]['exon']:
    #         fig.add_trace(go.Scatter(x = [x - x_axis_offset for x in exon], y = [offset,offset], line=dict(color='royalblue', width=4),showlegend=False), row=4, col=1)
        
    #     for cds in gene_traces[transcript]['CDS']:
    #         fig.add_trace(go.Scatter(x = [x - x_axis_offset for x in cds], y = [offset,offset], line=dict(color='royalblue', width=8),showlegend=False), row=4, col=1)
        
    #     for sc in gene_traces[transcript]['start_codon']:
    #         fig.add_trace(go.Scatter(x = [x - x_axis_offset for x in sc], y = [offset,offset], line=dict(color='royalblue', width=20),showlegend=False), row=4, col=1)

    #     for sc in gene_traces[transcript]['stop_codon']:
    #         fig.add_trace(go.Scatter(x = [x - x_axis_offset for x in sc], y = [offset,offset], line=dict(color='green', width=20),showlegend=False), row=4, col=1)

    #     offset += 0.5
# endregion


    fig2 = make_subplots(rows=1, cols=3, shared_xaxes=True,subplot_titles=(f"RPF boxplot" f"RPF/RNA boxplot"))
    fig2.add_trace(go.Box(y=cov_rpf, quartilemethod="linear"), row=1, col=1)
    fig2.add_trace(go.Box(y=cov_rna, quartilemethod="linear"), row=1, col=2)
    fig2.add_trace(go.Box(y=cov_rpf_over_rna, quartilemethod="linear"), row=1, col=3)

    fig.update_layout(
        margin=dict(l=20, r=20, t=30, b=0),
    )

    return fig,fig2

@app.callback(
    Output('feature-select', 'options'),
    Input('coverage-filter-slider', 'value')
)
def update_dropdown(selected_filter):
    transcript_ids_filtered = (df_rpf
                .lazy()
                .filter(pl.col('coverage_percentage') > selected_filter/100)
                .select('transcript_id')
            ).collect().to_series().to_list()

    return transcript_ids_filtered


@app.callback(
    Output('mean-window-size', 'disabled'),
    Input('zero-handler-method', 'value')
)
def enable_1(zero_handler):
    return zero_handler != 'Mean'


@app.callback(
    Output('rolling-step', 'disabled'),
    Output('rolling-window', 'disabled'),
    Input('pause-rolling', 'value')
)
def enable_2(pause_rolling_detection):
    is_fixed = pause_rolling_detection != 'Rolling' 
    return is_fixed, is_fixed

@app.callback(
    Output('percentile-slider', 'disabled'),
    Output('z-score-slider', 'disabled'),
    Input('pause-detection-method', 'value')
)
def enable_3(pause_detection_method):
    return pause_detection_method != "Percentile" ,pause_detection_method != "Z-Score"

@app.callback(
    Output('cds-rolling-step', 'disabled'),
    Output('cds-rolling-window', 'disabled'),
    Input('cds-pause-rolling', 'value')
)
def enable_4(pause_rolling_detection):
    is_fixed = pause_rolling_detection != 'Rolling' 
    return is_fixed, is_fixed

@app.callback(
    Output('cds-percentile-slider', 'disabled'),
    Output('cds-z-score-slider', 'disabled'),
    Input('cds-pause-detection-method', 'value')
)
def enable_5(pause_detection_method):
    return pause_detection_method != "Percentile" ,pause_detection_method != "Z-Score"



if __name__ == "__main__":
    app.run(debug=True)