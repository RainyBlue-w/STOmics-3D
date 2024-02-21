## In[] env

from dash import Dash
from dash import dcc, html, dash_table, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash import ALL, MATCH, ALLSMALLER
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components as fac
import feffery_utils_components as fuc
import dash_bootstrap_components as dbc
from dash_extensions.enrich import Output, Input, html, callback, DashProxy
from dash_extensions.enrich import MultiplexerTransform, Trigger, TriggerTransform, ServersideOutputTransform, Serverside, LogTransform

import plotly.express as px
import plotly.graph_objects as go
from plotnine import *

import flask
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import math
from functools import reduce

from typing import List, Dict
import diskcache

background_callback_manager = DiskcacheManager(diskcache.Cache("cache/"))

## In[] data

exp_data = { # adata
  'E7.5': sc.read_h5ad("data/E7.5_HC0.5_min400.h5ad"),
  'E7.75': sc.read_h5ad("data/E7.75_HC0.5_min400.h5ad"),
  'E8.0': sc.read_h5ad("data/E8.0_HC0.5_min400.h5ad")
}

ctp_cmap = pd.read_csv("/rad/wuc/dash_data/spatial/celltype_cmap.csv")
ctp_cmap = dict(zip(ctp_cmap['celltype'], ctp_cmap['Epiblast']))


def create_3dviewer(requests_pathname_prefix: str = None) -> Dash:

  # In[] functions:
  def show_expViolin(adata, feature, **kws):
    data = adata[:,feature].to_df()[feature]
    # data = data[data>0]
    fig = go.Figure(
      data = go.Violin(
        x=data, y0=f'{feature}({len(data)})', line_color='black',
        fillcolor='lightseagreen', opacity=0.6,
        orientation='h', side='positive', width=1.5, **kws,
      )
    )
    
    fig.update_layout(
      plot_bgcolor = 'rgba(200,200,200,0.1)', showlegend=False
    ).update_yaxes(
      gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
    ).update_xaxes(
      dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
    )
    return fig

  def show_ctpExpViolin(adata, feature, **kws):
    pdf = pd.concat([adata[:,feature].to_df(), adata.obs.celltype], axis=1)
    counts = pdf.celltype.value_counts()
    counts = counts[counts>0]
    sorted_ctp = counts.index.to_list()
    pdf['celltype'] = pd.Categorical(pdf['celltype'].to_list(),
                                    categories=sorted_ctp[::-1])
    fig = px.violin(
      pdf, x=feature, y='celltype', color = 'celltype', 
      color_discrete_map=ctp_cmap, orientation='h', height=800,
    ).update_traces(
      side='positive', width=1.5, **kws,
    ).update_layout(
      plot_bgcolor = 'rgba(200,200,200,0.1)',
    ).update_yaxes(
      gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
    ).update_xaxes(
      dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
    )
    return fig

  def show_multiFeatures_expViolin(adata, features_dict, **kws):
    
    fig = go.Figure()
    
    filt_dict = {}
    for color,feature  in features_dict.items():
        if feature:
            filt_dict[color] = feature

    for color in list(filt_dict.keys())[::-1]:
      data = adata[:,filt_dict[color]].to_df()[filt_dict[color]]
      fig.add_trace(
        go.Violin(
          x=data, y0=f'{filt_dict[color]}({len(data)})', box_visible=False, 
          line_color='black', meanline_visible=False,
          fillcolor=color, opacity=0.6,
          orientation='h', side='positive', width=1.5, **kws, 
        )
      )
    fig.update_layout(
      plot_bgcolor = 'rgba(200,200,200,0.1)', showlegend=False
    ).update_yaxes(
      gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
    ).update_xaxes(
      dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
    )
    
    return fig

  def show_multiFeatures_ctpExpViolin(adata, features_dict, **kws):
    
    from plotly.subplots import make_subplots
    
    filt_dict = {}
    for color,feature  in features_dict.items():
        if feature:
            filt_dict[color] = feature
    features = list(filt_dict.values())

    pdf = pd.concat([adata[:,features].to_df(), adata.obs.celltype], axis=1)
    pdf = pdf.melt(id_vars='celltype')
    pdf = pdf.rename(columns = {'variable': 'Gene', 'value': 'expression'})
    # pdf = pdf[pdf['expression']>0]

    pdf.celltype = pd.Categorical(pdf.celltype, ordered=True)
    # counts = pdf.groupby('Gene').apply(lambda x: x.value_counts())

    fig = px.violin(
      pdf, x='expression', y='celltype', color = 'celltype', 
      color_discrete_map=ctp_cmap, orientation='h', height=800,
      animation_frame='Gene', 
    ).update_traces(
      side='positive', width=1.5, **kws,
    ).update_layout(
      plot_bgcolor = 'rgba(200,200,200,0.1)',
    ).update_yaxes(
      gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
    ).update_xaxes(
      dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
    )
      
    return fig

  def vector_to_rgba(v):
    color = list(v.keys())
    color = [str(math.ceil(v[i])) if i in color else '244' for i in ['R', 'G', 'B'] ]
    if(all([ i=='244' for i in color])):
      rgba = 'rgba(244,244,244,1)'
    else:
      rgba = f'rgba({color[0]}, {color[1]}, {color[2]}, 1)'
      
    return rgba

  def multiGenes_show_color(adata, genes_dict):
    import numpy
    tmp = {}
    for key,value  in genes_dict.items():
        if value:
            tmp[key] = value
    genes_dict = tmp
    colors = list(genes_dict.keys())
    others = [i for i in ['R', 'G', 'B'] if i not in colors]
    genes = list(genes_dict.values())

    exp = adata[:, genes].to_df()
    exp.columns = colors

    delta = exp.div(exp.max(axis=0), axis=1)*244
    delta[others] = 0

    def delta_geoMean(a,b):
      geoMean = numpy.sqrt((a**2+b**2)/2)
      # geoMean = ((a**3+b**3)/2)**(1/3)
      return geoMean
    def mean(a,b, c=None):
      if c:
        return (a+b+c)/3
      else:
        return (a+b)/2

    if len(colors)==1:
      color = pd.DataFrame({
          colors[0] : 244,
          others[0] : 244-delta[colors[0]],
          others[1] : 244-delta[colors[0]],
      })
    elif len(colors)==2:
      color = pd.DataFrame({
          colors[0] : 244-delta[colors[1]],
          colors[1] : 244-delta[colors[0]],
          others[0] : 244-delta_geoMean(delta[colors[1]],delta[colors[0]]),
      })
    elif len(colors)==3:
      color = pd.DataFrame({
          'R' : 244-delta_geoMean(delta['G'], delta['B']),
          'G' : 244-delta_geoMean(delta['R'], delta['B']),
          'B' : 244-delta_geoMean(delta['R'], delta['G']),
      })
    
    color['RGBA'] = color.apply(vector_to_rgba, axis=1)
    return color['RGBA']

  def hex_to_rgbList(hex_color):
    hex_color = hex_color.replace(' ', '').replace('#', '')
    if len(hex_color) == 6:
      r = int(hex_color[0:2], 16)
      g = int(hex_color[2:4], 16)
      b = int(hex_color[4:6], 16)
    return [r,g,b]

  def mix_multipy(color, alpha):

    def multipy(x,y):
      return x*y/255

    def mix(x, y):
      alpha = x[3]+y[3]-x[3]*y[3]
      if alpha==0:
        return [244,244,244, 0]
      else:
        R = np.round( (x[3]*(1-y[3])*x[0]+x[3]*y[3]*multipy(x[0],y[0])+(1-x[3])*y[3]*y[0])/alpha).astype(int)
        G = np.round( (x[3]*(1-y[3])*x[1]+x[3]*y[3]*multipy(x[1],y[1])+(1-x[3])*y[3]*y[1])/alpha).astype(int) 
        B = np.round( (x[3]*(1-y[3])*x[2]+x[3]*y[3]*multipy(x[2],y[2])+(1-x[3])*y[3]*y[2])/alpha).astype(int)
        return [R,G,B,alpha]

    array = []
    for c,a in zip(color, alpha):
      array.append(c.copy())
      array[-1].append(a)

    res = reduce(mix, array)
    res = f'rgb{res[0],res[1],res[2]}'

    return res

  def color_mixer(adata, genes_dict):
    genes_dict_copy = genes_dict.copy()
    _ = [genes_dict_copy.pop(color) for color in genes_dict.keys() if not genes_dict[color]]
    colors = [hex_to_rgbList(c) for c in genes_dict_copy.keys()]
    genes = list(genes_dict_copy.values())

    exp = adata[:,genes].to_df()
    alpha = exp.div(exp.max(axis=0), axis=1)
    cell_colors = alpha.apply(axis=1, func=lambda row: mix_multipy(colors,row) )

    return cell_colors

  def cal_moran_3D(adata):
    tmp = adata.copy()
    sq.gr.spatial_neighbors(tmp, spatial_key='X_spatial')
    sq.gr.spatial_autocorr(tmp, mode='moran', n_jobs=1)
    df = tmp.uns['moranI'][['I']]
    df.columns = ["Moran's I"]
    return df

  # In[] app

  app = DashProxy(
    __name__, 
    external_stylesheets=[
      dbc.themes.BOOTSTRAP
    ],
    external_scripts = [
      {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
    ],
    transforms=[
      LogTransform(), ServersideOutputTransform(), MultiplexerTransform(), TriggerTransform()
    ],
    requests_pathname_prefix = requests_pathname_prefix
  )

  header = dbc.NavbarSimple(
      [
          dbc.DropdownMenu(
              children=[
                  # dbc.DropdownMenuItem('spatial', href='/'),
                  # dbc.DropdownMenuItem("atlas", href='/atlas'),
                  # dbc.DropdownMenuItem("reik", href='/reik'),
              ],
              nav=True,
              in_navbar=True,
              label="Dataset",
          ),
      ],
      brand="Omics-viewer",
      color="dark",
      dark=True,
      # sticky='top',
      style = {"height": "6vh"}
  )

  # In[] global vars

  colorPicker_swatches = [
    "#25262b", "#868e96", "#fa5252", "#e64980", "#be4bdb", "#7950f2", "#4c6ef5",
    '#225ea8', "#228be6", "#15aabf", "#12b886", "#40c057", "#82c91e", "#fab005", "#fd7e14",
  ]

  initColor_multiName = [
    "#fa5252", "#228be6", "#40c057", "#fd7e14", "#be4bdb", "#e64980", "#15aabf", "#fab005", "#868e96", 
  ]

  config_scatter3d = {
    'toImageButtonOptions': {
      'format': 'png', # one of png, svg, jpeg, webp,
      'scale': 3
    }
  }

  config_violin = {
    'toImageButtonOptions': {
      'format': 'png', # one of png, svg, jpeg, webp,
      'scale': 3
    }
  }


  # In[] widgets

  SET_STORE_JSONtoPlot_3D = html.Div(
    [
      dcc.Store(data={}, id='STORE_obs_3D'),
      dcc.Store(data={}, id='STORE_cellsObsFilter_3D'), # Server-side
      dcc.Store(data={}, id='STORE_cellsSliceFilter_3D'), # Server-side
      dcc.Store(data={}, id='STORE_cellsExpFilter_3D'), # Server-side
      dcc.Store(data={}, id='STORE_singleExp_3D'),
      dcc.Store(data={}, id='STORE_multiExp_3D'),
      dcc.Store(data={}, id='STORE_mixedColor_3D'),
      dcc.Store(data=False, id='STORE_ifmulti_3D'),
      dcc.Store(data=ctp_cmap, id='STORE_ctpCmap_3D', storage_type='local'),
      dcc.Store(id='STORE_cellsCtpFilter_3D'), # Server-side
      dcc.Store(id='STORE_cellsIntersection_3D'),
      dcc.Store(id='test'),
    ]
  )

  init_range = dict(
    x_min = np.floor(exp_data['E7.5'].obs.x.min()/10)*10, x_max = np.ceil(exp_data['E7.5'].obs.x.max()/10)*10,
    y_min = np.floor(exp_data['E7.5'].obs.y.min()/10)*10, y_max = np.ceil(exp_data['E7.5'].obs.y.max()/10)*10,
    z_min = np.floor(exp_data['E7.5'].obs.z.min()/10)*10, z_max = np.ceil(exp_data['E7.5'].obs.z.max()/10)*10,
  )

  SET_STORE_Ranges_3D = html.Div(
    [
      dcc.Store(
        data = init_range, 
        id='STORE_previewRange_3D'
      ),
      dcc.Store(id='STORE_sliceRange_3D'),
      dcc.Store(
        data = init_range,
        id='STORE_maxRange_3D'),
    ]
  )

  def iconHover_colorPicker(init_color: str, id: Dict, swatches: List[str], placement='left', trigger='click'):
    return fac.AntdPopover(
      # openDelay=200,
      placement = placement,
      trigger= trigger,
      children=[
        dmc.ActionIcon(
          DashIconify(icon = 'fluent:circle-48-filled', color=init_color, width=48),
          variant='transparent', id=id['dmc_ActionIcon'], mt=3
        ),
      ],
      content = [
        dmc.ColorPicker(id=id['dmc_ColorPicker'], format='hex', value=init_color, swatches=swatches),
        dmc.TextInput(value=init_color, id=id['dmc_TextInput']),
      ]
    )

  def drawerContent_ctpColorPicker(celltypes: List[str], cmap: Dict, swatches=colorPicker_swatches):
    stack = dmc.Stack(
      children=[
        dmc.Grid(
          gutter = 2,
          children=[
            dmc.Col(dmc.Text(ctp), span=10),
            dmc.Col(
              iconHover_colorPicker(
                id = {
                  'dmc_ActionIcon': {'type': 'ACTIONICON_colorCtp_3D', 'id': ctp},
                  'dmc_ColorPicker': {'type': 'COLORPICKER_colorCtp_3D', 'id': ctp},
                  'dmc_TextInput': {'type': 'TEXT_colorCtp_3D', 'id': ctp},
                },  placement='right', init_color=cmap[ctp], swatches = [ctp_cmap[ctp]]+swatches
              ),
              span=2
            )
          ],
        )
        for ctp in celltypes
      ],
    )

    return stack

  def selectData_newFilter(index: int=0):
    component = html.Div([
      dcc.Store(id={'type': 'STORE_filterBodyPreservedCells_3D', 'index': index}), # store preserved cells
      dmc.Grid([
        dmc.Col(
          [
            dmc.Text('Item', className='dmc-Text-label'),
            dcc.Dropdown(id={'type': 'DROPDOWN_filterItem_3D', 'index': index}, clearable=True, searchable=True,
                        persistence=True, persistence_type='local'),
          ],
          span=7),
        dmc.Col(
          [
            dmc.Text('Type', className='dmc-Text-label'),
            dcc.Dropdown(['numeric', 'categorical'], id={'type': 'DROPDOWN_filterType_3D', 'index': index}, 
                        clearable=False, searchable=False, persistence=True, persistence_type='local'),
          ],
          span=5
        ),
        dmc.Col(dmc.Switch(onLabel='ON', offLabel='OFF', id={'type': 'SWITCH_filterApply_3D', 'index': index}, size='lg', checked=False,), span=4),
        dmc.Col(
          fac.AntdPopover(
            placement='right',
            trigger='click',
            children=[
              dmc.Button('Configure', color='dark', size='sm', id={'type': 'BUTTON_filterConfigure_3D', 'index': index}, fullWidth=True),
            ],
            id={'type': 'POPOVER_filterBody_3D', 'index': index},
            style={'width': '30vh'}
          ),
          span=8
        ),
        dmc.Col(
          dmc.Text('Selected cells: ', id={'type': 'TEXT_filterBodyPreservedCells_3D', 'index': index}, className='dmc-Text-sidebar-tips'),
          span=12
        )
      ]),
      dmc.Space(h=5),
      dmc.Divider(variant='dashed'),
      dmc.Space(h=5),
    ])
    return component

  # In[] tab

  spatial_tab_plotFeature3D = dbc.Tab(
    [dmc.Grid([
      # options
      dmc.Col([
        fuc.FefferySticky([
          fac.AntdSpace(
            size=0,
            direction='vertical',
            className='fac-AntdSpace-sideBar',
            children=[
              # Select data
              fac.AntdCollapse(
                isOpen = True,
                forceRender = True,
                className = 'fac-AntdCollapse-sidebar',
                ghost=True,
                title = dmc.Text('Select data', className='dmc-Text-sidebar-title'),
                children = [
                  dmc.Grid([
                    dmc.Col(dmc.Text("Stage:", className='dmc-Text-label'), span=4),
                    dmc.Col([
                      dcc.Dropdown(
                        ['E7.5', 'E7.75', 'E8.0'],
                        'E7.5',
                        id="DROPDOWN_stage_3D",
                        clearable=False,
                        searchable=True,
                      ),
                    ], span=8),
                    dmc.Col(dmc.Text(id='TEXT_dataSummary_3D', color='gray'), span=12),
                    dmc.Col(
                      fac.AntdCollapse(
                        title = 'Filter(metadata)', className='fac-AntdCollapse-inline',
                        forceRender=True, isOpen=False, ghost=True,
                        children = [
                          html.Div(
                            children= [selectData_newFilter(index=i) for i in [0]],
                            id='CONTAINER_filterList_3D',
                          ),
                          dcc.Store(id='STORE_filterCurNumber_3D', data=1),
                          dmc.Grid(
                            [
                              dmc.Col(dmc.Button(
                                id='BUTTON_addFilter_3D', color='teal', fullWidth=True,
                                children = DashIconify(icon="fluent:add-square-20-regular", width=20)
                              ), span=6),
                              dmc.Col(dmc.Button(
                                id='BUTTON_deleteFilter_3D', color='red', fullWidth=True,
                                children = DashIconify(icon="fluent:subtract-square-20-regular", width=20)
                              ), span=6),
                            ],
                          ),
                        ]
                      ),
                      span=12
                    )
                  ], gutter='xs'),
                ]
              ),
              # Plot options
              fac.AntdCollapse(
                isOpen = True,
                forceRender = True,
                className = 'fac-AntdCollapse-sidebar',
                ghost=True,
                title = dmc.Text('Plot options', className='dmc-Text-sidebar-title'),
                children = [          
                  dmc.Tabs(
                    [
                      dmc.TabsList([
                        dmc.Tab('Settings', value='settings'),
                        dmc.Tab('Single', value='single'),
                        dmc.Tab('Multiple', value='multi'),
                      ], grow=True),
                      # settings
                      dmc.TabsPanel(
                        [
                          dmc.Tabs(
                            [
                              dmc.TabsList([
                                dmc.Tab('Scatter-3D', value='Scatter-3D'),
                                dmc.Tab('Violin', value='Violin')
                              ], grow=False),
                              # scatter-3d
                              dmc.TabsPanel(
                                [
                                  
                                  dmc.Divider(label = 'Points', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.Grid([
                                    dmc.Col(dmc.Text('Point size:', className='dmc-Text-label'), span=5),
                                    dmc.Col(dmc.NumberInput(
                                      value=3, step=0.5, min=0.1, id='NUMBERINPUT_scatter3dPointsize_3D', precision=1,
                                      persistence = True, persistence_type = 'local'
                                    ), span=7),
                                  ], justify='center', gutter=3, className='dmc-Grid-center'),
                                  dmc.Space(h=5),

                                  dmc.Switch(label='Hide non-expressing cells', id='SWITCH_hideZero_3D',  size='md',
                                            onLabel=DashIconify(icon='solar:eye-closed-linear', width=14), 
                                            offLabel=DashIconify(icon='solar:eye-linear', width=14),
                                            persistence = False, persistence_type = 'local'),
                                  dmc.Space(h=5),
                                  
                                  dmc.Divider(label = 'Axes', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.Text('Projection type:', className='dmc-Text-label'),
                                  dmc.SegmentedControl(
                                    value='orthographic', 
                                    data=[
                                      {'value': 'perspective', 'label': 'Perspective'},
                                      {'value': 'orthographic', 'label': 'Orthographic'},
                                    ], 
                                    fullWidth=True, id='SEGMENTEDCONTROL_projection_3D',
                                    persistence = True, persistence_type = 'local',
                                  ),
                                  dmc.Space(h=5),
                                  dmc.Switch(label='Hide axes', id='SWITCH_hideAxes_3D', size='md',
                                    onLabel=DashIconify(icon='solar:eye-closed-linear', width=14), 
                                    offLabel=DashIconify(icon='solar:eye-linear', width=14),
                                    persistence = True, persistence_type = 'local'),
                                  
                                  dmc.Divider(label='Download', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.Text('tip: replot to take effect', className='dmc-Text-sidebar-tips'),
                                  dmc.Grid([
                                    dmc.Col(dmc.Select( label = 'type', id='NUMBERINPUT_scatter3dFigtype_3D',
                                      value='png', data = ['svg', 'png', 'jpeg', 'webp'],
                                      persistence = True, persistence_type = 'local', 
                                    ), span=6),
                                    dmc.Col(dmc.NumberInput( label = 'scale(resolution)', id='NUMBERINPUT_scatter3dFigscale_3D',
                                      value=3, step=1, min=1, icon=DashIconify(icon='uim:multiply', width=16),
                                      persistence = True, persistence_type = 'local', 
                                    ), span=6),
                                  ], justify='center', gutter=3, className='dmc-Grid-center'),
                                ],
                                value = 'Scatter-3D'
                              ),
                              # violin
                              dmc.TabsPanel(
                                [
                                  dmc.Divider(label = 'Points', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.SegmentedControl(
                                    value='outliers',
                                    data = [
                                      {'value': 'none', 'label': 'none'},
                                      {'value': 'outliers', 'label': 'outliers'},
                                      {'value': 'all', 'label': 'all'}
                                    ],
                                    fullWidth=True, id='SEGMENTEDCONTROL_violinPoints_3D',
                                    persistence = True, persistence_type = 'local',
                                  ),
                                  dmc.Grid(
                                    [
                                      dmc.Col(dmc.NumberInput(label='position', value=0, step=0.1, min=-2, max=2, 
                                                      id='NUMBERINPUT_violinPointpos_3D', precision=2,
                                                      persistence = True, persistence_type = 'local',), span=4),
                                      dmc.Col(dmc.NumberInput(label='size', value=2.5, step=0.5, min=0, max=10,
                                                      id='NUMBERINPUT_violinPointsize_3D', precision=1,
                                                      persistence = True, persistence_type = 'local',), span=4),
                                      dmc.Col(dmc.NumberInput(label='jitter', value=0.15, step=0.05, min=0, max=1,
                                                      id='NUMBERINPUT_violinPointjitter_3D', precision=2,
                                                      persistence = True, persistence_type = 'local',), span=4),
                                    ],
                                  ),
                                  dmc.Divider(label = 'Box', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.SegmentedControl(
                                    value='all',
                                    data = [
                                      {'value': 'none', 'label': 'none'},
                                      {'value': 'box', 'label': 'box'},
                                      {'value': 'meanline', 'label': 'mean'},
                                      {'value': 'all', 'label': 'all'}
                                    ],
                                    id='SEGMENTEDCONTROL_violinBox_3D', fullWidth=True,
                                    persistence = True, persistence_type = 'local',
                                  ),
                                  dmc.NumberInput(label='Box width', value=0.5, step=0.1, min=0, max=1,
                                                  id='NUMBERINPUT_violinBoxwidth_3D', precision=1,
                                                  persistence = True, persistence_type = 'local',),
                                  
                                  dmc.Divider(label='Download', labelPosition='center', variant='dashed', className='dmc-divider-sidebar-inline'),
                                  dmc.Text('tip: replot to take effect', className='dmc-Text-sidebar-tips'),
                                  dmc.Grid([
                                    dmc.Col(dmc.Select( label = 'type', id='NUMBERINPUT_violinFigtype_3D',
                                      value='png', data = ['svg', 'png', 'jpeg', 'webp'],
                                      persistence = True, persistence_type = 'local', 
                                    ), span=6),
                                    dmc.Col(dmc.NumberInput( label = 'scale(resolution)', id='NUMBERINPUT_violinFigscale_3D',
                                      value=3, step=1, min=1, icon=DashIconify(icon='uim:multiply', width=16),
                                      persistence = True, persistence_type = 'local', 
                                    ), span=6),
                                  ], justify='center', gutter=3, className='dmc-Grid-center'),
                                ],
                                value = 'Violin'
                              ),
                            ],
                            value = 'Scatter-3D',
                            variant = 'pills',
                            color = 'grape'
                          ),
                        ],
                        value = 'settings',
                      ),
                      # single
                      dmc.TabsPanel(
                        [
                          dmc.Grid([
                            dmc.Col([
                              dcc.Dropdown(
                                options = exp_data['E7.5'].var_names,
                                value = 'Cdx1',
                                id="DROPDOWN_singleName_3D",
                                clearable=False
                              ),
                            ], span=10),
                            dmc.Col([
                              iconHover_colorPicker(
                                id={
                                  'dmc_ActionIcon': 'ACTIONICON_colorSingle_3D', 
                                  'dmc_ColorPicker': 'COLORPICKER_single_3D', 
                                  'dmc_TextInput': 'TEXT_colorSingle_3D',
                                  },
                                init_color='#225ea8', swatches=colorPicker_swatches,
                              )
                            ], span=2),
                            dmc.Col([
                              dmc.Button('Plot', id='BUTTON_singlePlot_3D', color='dark', fullWidth=True,
                                        leftIcon=DashIconify(icon="gis:cube-3d", width=24)),
                            ], span=12),
                          ], gutter='xs')  
                        ],
                        value='single',
                      ),
                      # multi
                      dmc.TabsPanel(
                        [
                          # extendable selector
                          html.Div(
                            [
                              dmc.Grid([
                                dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': 0}), span=10),
                                dmc.Col(
                                  iconHover_colorPicker(
                                    id={
                                      'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': 0}, 
                                      'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': 0}, 
                                      'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': 0},
                                      },
                                    init_color=initColor_multiName[0], swatches=colorPicker_swatches,
                                  ),
                                  span=2
                                ),
                              ]),
                              dmc.Grid([
                                dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': 1}), span=10),
                                dmc.Col(
                                  iconHover_colorPicker(
                                    id={
                                      'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': 1}, 
                                      'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': 1}, 
                                      'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': 1},
                                      },
                                    init_color=initColor_multiName[1], swatches=colorPicker_swatches,
                                  ),
                                  span=2
                                ),
                              ]),
                            ],
                            id = 'DIV_multiNameDynamic_3D',
                          ),
                          dcc.Store(data=2, id='STORE_multiNameCurNumber'),
                          # buttons
                          dmc.Grid(
                            [
                              dmc.Col(dmc.Button(
                                'Add', id='BUTTON_addFeature_3D', color='teal', fullWidth=True,
                                leftIcon=DashIconify(icon="fluent:add-square-20-regular", width=20)
                              ), span=23),
                              dmc.Col(dmc.Button(
                                'Delete', id='BUTTON_deleteFeature_3D', color='red', fullWidth=True,
                                leftIcon=DashIconify(icon="fluent:subtract-square-20-regular", width=20)
                              ), span=27),
                              dmc.Col(dmc.Button(
                                'Plot', id='BUTTON_multiPlot_3D', color='dark', fullWidth=True,
                                leftIcon=DashIconify(icon="gis:cube-3d", width=24),
                              ), span=50),
                            ],
                            columns=50,
                          ),
                          dcc.Store(id='STORE_multiNameInfo_3D'),
                        ],
                        value='multi',
                      ),
                    ], 
                    orientation = 'horizontal',
                    className = 'dmc-Tabs-inline',
                    # variant = 'pills',
                    value = 'single',
                  ),
                ]
              ),
              # Slicer
              fac.AntdCollapse(
                isOpen = False,
                forceRender = True,
                className = 'fac-AntdCollapse-sidebar',
                ghost=True,
                title = dmc.Text('Slicer', className='dmc-Text-sidebar-title'),
                children = [
                  dmc.Grid([
                    dmc.Col(dmc.Text('x', className='.dmc-Text-label-center'), span=2),
                    dmc.Col(
                      dcc.RangeSlider(
                        step=10, id='SLIDER_Xrange_3D',
                        marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                      ),
                      span=10
                    ),
                    dmc.Col(dmc.Text('y', className='.dmc-Text-label-center'), span=2),
                    dmc.Col(
                      dcc.RangeSlider(
                        step=10, id='SLIDER_Yrange_3D',
                        marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                      ),
                      span=10
                    ),
                    dmc.Col(dmc.Text('z', className='.dmc-Text-label-center'), span=2),
                    dmc.Col(
                      dcc.RangeSlider(
                        step=10, id='SLIDER_Zrange_3D',
                        marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                      ),
                      span=10
                    ),
                    dmc.Col(
                      dmc.Switch(size='md', id='SWITCH_previewBox_3D', label='Preview', checked=False),
                      span=12
                    ),
                    dmc.Col(
                      dmc.Button(
                        'Slice', color='red', id='BUTTON_slice_3D', fullWidth=True, 
                        leftIcon=DashIconify(icon='fluent:screen-cut-20-regular', width=20),
                      ),
                      span=6
                    ),
                    dmc.Col(
                      dmc.Button(
                        'Recover', color='teal', id='BUTTON_recover_3D', fullWidth=True, 
                        leftIcon=DashIconify(icon='fluent:arrow-sync-circle-20-regular', width=20),
                      ),
                      span=6
                    )
                  ]),
                  SET_STORE_Ranges_3D,
                ],
              ),
              # Moran
              fac.AntdCollapse(
                isOpen = False,
                forceRender = True,
                className = 'fac-AntdCollapse-sidebar',
                ghost=True,
                title = dmc.Text('Compute SVG(moran)', className='dmc-Text-sidebar-title'),
                children = [
                  dmc.Grid([
                    dmc.Col(
                      dmc.Button('Compute', id='BUTTON_calMoran_3D', color='dark', fullWidth=True,
                          leftIcon = DashIconify(icon='fluent:clipboard-math-formula-20-regular', width=20) ),
                      span=6,
                    ),
                    dmc.Col(
                      dmc.Button('Result', id='BUTTON_showMoran_3D', fullWidth=True,
                          leftIcon = DashIconify(icon='fluent:clipboard-checkmark-20-regular', width=20) ),
                      span=6,
                    ),
                    dmc.Text('Using current cells to compute SVGs', className='dmc-Text-sidebar-tips'),
                  ]),
                  dbc.Offcanvas(
                    [dash_table.DataTable(
                      id='DATATABLE_moranRes_3D',
                      sort_action="native", page_action='native', filter_action="native",
                      page_current= 0, page_size= 20, fill_width=True,
                      style_cell={'textAlign': 'center'},
                      style_table={'overflowX': 'auto'},
                    )],
                    title = 'SVGs:',
                    placement='end', scrollable=True, backdrop=False, is_open=False,
                    id = 'OFFCANVAS_moranRes_3D',
                  ),
                ],
              ),
            ],
          ),
        ], top=10),
      ], span=9),
      # viewer
      dmc.Col([
        SET_STORE_JSONtoPlot_3D,
        # scatter3d
        dmc.Grid([
          dmc.Col([
            dcc.Graph(figure={}, id="FIGURE_3Dexpression", 
                      className='dcc-Graph-scatter3d', config=config_scatter3d),
          ], span=20),
          dmc.Col([
            dcc.Graph(figure={}, id="FIGURE_3Dcelltype",
                      className='dcc-Graph-scatter3d', config=config_scatter3d),
          ], span=20),
          dmc.Col([
            # DIY-legend
            dmc.Grid([
              # set colors
              dmc.Col(dmc.Button(
                'Setting Colors', variant="gradient", gradient={"from": "grape", "to": "pink", "deg": 35},
                id='BUTTON_setColors_3D', fullWidth=True,
                leftIcon=DashIconify(icon='fluent:color-20-regular', width=20)
              ), span=12),
              # invert selection
              dmc.Col(
                dmc.Button(
                  DashIconify(icon='system-uicons:reverse', width=21), variant='light', color='gray',
                  id='BUTTON_invertSelectionCtp_3D', fullWidth=True,),
                span=4),
              # clear selection
              dmc.Col(dmc.Button(
                DashIconify(icon='fluent:border-none-20-regular', width=20), variant="light", color='gray',
                id='BUTTON_clearSelectionCtp_3D', fullWidth=True,
              ), span=4),
              # all selection
              dmc.Col(dmc.Button(
                DashIconify(icon='fluent:checkbox-indeterminate-20-regular', width=20), variant="light", color='gray',
                id='BUTTON_allSelectionCtp_3D', fullWidth=True,
              ), span=4),
            ], gutter=2),
            # tooltips for buttons
            html.Div(
              [
                dbc.Tooltip( i.capitalize(), target=f'BUTTON_{i}SelectionCtp_3D', placement='top')
                for i in ['invert', 'clear', 'all']
              ],
            ),
            dmc.ChipGroup(
              children=[], value=[], multiple=True, align='center', spacing=1, 
              id = 'CHIPGROUP_celltype_3D', className='dmc-ChipGroup-legend'
            ),
            dcc.Store(id='STORE_allCelltypes_3D'),
            fac.AntdDrawer(
              children=[], id='DRAWER_setColorCtp_3D',
              title=dmc.Stack([
                dmc.Text('Setting colors', className='dmc-Text-drawerTitle'),
                dmc.Text("tip: colors will be saved locally", className='dmc-Text-drawerSubTitle')
              ], spacing=1),
              width=300,
            )
          ], span=10)
        ], columns=50),
        # violin
        dbc.Row([
          dbc.Col([
            dbc.Label( 'Normalized expression in all celltypes(left)'),
            dbc.Label('and in each celltype(right):'),
            dcc.Graph(figure={}, id="FIGURE_expViolin_3D", className='dcc-Graph-violin-exp', config=config_violin)
          ], align='center', width=4),
          dbc.Col([
            dcc.Graph(figure={}, id="FIGURE_ctpViolin_3D", className='dcc-Graph-violin-ctp', config=config_violin)
          ], align='center', width=8)
        ], style={'overflow-y': 'auto'}, id='test_sticky')
      ],span=41),
    ], columns=50)],
    label = "Plot feature(3D)",
    tab_id = "spatial_tab_plotFeature3D",
  )

  spatial_tabs = dbc.Tabs(
      children=[spatial_tab_plotFeature3D],
      active_tab = 'spatial_tab_plotFeature3D',
      id='tabs'
    ),

  # In[] callbacks

  # select data -> filter

  @app.callback( # update item options
    Output({'type': 'DROPDOWN_filterItem_3D', 'index': ALL}, 'options'),
    Output({'type': 'SWITCH_filterApply_3D', 'index': ALL}, 'checked'),
    
    Input('DROPDOWN_stage_3D', 'value'),
    State('STORE_filterCurNumber_3D', 'data'),
    # prevent_initial_call=True
  )
  def update_selectData_FilterItemOptions_3D(stage, curNumber):
    obs = exp_data[stage].obs
    options = obs.columns.to_list()
    return [options]*curNumber, [False]*curNumber

  @app.callback( # update item type (numeric/categorical)
    Output({'type': 'DROPDOWN_filterType_3D', 'index': MATCH}, 'value'),
    Output({'type': 'SWITCH_filterApply_3D', 'index': MATCH}, 'checked'),
    
    Input({'type': 'DROPDOWN_filterItem_3D', 'index': MATCH}, 'value'),
    State('DROPDOWN_stage_3D', 'value'),
  )
  def update_selectData_FilterType_3D(item, stage):

    if item is None:
      return no_update, False
    
    dtype = exp_data[stage].obs[item].dtype
    itemType = 'categorical' if ( dtype in [np.dtype('O'), 'category']) else 'numeric'
    
    return itemType, False

  @app.callback( # generate filter body
    Output({'type': 'POPOVER_filterBody_3D', 'index': MATCH}, 'content'),
    Output({'type': 'SWITCH_filterApply_3D', 'index': MATCH}, 'checked'),
    
    State({'type': 'DROPDOWN_filterItem_3D', 'index': MATCH}, 'value'),
    Input({'type': 'DROPDOWN_filterType_3D', 'index': MATCH}, 'value'),
    State({'type': 'DROPDOWN_filterType_3D', 'index': MATCH}, 'id'), # MATCH <-> id['index']
    State('DROPDOWN_stage_3D', 'value'),
    # prevent_initial_call=True,
    # _allow_dynamic_callbacks=True, # not suitable for multi-user/multi-process apps, which will cause growing dynamic callbacks
  )
  def update_selectData_FilterContainer_3D(item, itemType, id, stage):
    
    if item is None:
      raise PreventUpdate
    
    obs = exp_data[stage].obs
    
    index = id['index']
    
    if itemType == 'numeric':
      component = html.Div(
        [
          dmc.Grid(
            [
              dmc.Col(dmc.Text(f'{item}'), span=6),
              dmc.Col(dmc.Text(' ≥ '), span=2),
              dmc.Col(
                dmc.NumberInput(precision=4, step=0.1, id={'type': 'NUMBERINPUT_filterBodyNumberLeft', 'index': index}),
                span=4
              ),
              dmc.Col(dmc.Text(f'{item}'), span=6),
              dmc.Col(dmc.Text(' ≤ '), span=2),
              dmc.Col(
                dmc.NumberInput(precision=4, step=0.1, id={'type': 'NUMBERINPUT_filterBodyNumberRight', 'index': index}),
                span=4
              ),
            ],
          ),
          # solve nonexistent id problem(callback)
          html.Div(id={'type': 'TRANSFER_filterBodyCategorical_3D', 'index': index})
        ],
        className='div-filterBody-numeric'
      )
      
    elif itemType == 'categorical':
      options = [{'value': i, 'label': i} for i in obs[item].unique()]
      component = html.Div(
        [
          dmc.TransferList(
            value=[options, []], titles=['Unselected', 'Selected'],
            id={'type': 'TRANSFER_filterBodyCategorical_3D', 'index': index}
          ),
          # solve nonexistent id problem(callback)
          html.Div(id={'type': 'NUMBERINPUT_filterBodyNumberLeft', 'index': index}),
          html.Div(id={'type': 'NUMBERINPUT_filterBodyNumberRight', 'index': index}),
        ]
      )
      
    else:
      raise PreventUpdate
    
    return component, False

  @app.callback( # calculate preserved cells
    Output({'type': 'STORE_filterBodyPreservedCells_3D', 'index': MATCH}, 'data'),
    Output({'type': 'TEXT_filterBodyPreservedCells_3D', 'index': MATCH}, 'children'),
    
    Input({'type': 'TRANSFER_filterBodyCategorical_3D', 'index': MATCH}, 'value'),
    Input({'type': 'NUMBERINPUT_filterBodyNumberLeft', 'index': MATCH}, 'value'),
    Input({'type': 'NUMBERINPUT_filterBodyNumberRight', 'index': MATCH}, 'value'),
    State('DROPDOWN_stage_3D', 'value'),
    State({'type': 'DROPDOWN_filterItem_3D', 'index': MATCH}, 'value'),
    State({'type': 'DROPDOWN_filterType_3D', 'index': MATCH}, 'value'),
    Input({'type': 'SWITCH_filterApply_3D', 'index': MATCH}, 'checked'),
    prevent_initial_call=True
  )
  def store_categoricalObsFilterInfo_forFilterBody_3D(transfer_value, numberLeft, numberRight, stage, item, itemType, checked):
    
    if item is None:
      return None, 'Selected cells:' # clear filter when items are deleted

    item_column = exp_data[stage].obs[item]
    
    if itemType == 'categorical':
      if checked:
        selected =[ d['value'] for d in transfer_value[1]]
        cells = item_column[item_column.isin(selected)].index.to_list()
        text = f'Selected cells: {len(cells)}'
      else:
        cells = None
        text = 'Selected cells:'
      return Serverside(cells), text

    elif itemType == 'numeric':
      if checked and (numberLeft!='' or numberRight!=''):
        boolLeft = item_column >= numberLeft if numberLeft!='' else [True]*len(item_column)
        boolRight = item_column <= numberRight if numberRight!='' else [True]*len(item_column)
        cells = item_column[boolLeft & boolRight].index.to_list()
        text = f'Selected cells: {len(cells)}'
      else:
        cells = None
        text = 'Selected cells:'
      return Serverside(cells), text
    
    else:
      raise PreventUpdate

  @app.callback( # update transfer value (stage changed, categorical)
    Output({'type': 'TRANSFER_filterBodyCategorical_3D', 'index': ALL}, 'value'),
    Input('DROPDOWN_stage_3D', 'value'),
    State({'type': 'DROPDOWN_filterItem_3D', 'index': ALL}, 'value'),
    State({'type': 'TRANSFER_filterBodyCategorical_3D', 'index': ALL}, 'id'),
  )
  def update_categoricalTransferValue_forFilterBody_3D(stage, item_list, id_list):
    
    indices = [ id['index'] for id in id_list ]
    items = [item_list[i] for i in indices]

    value_list = [
      [[{'value': i, 'label': i} for i in exp_data[stage].obs[item].unique()], []]
      for item in items
    ]

    return value_list

  @app.callback( # intersection of filter-cells
    Output('STORE_cellsObsFilter_3D', 'data'),
    Input({'type': 'STORE_filterBodyPreservedCells_3D', 'index': ALL}, 'data'),
    prevent_initial_call=True,
  )
  def update_cellsObsFilterStore_forJSONtoPlot(cells_list):
    
    def func(a,b):
      if (a is not None) and (b is None):
        return set(a)
      elif (a is None) and (b is not None):
        return set(b)
      elif (a is not None) and (b is not None):
        return set(a) & set(b)
      else:
        return None
    
    if len(cells_list) > 1:
      cells = reduce(func, cells_list)
    else:
      cells = cells_list[0]
    # at least one filter with cells: [None]
    
    return Serverside(cells)

  @app.callback( # add & delete selectData-filter
    Output('CONTAINER_filterList_3D', 'children'),
    Output('STORE_filterCurNumber_3D', 'data'),
    Input('BUTTON_addFilter_3D', 'n_clicks'),
    Input('BUTTON_deleteFilter_3D', 'n_clicks'),
    State('STORE_filterCurNumber_3D', 'data'),
    prevent_initial_call=True,
  )
  def addAndDelte_forSeleteDataFilter_3D(add, delete, curNumber):
    
    tid = ctx.triggered_id
    nextIndex = curNumber
    patch_children = Patch()
    
    if 'BUTTON_addFilter_3D' in tid:
      patch_children.append(selectData_newFilter(index=nextIndex))
      nextNumber = curNumber+1
    elif 'BUTTON_deleteFilter_3D' in tid:
      if nextIndex >= 2:
        del patch_children[nextIndex-1]
        nextNumber = curNumber-1 if curNumber>0 else 0
      else:
        raise PreventUpdate
    else:
      raise PreventUpdate
    
    return patch_children, nextNumber

  # download figure format & resolution
  @app.callback(
    Output('FIGURE_3Dcelltype', 'config'),
    Output('FIGURE_3Dexpression', 'config'),
    Input('NUMBERINPUT_scatter3dFigtype_3D', 'value'),
    Input('NUMBERINPUT_scatter3dFigscale_3D', 'value'),
  )
  def update_scatter3dDownloadConfig_3D(type, scale):
    
    patch=Patch()
    patch['toImageButtonOptions']['format'] = type
    patch['toImageButtonOptions']['scale'] = scale
    
    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'config'),
    Output('FIGURE_ctpViolin_3D', 'config'),
    Input('NUMBERINPUT_violinFigtype_3D', 'value'),
    Input('NUMBERINPUT_violinFigscale_3D', 'value'),
  )
  def update_violinDownloadConfig_3D(type, scale):
    
    patch=Patch()
    patch['toImageButtonOptions']['format'] = type
    patch['toImageButtonOptions']['scale'] = scale
    
    return patch, patch

  # violin options hot-update
  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('SEGMENTEDCONTROL_violinPoints_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinPointStyle_3D(points, allCelltypes, minfo):
    
    points = False if points=='none' else points
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['points'] = points

    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('NUMBERINPUT_violinPointpos_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinPointpos_3D(pointpos, allCelltypes, minfo):
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['pointpos'] = pointpos

    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('NUMBERINPUT_violinPointsize_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinPointsize_3D(pointsize, allCelltypes, minfo):
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['marker']['size'] = pointsize

    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('SEGMENTEDCONTROL_violinBox_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinBox_3D(box, allCelltypes, minfo):
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    box_visible = True if box=='box' or box=='all' else False
    meanline_visible = True if box=='meanline' or box=='all' else False
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['box']['visible'] = box_visible
      patch['data'][i]['meanline']['visible'] = meanline_visible

    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('NUMBERINPUT_violinPointjitter_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinPointpos_3D(jitter, allCelltypes, minfo):
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['jitter'] = jitter

    return patch, patch

  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('NUMBERINPUT_violinBoxwidth_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    State('STORE_multiNameInfo_3D', 'data'),
  )
  def update_violinPointpos_3D(boxwidth, allCelltypes, minfo):
    
    n_gene = len(minfo)
    n_ctp = len(allCelltypes)
    
    patch = Patch()
    for i in range(0, max(n_gene, n_ctp)):
      patch['data'][i]['box']['width'] = boxwidth

    return patch, patch

  # update_dataSummary
  @app.callback(
    Output('TEXT_dataSummary_3D', 'children'),
    Input('DROPDOWN_stage_3D', 'value')
  )
  def update_dataSummary_3D(stage):
    adata = exp_data[stage]
    str = f'{adata.shape[0]}(cells) × {adata.shape[1]}(features)'
    return str

  # update_nameOptions
  @app.callback(
    Output('DROPDOWN_singleName_3D', 'options'),
    Input('DROPDOWN_singleName_3D', 'search_value'),
    Input('DROPDOWN_stage_3D', 'value')
  )
  def update_nameOptions_single_3D(search, stage):
    if not search:
      raise PreventUpdate

    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()

    return opts

  @app.callback(
    Output({'type': 'DROPDOWN_multiName_3D', 'index': MATCH}, 'options'),
    Input({'type': 'DROPDOWN_multiName_3D', 'index': MATCH}, 'search_value'),
    Input('DROPDOWN_stage_3D', 'value'),
    prevent_initial_call=True,
  )
  def update_nameOptions_multi_3D(search,stage):
    if not search:
      raise PreventUpdate
    
    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
    
    return opts

  # add & delte components for multiName
  @app.callback(
    Output('DIV_multiNameDynamic_3D', 'children'),
    Output('STORE_multiNameCurNumber', 'data'),
    Input('BUTTON_addFeature_3D', 'n_clicks'),
    Input('BUTTON_deleteFeature_3D', 'n_clicks'),
    State('STORE_multiNameCurNumber', 'data'),
    State('DROPDOWN_stage_3D', 'value'),
    prevent_initial_call = True,
  )
  def add_components_multiName_3D(add, delete, curNumber, stage):

    id = ctx.triggered_id

    nextIndex = curNumber
    nextColor = initColor_multiName[int(nextIndex % len(initColor_multiName))]
    
    patch_children = Patch()
    if 'BUTTON_addFeature_3D' in id:
      patch_children.append(
        dmc.Grid([
          dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': nextIndex}), span=10),
          dmc.Col(
            iconHover_colorPicker(
              id={
                'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': nextIndex}, 
                'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': nextIndex}, 
                'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': nextIndex},
                },
              init_color=nextColor, swatches=colorPicker_swatches,
            ),
            span=2
          ),
        ])
      )
      nextNumber = curNumber+1
    elif 'BUTTON_deleteFeature_3D' in id:
      if nextIndex >= 3 :
        del patch_children[nextIndex-1]
        nextNumber = curNumber-1 if curNumber>0 else 0
      else:
        nextNumber = curNumber

    return patch_children, nextNumber

  # store_previewRange
  app.clientside_callback(
    ClientsideFunction(
      namespace='plotFunc_3Dtab',
      function_name='store_previewRange',
    ),
    Output('STORE_previewRange_3D', 'data'),
    Input('SLIDER_Xrange_3D', 'value'),
    Input('SLIDER_Yrange_3D', 'value'),
    Input('SLIDER_Zrange_3D', 'value'),
  )

  # store_sliceRange
  app.clientside_callback(
    ClientsideFunction(
      namespace='plotFunc_3Dtab',
      function_name='store_sliceRange'),
    Output('STORE_sliceRange_3D', 'data'),
    Input('BUTTON_slice_3D', 'n_clicks'),
    Input('BUTTON_recover_3D', 'n_clicks'),
    Input('STORE_maxRange_3D', 'data'),
    State('STORE_previewRange_3D', 'data'),
  )

  # max range
  @app.callback(
    Output('STORE_maxRange_3D', 'data'),
    Input('DROPDOWN_stage_3D', 'value'),
  )
  def update_maxRange_3D(stage):
    obs = exp_data[stage].obs
    maxRange = dict(
      x_min = np.floor(obs.x.min()/10)*10, x_max = np.ceil(obs.x.max()/10)*10,
      y_min = np.floor(obs.y.min()/10)*10, y_max = np.ceil(obs.y.max()/10)*10,
      z_min = np.floor(obs.z.min()/10)*10, z_max = np.ceil(obs.z.max()/10)*10,
    )
    return maxRange

  @app.callback(
    output=[
      ( Output('SLIDER_Xrange_3D', 'min'), Output('SLIDER_Xrange_3D', 'max'), Output('SLIDER_Xrange_3D', 'value') ),
      ( Output('SLIDER_Yrange_3D', 'min'), Output('SLIDER_Yrange_3D', 'max'), Output('SLIDER_Yrange_3D', 'value') ),
      ( Output('SLIDER_Zrange_3D', 'min'), Output('SLIDER_Zrange_3D', 'max'), Output('SLIDER_Zrange_3D', 'value') ),
    ],
    inputs = Input('STORE_maxRange_3D', 'data'),
  )
  def update_sliderRange_3D(maxRange):
    
    res = [
      ( maxRange[f'{c}_min'], maxRange[f'{c}_max'], (maxRange[f'{c}_min'], maxRange[f'{c}_max']) )
      for c in ['x', 'y', 'z']
    ]
    
    return  res

  # store_cellsObs_forJSONtoPlot
  @app.callback(
    Output('STORE_obs_3D', 'data'),
    Input('DROPDOWN_stage_3D', 'value'),
  )
  def store_cellsObs_forJSONtoPlot_3D(stage):
    
    adata = exp_data[stage]
    
    obs = adata.obs[['x','y','z','celltype']]
    return obs.to_dict('index')

  # store sliceInfo for JSONtoPLot
  @app.callback(
    Output('STORE_cellsSliceFilter_3D', 'data'),
    
    Input('STORE_sliceRange_3D', 'data'),
    Input('DROPDOWN_stage_3D', 'value'),
    
    Trigger('BUTTON_slice_3D', 'n_clicks'),
    prevent_initial_call=True
  )
  def store_sliceInfo_forJSONtoPlot_3D(sliceRange, stage):

    adata = exp_data[stage]
    
    obs = adata.obs[['x','y','z']]

    tid = ctx.triggered_id
    if tid and 'BUTTON_slice_3D' in tid:
      if_inSliceRange = ( 
        (obs['x'] <= sliceRange['x_max']) & 
        (obs['x'] >= sliceRange['x_min']) & 
        (obs['y'] <= sliceRange['y_max']) & 
        (obs['y'] >= sliceRange['y_min']) & 
        (obs['z'] <= sliceRange['z_max']) & 
        (obs['z'] >= sliceRange['z_min'])
      )
      obsnames_filt = adata.obs_names[if_inSliceRange]
    else:
      obsnames_filt = adata.obs_names
    
    return Serverside(obsnames_filt.to_list())

  # store_expInfo_forJSONtoPlot
  @app.callback(
    Output('STORE_cellsExpFilter_3D', 'data'),
    Output('STORE_singleExp_3D', 'data'),
    Output('STORE_multiExp_3D', 'data'),
    Output('STORE_ifmulti_3D', 'data'),
    Output('STORE_mixedColor_3D', 'data'),
    
    Input('BUTTON_singlePlot_3D', 'n_clicks'),
    Input('BUTTON_multiPlot_3D', 'n_clicks'),
    Input('SWITCH_hideZero_3D', 'checked'),
    Input('DROPDOWN_stage_3D', 'value'),
    
    State('DROPDOWN_singleName_3D', 'value'),
    State('STORE_multiNameInfo_3D', 'data'),
    State('STORE_ifmulti_3D', 'data'),
  )
  def store_expInfo_forJSONtoPlot_3D(sclick, mclick, hideZero, stage, sname, minfo, ifmulti):

    adata = exp_data[stage]

    def return_single():
      ifmulti = False
      exp = adata[:,sname].to_df()
      if hideZero:
        cellsExpFilter = exp[(exp>0)[sname]].index
      else:
        cellsExpFilter = exp.index
      exp = exp.loc[cellsExpFilter,:]
      cellsExpFilter = cellsExpFilter.to_list()
      return (ifmulti, exp, cellsExpFilter)
    
    def return_multi():
      ifmulti = True
      mixColor = color_mixer(adata, minfo)
      if hideZero:
        cellsExpFilter = mixColor[mixColor!='rgb(244, 244, 244)'].index
      else:
        cellsExpFilter = mixColor.index
      mixColor = mixColor[cellsExpFilter]
      cellsExpFilter = cellsExpFilter.to_list()
      return (ifmulti, [], cellsExpFilter, mixColor.to_dict()) 
    
    def return_multiExp():
      tmp = {}
      for key,value in minfo.items():
        if value:
          tmp[key] = value
      colors = list(tmp.keys())
      genes = list(tmp.values())

      exp = adata[:, genes].to_df()
      exp.columns = colors
      exp = exp.to_dict('index')

      return exp
    
    btn_id = ctx.triggered_id
    if btn_id:
      if 'DROPDOWN_stage_3D' in btn_id:
        if not ifmulti:
          ifmulti,exp,cellsExpFilter = return_single()
          exp = exp.to_dict('index')
          return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)
        else:
          ifmulti,_,cellsExpFilter,mixcolor = return_multi()
          exp_multi = return_multiExp()
          return (Serverside(cellsExpFilter), no_update, exp_multi, ifmulti, mixcolor)

      elif 'BUTTON_singlePlot_3D' in btn_id:
        ifmulti,exp,cellsExpFilter = return_single()
        exp = exp.to_dict('index')
        if hideZero:
          return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)
        else:
          return (no_update, exp, no_update, ifmulti, no_update)
      
      elif 'BUTTON_multiPlot_3D' in btn_id:
        ifmulti,_,cellsExpFilter,mixcolor = return_multi()
        exp_multi = return_multiExp()
        if hideZero:
          return (Serverside(cellsExpFilter), no_update, exp_multi, ifmulti, mixcolor)
        else:
          return (no_update, no_update, exp_multi, ifmulti, mixcolor)
      
      elif 'SWITCH_hideZero_3D' in btn_id:
        
        if not hideZero:
          cellsExpFilter = adata.obs_names.to_list()
          return (Serverside(cellsExpFilter), no_update, no_update, no_update, no_update)
        
        else:
          if not ifmulti:
            _,_,cellsExpFilter = return_single()
            return (Serverside(cellsExpFilter), no_update, no_update, no_update, no_update)
          else:
            _,_,cellsExpFilter,_ = return_multi()
            exp_multi = return_multiExp()
            return (Serverside(cellsExpFilter), no_update, exp_multi, no_update, no_update)

    else:
        ifmulti,exp,cellsExpFilter = return_single()
        exp = exp.to_dict('index')
        return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)

  # update ChipGroup-celltype chips
  @app.callback(
    Output('CHIPGROUP_celltype_3D', 'children'),
    Output('CHIPGROUP_celltype_3D', 'value'),
    Output('STORE_allCelltypes_3D', 'data'),
    Input('STORE_cellsSliceFilter_3D', 'data'),
    Input('STORE_cellsExpFilter_3D', 'data'),
    State('DROPDOWN_stage_3D', 'value'),
    State('STORE_ctpCmap_3D', 'data'),
  )
  def update_chipGroupCelltype_3D(obsFilter, expFilter, stage, cmap):
    cells = list(set(obsFilter)&set(expFilter))
    cells.sort()
    celltypes = list(exp_data[stage].obs['celltype'][cells].unique())
    
    chips = [
      dmc.Chip(
        children=ctp, value=ctp, size='xs', color='gray', variant='filled', type='radio',
        styles = {
          'label': {
            'color': cmap[ctp],
            'font-size': 12,
            'font-weight': 'bold'
          },
          'checkIcon': {
            'color': cmap[ctp],
          }
        },
        id = {'type': 'CHIP_ctpColorLegend_3D', 'id': ctp}
      ) 
      for ctp in celltypes
    ]
    
    return chips, celltypes, celltypes

  @app.callback(
    Output('CHIPGROUP_celltype_3D', 'value'),
    Input('BUTTON_invertSelectionCtp_3D', 'n_clicks'),
    State('CHIPGROUP_celltype_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    prevent_initial_call=True,
  )
  def invertSelection_celltypesButton_3D(click, curValue, allCelltypes):
    return list(set(allCelltypes) - set(curValue))

  @app.callback(
    Output('CHIPGROUP_celltype_3D', 'value'),
    Input('BUTTON_clearSelectionCtp_3D', 'n_clicks'),
    prevent_initial_call=True
  )
  def clearSelection_celltypesButton_3D(click):
    return []

  @app.callback(
    Output('CHIPGROUP_celltype_3D', 'value'),
    Input('BUTTON_allSelectionCtp_3D', 'n_clicks'),
    State('CHIPGROUP_celltype_3D', 'value'),
    State('STORE_allCelltypes_3D', 'data'),
    prevent_initial_call=True
  )
  def allSelection_celltypesButton_3D(click, curValue, allCelltypes):
    if set(curValue) == set(allCelltypes):
      return no_update
    else:
      return list(set(allCelltypes))

  # store_ctpInfo_forJSONtoPLot
  @app.callback(
    Output('STORE_cellsCtpFilter_3D', 'data'),
    Input('CHIPGROUP_celltype_3D', 'value'),
    State('DROPDOWN_stage_3D', 'value')
  )
  def store_ctpInfo_forJSONtoPlot_3D(selectedCtps, stage):
      
    series = exp_data[stage].obs['celltype']
    series = series[series.isin(selectedCtps)]
    
    return Serverside(series.index.to_list())
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
  # plot_3Dfigure_exp
  app.clientside_callback(
    ClientsideFunction(
      namespace='plotFunc_3Dtab',
      function_name='exp_3Dscatter',
    ),
    Output("FIGURE_3Dexpression", "figure"),
    Input('STORE_obs_3D', 'data'),
    Input('STORE_cellsIntersection_3D', 'data'),
    Input('STORE_singleExp_3D', 'data'),
    Input('STORE_ifmulti_3D', 'data'),
    Input('STORE_mixedColor_3D', 'data'),
    State('SWITCH_hideAxes_3D', 'checked'),
    State('SWITCH_previewBox_3D', 'checked'),
    State('STORE_previewRange_3D', 'data'),
    State('SEGMENTEDCONTROL_projection_3D', 'value'),
    State('NUMBERINPUT_scatter3dPointsize_3D', 'value')
  )

  # colorpicker for singleExp
  @app.callback(
    Output('ACTIONICON_colorSingle_3D', 'children'),
    Output('FIGURE_3Dexpression', 'figure'),
    Input('COLORPICKER_single_3D', 'value'),
  )
  def colorpicker_for_singleExp_3D(color):
    patch = Patch()
    patch['layout']['coloraxis']['colorscale'][1][1] = color
    icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
    return icon, patch

  @app.callback(
    Output('TEXT_colorSingle_3D', 'value'),
    Output('COLORPICKER_single_3D', 'value'),
    Input('TEXT_colorSingle_3D', 'value'),
    Input('COLORPICKER_single_3D', 'value'),
    prevent_initial_call=True,
  )
  def linkage_colorPickerAndTextSingle_3D(value1, value2):
    id = ctx.triggered_id
    if id == 'TEXT_colorSingle_3D':
      if((len(value1)==4) or (len(value1)==7)):
        return no_update, value1
      else:
        raise PreventUpdate 
    else:
      return value2, no_update

  # colopicker for multiExp

  @app.callback(
    Output('STORE_multiNameInfo_3D', 'data'),
    Input({'type': 'COLORPICKER_multi_3D', 'index': ALL}, 'value'),
    Input({'type': 'DROPDOWN_multiName_3D', 'index': ALL}, 'value'),
  )
  def store_multiNameInfo_3D(colors, genes):
    return dict(zip(colors, genes))

  @app.callback(
    Output({'type':'ACTIONICON_colorMulti_3D', 'index': MATCH}, 'children'),
    Output({'type': 'TEXT_colorMulti_3D', 'index': MATCH}, 'value'),
    Output({'type': 'COLORPICKER_multi_3D', 'index': MATCH}, 'value'),
    Input({'type': 'TEXT_colorMulti_3D', 'index': MATCH}, 'value'),
    Input({'type': 'COLORPICKER_multi_3D', 'index': MATCH}, 'value'),
    prevent_initial_call=True,
  )
  def linkage_colorPickerAndTextMulti_3D(value1, value2):
    id = ctx.triggered_id
    
    if id['type'] == 'TEXT_colorSingle_3D':
      if((len(value1)==4) or (len(value1)==7)):
        color = value1
        icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
        return icon, no_update, color
      else:
        raise PreventUpdate   
    else:
      color = value2
      icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
      return icon, color, no_update

  # colorpicker for ctpLegend
  @app.callback(
    Output('DRAWER_setColorCtp_3D', 'visible'),
    Input('BUTTON_setColors_3D', 'n_clicks'),
    prevent_initial_call=True,
  )
  def setCelltypeColorsInDrawer_3D(click):
    return True

  @app.callback(
    Output('DRAWER_setColorCtp_3D', 'children'),
    Input('STORE_allCelltypes_3D', 'data'),
    State('STORE_ctpCmap_3D', 'data'),
    # prevent_initial_call=True, 
  )
  def generate_drawerLegendContent_3D(curAllCtps, cmap):
    return drawerContent_ctpColorPicker(curAllCtps, cmap)

  @app.callback(
    Output({'type': 'ACTIONICON_colorCtp_3D', 'id': MATCH}, 'children'),
    Output({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
    Output({'type': 'COLORPICKER_colorCtp_3D', 'id': MATCH}, 'value'),
    Input({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
    Input({'type': 'COLORPICKER_colorCtp_3D', 'id': MATCH}, 'value'),
    prevent_initial_call=True,
  )
  def syncAndReturn_colorValue_3D(text, picker):


    tid = ctx.triggered_id
    celltype = tid['id']
    
    if tid['type'] == 'TEXT_colorCtp_3D':
      if((len(text)==4) or (len(text)==7)):
        color = text  
        icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
        return icon, no_update, color
      else:
        raise PreventUpdate
    else:
      color = picker
      icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
      return icon, color, no_update
    
  @app.callback(
    Output('STORE_ctpCmap_3D', 'data'),
    Input({'type': 'TEXT_colorCtp_3D', 'id': ALL}, 'value'),
    prevent_initial_call=True
  )
  def update_storeCtpCmap_3D(colors):
      triggered = ctx.triggered
      triggered_id = ctx.triggered_id
      if(len(triggered) > 1):
        raise PreventUpdate
      
      color = triggered[0]['value']
      ctp = triggered_id['id']
      patch = Patch()
      patch[ctp] = color
      return patch

  @app.callback(
    Output('FIGURE_3Dcelltype', 'figure'),
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('STORE_ctpCmap_3D', 'data'),
    State('STORE_allCelltypes_3D', 'data'),
    prevent_initial_call=True,
  )
  def update_figureCtpCmap_3D(cmap, curCtps):
    
    patch_fig=Patch()
    for i in range(0, len(curCtps)):
      patch_fig['data'][i]['marker']['color'] =  cmap[curCtps[i]]
    return patch_fig, patch_fig

  @app.callback(
    Output({'type': 'CHIP_ctpColorLegend_3D', 'id': MATCH}, 'styles'),
    Input({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
    prevent_initial_call=True
  )
  def update_chipColor_3D(color):
    patch = Patch()
    patch['label']['color'] = color
    patch['checkIcon']['color'] = color
    return patch

  # plot_3Dfigure_ctp
  app.clientside_callback(
    ClientsideFunction(
      namespace='plotFunc_3Dtab',
      function_name='ctp_3Dscatter',
    ),
    Output("FIGURE_3Dcelltype", "figure"),
    Input('STORE_obs_3D', 'data'),
    Input('STORE_cellsIntersection_3D', 'data'),
    State('SWITCH_hideAxes_3D', 'checked'),
    State('SWITCH_previewBox_3D', 'checked'),
    State('STORE_previewRange_3D', 'data'),
    State('STORE_ctpCmap_3D', 'data'),
    State('SEGMENTEDCONTROL_projection_3D', 'value'),
    State('NUMBERINPUT_scatter3dPointsize_3D', 'value')
  )

  # sync layout between exp and ctp figure
  @app.callback(
    Output("FIGURE_3Dexpression", "figure"),
    Output("FIGURE_3Dcelltype", "figure"),
    Input("FIGURE_3Dexpression", "relayoutData"),
    Input("FIGURE_3Dcelltype", "relayoutData"),
    State('SEGMENTEDCONTROL_projection_3D', 'value'),
    # prevent_initial_call=True,
    # background=True,
    # manager=background_callback_manager
  )
  def update_relayout(expLayout, ctpLayout, proj):
    tid = ctx.triggered_id
    patch = Patch()
    
    if tid == 'FIGURE_3Dexpression':
      
      if 'scene.camera' in expLayout:
        patch['layout']['scene']['camera'] = expLayout['scene.camera']
      if 'scene.aspectratio' in expLayout:
        patch['layout']['scene']['aspectmode'] = 'manual'
        patch['layout']['scene']['aspectratio'] = expLayout['scene.aspectratio']

      return patch, patch

    elif tid == 'FIGURE_3Dcelltype':
      
      if 'scene.camera' in ctpLayout:
        patch['layout']['scene']['camera'] = ctpLayout['scene.camera']
      if 'scene.aspectratio' in ctpLayout:
        patch['layout']['scene']['aspectmode'] = 'manual'
        patch['layout']['scene']['aspectratio'] = ctpLayout['scene.aspectratio']

      return patch, patch
    
    else:
      raise PreventUpdate

  # update scatter-3d point size
  @app.callback(
    Output('FIGURE_3Dexpression', 'figure'),
    Input('NUMBERINPUT_scatter3dPointsize_3D', 'value'),
    prevent_initial_call = True
  )
  def update_expPointSize_3D(size):
    
    patch = Patch()
    patch['data'][0]['marker']['size'] = size
    
    return patch

  @app.callback(
    Output('FIGURE_3Dcelltype', 'figure'),
    Input('NUMBERINPUT_scatter3dPointsize_3D', 'value'),
    State('STORE_cellsIntersection_3D', 'data'),
    State('DROPDOWN_stage_3D', 'value'),
    prevent_initial_call = True,
  )
  def update_ctpPointSize_3D(size, cells, stage):
    
    adata = exp_data[stage]
    celltypes = adata.obs.loc[cells, 'celltype'].unique()
    patch = Patch()
    for i in range(0,len(celltypes)):
      patch['data'][i]['marker']['size'] = size
    return patch

  # switch projection type
  @app.callback(
    Output('FIGURE_3Dcelltype', 'figure'),
    Output('FIGURE_3Dexpression', 'figure'),
    Input('SEGMENTEDCONTROL_projection_3D', 'value'),
  )
  def switch_projectionType(type):
    patch=Patch()
    patch['layout']['scene']['camera']['projection'] = {'type': type}
    return patch, patch

  # find intersection of all cellsFilters
  @app.callback(
    Output('STORE_cellsIntersection_3D', 'data'),
    Input('STORE_cellsSliceFilter_3D', 'data'),
    Input('STORE_cellsExpFilter_3D', 'data'),
    Input('STORE_cellsCtpFilter_3D', 'data'),
    Input('STORE_cellsObsFilter_3D', 'data'),
    prevent_initial_call=True,
  )
  def intersection_of_filter(sliceFilter, expFilter, ctpFilter, obsFilter):
    filter_list = [ l for l in [sliceFilter, expFilter, ctpFilter, obsFilter] if l]
    
    if len(filter_list) == 0:
      tmp = []
    elif len(filter_list) == 1 :
      tmp = filter_list[0]
    elif len(filter_list) > 1 :
      tmp = list(reduce(lambda a,b: set(a) & set(b), filter_list))
    tmp.sort()
    
    return tmp

  # hide axes
  @app.callback(
    Output("FIGURE_3Dexpression", "figure", allow_duplicate=True),
    Output("FIGURE_3Dcelltype", "figure", allow_duplicate=True),
    Input('SWITCH_hideAxes_3D', 'checked'),
    prevent_initial_call=True
  )
  def hideAxes_3D(hideAxes):
    patch = Patch()
    if hideAxes:
      patch['layout']['scene']['xaxis']['visible'] = False
      patch['layout']['scene']['yaxis']['visible'] = False
      patch['layout']['scene']['zaxis']['visible'] = False
    else: 
      patch['layout']['scene']['xaxis']['visible'] = True
      patch['layout']['scene']['yaxis']['visible'] = True
      patch['layout']['scene']['zaxis']['visible'] = True
    return patch, patch

  # show preview Box
  @app.callback(
    Output('FIGURE_3Dexpression', 'figure', allow_duplicate=True),
    Output('FIGURE_3Dcelltype', 'figure', allow_duplicate=True),
    Input('SWITCH_previewBox_3D', 'checked'),
    Input('STORE_previewRange_3D', 'data'),
    prevent_initial_call=True,
  )
  def update_previewBox(showBox, preRange):
    patch = Patch()
    if showBox:
      patch['data'][-1] = {
                      'x': [preRange['x_min'], preRange['x_min'], preRange['x_min'], preRange['x_min'],
                            preRange['x_max'], preRange['x_max'], preRange['x_max'], preRange['x_max']],
                      'y': [preRange['y_min'], preRange['y_max'], preRange['y_min'], preRange['y_max'],
                            preRange['y_min'], preRange['y_max'], preRange['y_min'], preRange['y_max']],
                      'z': [preRange['z_min'], preRange['z_min'], preRange['z_max'], preRange['z_max'],
                            preRange['z_min'], preRange['z_min'], preRange['z_max'], preRange['z_max']],
                      'i': [0, 1, 0, 0, 0, 0, 2, 2, 7, 7, 7, 7],
                      'j': [1, 2, 4, 1, 4, 2, 3, 6, 4, 4, 1, 1],
                      'k': [2, 3, 5, 5, 6, 6, 7, 7, 6, 5, 3, 5],
                      'color': 'black', 'opacity': 0.60, 'type': 'mesh3d'
                    }
    else:
      patch['data'][-1] = {
                      'x': [], 'y': [], 'z': [], 'i': [], 'j': [], 'k': [],
                      'type': 'mesh3d', 'color': 'black', 'opacity': 0.60
                    }

    return patch, patch

  # violin plot
  @app.callback(
    Output('FIGURE_expViolin_3D', 'figure'),
    Input('DROPDOWN_stage_3D', 'value'),
    Input('STORE_cellsIntersection_3D', 'data'),
    Input('STORE_ifmulti_3D', 'data'),
    Input('BUTTON_singlePlot_3D', 'n_clicks'),
    Input('BUTTON_multiPlot_3D', 'n_clicks'),
    State('DROPDOWN_singleName_3D', 'value'),
    State('STORE_multiNameInfo_3D', 'data'),
    State('SEGMENTEDCONTROL_violinPoints_3D', 'value'),
    State('NUMBERINPUT_violinPointpos_3D', 'value'),
    State('NUMBERINPUT_violinPointsize_3D', 'value'),
    State('NUMBERINPUT_violinPointjitter_3D', 'value'),
    State('SEGMENTEDCONTROL_violinBox_3D', 'value'),
    State('NUMBERINPUT_violinBoxwidth_3D', 'value'),
    # background = True,
    # manager = background_callback_manager,
  )
  def update_spatial_plotFeature3D_expViolin(stage, cells, ifmulti, splot, mplot, sname, minfo, 
                                            points, pointpos, pointsize,jitter, box, boxwidth):
    
    adata = exp_data[stage]

    adata = adata[cells]

    points = False if points=='none' else points
    
    box_visible = True if box=='box' or box=='all' else False
    meanline_visible = True if box=='meanline' or box=='all' else False

    if not ifmulti:
      fig = show_expViolin(adata, sname, points=points, pointpos=pointpos, marker_size=pointsize, 
                          meanline_visible=meanline_visible,  box_visible=box_visible, jitter=jitter, box_width=boxwidth)
    else:
      fig = show_multiFeatures_expViolin(adata, minfo, points=points, pointpos=pointpos, marker_size=pointsize, 
                                        meanline_visible=meanline_visible,  box_visible=box_visible, jitter=jitter, box_width=boxwidth)

    return fig

  @app.callback(
    Output('FIGURE_ctpViolin_3D', 'figure'),
    Input('DROPDOWN_stage_3D', 'value'),
    Input('STORE_cellsIntersection_3D', 'data'),
    Input('STORE_ifmulti_3D', 'data'),
    Input('BUTTON_singlePlot_3D', 'n_clicks'),
    Input('BUTTON_multiPlot_3D', 'n_clicks'),
    State('DROPDOWN_singleName_3D', 'value'),
    State('STORE_multiNameInfo_3D', 'data'),
    State('SEGMENTEDCONTROL_violinPoints_3D', 'value'),
    State('NUMBERINPUT_violinPointpos_3D', 'value'),
    State('NUMBERINPUT_violinPointsize_3D', 'value'),
    State('NUMBERINPUT_violinPointjitter_3D', 'value'),
    State('SEGMENTEDCONTROL_violinBox_3D', 'value'),
    State('NUMBERINPUT_violinBoxwidth_3D', 'value'),
    # background = True,
    # manager = background_callback_manager,
  )
  def update_spatial_plotFeature3D_ctpExpViolin(stage, cells, ifmulti, splot, mplot, sname, minfo, 
                                                points, pointpos, pointsize, jitter, box, boxwidth):
    adata = exp_data[stage]
    adata = adata[cells]

    points = False if points=='none' else points
    
    box_visible = True if box=='box' or box=='all' else False
    meanline_visible = True if box=='meanline' or box=='all' else False

    if not ifmulti:
      fig = show_ctpExpViolin(adata, sname, points=points, pointpos=pointpos, marker_size=pointsize, 
                              meanline_visible=meanline_visible, box_visible=box_visible, jitter=jitter, box_width=boxwidth)
    else:
      fig = show_multiFeatures_ctpExpViolin(adata, minfo, points=points, pointpos=pointpos, marker_size=pointsize, 
                                            meanline_visible=meanline_visible, box_visible=box_visible, jitter=jitter, box_width=boxwidth)

    return fig


  # moran SVG offcanvas
  @app.callback(
    Output('OFFCANVAS_moranRes_3D', 'is_open'),
    Input('BUTTON_showMoran_3D', 'n_clicks'),
    prevent_initial_call = True
  )
  def show_moranRes_offcanvas(click):
    if click:
      return True

  @app.callback(
    Output('DATATABLE_moranRes_3D', 'data'),
    Output('DATATABLE_moranRes_3D', 'columns'),
    
    Input('BUTTON_calMoran_3D', 'n_clicks'),
    State('STORE_cellsIntersection_3D', 'data'),
    State('DROPDOWN_stage_3D', 'value'),
    prevent_initial_call=True,
    background = True,
    manager = background_callback_manager,
    running = [
      (Output('BUTTON_showMoran_3D', 'disabled'), True, False),
      (Output('BUTTON_calMoran_3D', 'children'), '< 1min', 'Compute'),
      (Output('BUTTON_calMoran_3D', 'loading'), True, False),
      (Output('OFFCANVAS_moranRes_3D', 'is_open'), False, True),
    ]
  )
  def cal_moranRes(click, cells, stage):
    
    adata = exp_data[stage]
    
    df = cal_moran_3D(adata[cells])
    df = df.reset_index(names='Feature')
    return (df.to_dict('records'),
            [
              {"name": i, "id": i, "deletable": False, 'type': 'numeric', 
                'format':Format(precision=4)} 
              for i in df.columns
            ]
          )

  # In[] run:

  tabs = html.Div(
    spatial_tabs,
  )

  app.layout = dbc.Container(
    [
      header,
      dbc.Row([
        dbc.Col([
          tabs,
        ], width=12)
      ],)
    ],
    fluid=True,
  )
  
  return app


