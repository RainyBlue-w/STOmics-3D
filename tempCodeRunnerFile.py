## In[] env

from math import isnan
import math
from functools import reduce
from dash import Dash, dcc, html, dash_table, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash import ALL, MATCH, ALLSMALLER
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate
import dash_daq as daq
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components as fac
import feffery_utils_components as fuc
import dash_bootstrap_components as dbc
from dash_extensions.enrich import Output, Input, html, callback, DashProxy, LogTransform, DashLogger, Serverside, ServersideOutputTransform
from dash_extensions.enrich import MultiplexerTransform

import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
from plotnine import *
import plotnine.options

from PIL import Image
import scanpy as sc
import os
import pandas as pd
import dask.dataframe as dd
import numpy as np
# import loompy as lp
import h5py
import json
import time

import squidpy as sq

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import re
import seaborn as sns
from concurrent import futures
from typing import List, Dict, Tuple
import diskcache
background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/wuc/dash_data/spatial/cache_test"))


## In[] data