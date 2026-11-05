#!/usr/bin/env python

import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.models import LogTicker, BasicTicker,ColumnDataSource,FixedTicker, Range1d, Slope, PreText, Select, CustomJS, Legend, RadioButtonGroup, HoverTool, LinearAxis, LogAxis, DataRange1d, Toggle
from bokeh.layouts import column, row


from bokeh.palettes import Category10, Category20, Colorblind, Turbo256
import sys
import math
import numpy as np

average=False
#avg_point = -15
avg_point = -5

font_size=20
marker_size=18
file_in = sys.argv[1]
unit = sys.argv[2] # MiB, GiB
#config = sys.argv[3]


metrics = [
    'Build time (sec)',
    f'Build peak memory ({unit})',
    f'Index size ({unit})',
    'Query time (sec)',
    f'Query peak memory ({unit})',
    'Build+Query time (sec)',
    f'Build+Query peak memory ({unit})',
    'TP', 'FP', 'FN', 'TN', 'avg_query'
]

file_out = ".".join(file_in.split(".")[:-1]) + ".html"
y_axis_type="log"  # "linear"

def make_markers(n_markers):
    # Markers in a "meaningful" order, considering their differences
    markers = [
        "circle",
        "plus",
        "triangle_pin",
        "square_dot",
        "square",
        "hex",
        "triangle",
        "diamond",
        "inverted_triangle",
        "square_pin",
        "circle_cross",
        "square_cross",
        "triangle_dot",
        "star_dot",
        "diamond_cross",
        "hex_dot",
        "circle_dot",
        "diamond_dot",
        "square_x",
        "triangle_pin",
        "circle_x",
        "circle_y",
        "star",
    ]
    if n_markers > len(markers):
        markers = markers * math.ceil(n_markers / len(markers))
    smarkers = markers[:n_markers]  # random.sample(markers, n_elements)
    return smarkers


def make_color_palette(n_colors):
    if n_colors == 1:
        # Only one raptor, keep green
        palette = [Colorblind[8][3]]
    elif n_colors <= 7:
        #skip double blue
        #palette = Colorblind[8][0:4] + Colorblind[8][6:]
        #             ('blue   ', 'orange ', 'yellow ', 'blugren', 'skyblue', 'vermill', 'redprpl', 'black  ') # key
        palette = ('#0072B2', '#E69F00', '#CC79A7', '#F0E442', '#009E73', '#56B4E9', '#000000') # reordered
    elif n_colors <= 10:
        palette = Category10[10]
    else:
        palette = Category10[10] * math.ceil(n_colors/10)
    return palette[:n_colors]


# Data
data = pd.read_csv(file_in, sep="\t")

data[metrics[1:]] = data[metrics[1:]].fillna(0)

# Replace mixed rebuild with full rebuild
data["Obs"] = data["Obs"].str.replace('Mixed Rebuild', 'Full Rebuild')


data[["Label", "Obs"]] = data[["Label", "Obs"]].fillna('')

data.sort_values(by="Label", inplace=True)
 
#data["tool"] = data["Label"].apply(lambda x: x.split(" ")[0])
data["tool"] = data["Label"]
data["tool"] = data["tool"].str.replace('raptor', 'raptor_dynamic')
data["tool"] = data["tool"].str.replace('raptor_dynamic_static', 'raptor_static')

data["LabelObs"] = data.Label + " " + data.Obs
data["ToolObs"] = data["tool"] + " " + data.Obs

labels = list(data["Label"].unique())
toolobs = sorted(list(data["ToolObs"].unique()))

markers = make_markers(len(toolobs))
colors = make_color_palette(len(labels))

x_start = "Run Number"
y_start = "Build time (sec)"

data["None"] = np.nan
data["x_data"] = data[x_start]
data["y_data"] = data[y_start]
data["marker_data"] = data["ToolObs"].apply(lambda x: markers[toolobs.index(x)])
data["color_data"] = data["Label"].apply(lambda x: colors[labels.index(x)])


data["alpha"] = data["Label"].apply(lambda x: 0.5 if not x.startswith("raptor ") else 1)
data["line_color"] = data["Label"].apply(lambda x: colors[labels.index(x)] if not x.startswith("raptor ") else "black")

# Put raptor as last to plot over the others
order = {"cobs":0, "fulgor": 1, "metagraph": 2, "raptor_static": 3}
data["sorter"] = data["tool"].apply(lambda x: order[x] if x in order else 999)

data.sort_values(by=["sorter"], inplace=True)

cds = ColumnDataSource(data)

p = figure(
    x_axis_label=x_start,
    y_axis_label=y_start,
    width=1600,
    height=800,
    y_axis_type=y_axis_type,
)

#p.y_range.start = 1
#p.y_range.start = 0.1

#p.yaxis.ticker = BasicTicker(desired_num_ticks=10)
#p.yaxis.ticker = LogTicker(base=10)

max_data = data[metrics].max().apply(lambda x: 10 ** np.ceil(np.log10(x)))
max_data = max_data + (max_data * 0.1)
min_data = data[metrics].min().apply(lambda x: 10 ** np.floor(np.log10(x)))
min_data = min_data - (min_data * 0.1)

max_cds = ColumnDataSource(pd.DataFrame(max_data).transpose())
min_cds = ColumnDataSource(pd.DataFrame(min_data).transpose())

p.y_range.start = min_cds.data[y_start][0]
p.y_range.end = max_cds.data[y_start][0]

# Add legend before plotting to be place outside the figure
p.add_layout(Legend(location="top"), "above")

p.scatter(
    x="x_data",
    y="y_data",
    source=cds,
    name="Label",
    legend_group="ToolObs", # LabelObs -> with params
    color="color_data",
    line_color="line_color",
    marker="marker_data",
    size=marker_size,
    fill_alpha="alpha"
)

mean_data = data.groupby("Label")[metrics].mean()

mean_data["x_data"] = avg_point
mean_data["y_data"] = mean_data[y_start]
mean_data["color_data"] = [colors[labels.index(x)] for x in mean_data.index.tolist()]
mean_data["alpha"] = [0.5 if not x.startswith("raptor ") else 1 for x in mean_data.index.tolist()]
mean_data["line_color"] = [colors[labels.index(x)] if not x.startswith("raptor ") else "black" for x in mean_data.index.tolist()]
mean_data["sorter"] = [1 if not x.startswith("raptor ") else 999 for x in mean_data.index.tolist()]
mean_data.sort_values(by=["sorter"], inplace=True)

mean_cds = ColumnDataSource(mean_data)

if average:
    p.scatter(
        x="x_data",
        y="y_data",
        source=mean_cds,
        name="Label",
        color="color_data",
        line_color="line_color",
        size=marker_size,
        marker="star",
        fill_alpha="alpha"
    )
    p.xaxis.major_label_overrides = {avg_point: "AVG"}
    #p.xaxis.ticker = FixedTicker(ticks=[avg_point] + [0, 20, 40, 60, 80])
    p.xaxis.ticker = FixedTicker(ticks=[avg_point] + [0, 20, 40, 60, 80, 100, 120, 140, 160, 180])
    #p.xaxis.ticker = FixedTicker(ticks=[avg_point] + [0, 100, 200, 300, 400])

# Log grid
p.ygrid.minor_grid_line_color = 'black'
p.ygrid.minor_grid_line_alpha = 0.1

tooltips = [
            ('Tool', '@tool'),
            ('Obs', '@Obs'),
            ('x', '@x_data'),
            ('y', '@y_data{0.00}'),
           ]
p.add_tools(HoverTool(tooltips=tooltips))

#p.title.text = config
#p.title.text_font_size = str(font_size + 3) + "pt"
p.legend.label_text_font_size = str(font_size - 2.5) + "pt"
p.legend.nrows = 4
p.legend.glyph_height = marker_size*3
p.legend.glyph_width = marker_size*3

#p.legend.location = "top_left"
# does not work: https://github.com/bokeh/bokeh/issues/10998
# https://discourse.bokeh.org/t/legend-field-multiple-data-sources-click-policy-is-click-policy-supposed-to-apply-to-entire-data-source-rather-than-just-one-legend-item/11342/2
#p.legend.click_policy="hide"
#p.legend.click_policy="mute"

p.xaxis.axis_label_text_font_size = str(font_size) + "pt"
p.yaxis.axis_label_text_font_size = str(font_size) + "pt"
p.yaxis.major_label_text_font_size = str(font_size) + "pt"
p.xaxis.major_label_text_font_size = str(font_size) + "pt"

toggle = Toggle(label="Hide legend", button_type="default")
toggle.js_on_event('button_click', CustomJS(args=dict(btn=toggle, leg=p.legend, p=p), code="""
    if(btn.active==true){
        leg[0].visible = false;
        p.height=665;
    }else{
        leg[0].visible = true;
        p.height=800;
    }
"""))

select_y = Select(title="y left:", value=y_start, options=metrics)
select_y.js_on_change("value", CustomJS(
    args=dict(
        cds=cds,
        mean_cds=mean_cds,
        min_cds=min_cds,
        max_cds=max_cds,
        yaxis=p.yaxis[0],
        y_range=p.y_range,
        select_y=select_y
    ),
    code="""

    yaxis.visible = true
    yaxis.axis_label = select_y.value

    if (min_cds.data[select_y.value][0]==0){
        y_range._initial_start = null
        y_range._initial_end = null
    }else{
        y_range._initial_start = min_cds.data[select_y.value][0]
        y_range._initial_end = max_cds.data[select_y.value][0]
    }

    cds.data["y_data"] = cds.data[select_y.value]
    cds.change.emit();

    mean_cds.data["y_data"] = mean_cds.data[select_y.value]
    mean_cds.change.emit();
    """,
))

output_file(file_out, mode="inline")

#with open(config, "r") as f:
#    pt = PreText(text="".join(f.readlines()), width=600, height=800)

save(column([p, select_y, toggle]))