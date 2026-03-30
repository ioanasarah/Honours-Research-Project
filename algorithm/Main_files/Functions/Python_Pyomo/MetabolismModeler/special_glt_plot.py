import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import os
import gurobi_logtools as glt_gurb
from io import StringIO
import json

from Functions.Python_Pyomo.MetabolismModeler import Analysis as an









def create_layout():
    pass




def run_plotly_dash(path, df_name = "*"):
    ### add task and sample selection
    ### > listselection (2 sets of dropdowns  with and without inclusion)
    ### custom seleciton of tasks
    ### set default values and make graphs work properly all the time
    ### add option to add more logs
    ### add option for trendline in case of scatter
    ### auto export



    os.chdir(path)
    summary, timelines = glt_gurb.get_dataframe([f"{df_name}.log"], timelines=True)
    df = summary
    df = timelines["nodelog"]

    app = dash.Dash()
    # add state var to store df

    y2_columns_with_None = df.columns.to_list()
    y2_columns_with_None.insert(0, "None")
    app.layout = html.Div([
        html.Label('Data'),
        # dcc.Store(id='df', data=df.to_json()),
        # dcc.Store(id='df', data=df.to_dict()),
        dcc.Dropdown(
            id='data',
            options=[{'label': "Summary", 'value': "summary"},
                     {'label': "Timelines", 'value': "timelines"}],
            value='timelines',
            maxHeight = 310,
            optionHeight= 25
        ),

        html.Label('X-axis'),
        dcc.Dropdown(
            id='x-axis',
            options=[{'label': html.Span(i, style = {'font-size': 14 }),
                      'value': i} for i in df.columns],
            value= "Time" if "Time" in df.columns else df.columns[0],
            maxHeight = 310,
            optionHeight= 25
        ),
        html.Label('Y-axis'),
        dcc.Dropdown(
            id='y-axis',
            options=[{'label': html.Span(i, style = {'font-size': 14 }),
                      'value': i} for i in df.columns],
            value= "Gap" if "Gap" in df.columns else df.columns[1],
            maxHeight = 310,
            optionHeight= 25
        ),
        html.Label('Second Y-axis'),
        dcc.Dropdown(
            id='y-axis2',
            options=[{'label': html.Span(i, style = {'font-size': 14 }),
                      'value': i} for i in y2_columns_with_None],
            value= "None",
            maxHeight = 310,
            optionHeight= 25
        ),

        html.Label('Color'),
        dcc.Dropdown(
            id='color',
            options=[{'label': html.Span(i, style = {'font-size': 14 }),
                      'value': i} for i in df.columns],
            value= "LogFilePath" if "LogFilePath" in df.columns else df.columns[2],
            maxHeight = 310,
            optionHeight= 25
        ),
        html.Label('Chart Type'),
        dcc.Dropdown(
            id='chart-type',
            options=[{'label': html.Span(i, style = {'font-size': 14 }),
                      'value': i} for i in ['box', 'bar', 'scatter', 'line']],
            value='line',
            maxHeight = 310,
            optionHeight= 25
        ),
        dbc.Checkbox(
            id='log-x',
            label='Log X',
            value=False
        ),
        dbc.Checkbox(
            id='log-y',
            label='Log Y',
            value=False
        ),
        dcc.Graph(id='graph')
    ])

    # update all dropdowns based on selected data and update df for graph
    @app.callback(
        Output('x-axis', 'options'),
        Output('x-axis', 'value'),
        Output('y-axis', 'options'),
        Output('y-axis', 'value'),
        Output('color', 'options'),
        Output('color', 'value'),
        Output('y-axis2', 'options'),
        Output('y-axis2', 'value'),
        Output('graph', 'figure'),
        # Output('df', 'data'),

        Input('data', 'value'),
        State('chart-type', 'value'),
        State('log-x', 'value'),
        State('log-y', 'value')
    )
    def update_dropdowns(data, chart_type, log_x, log_y):
        if data == "summary":
            df = summary
        else:
            df = timelines["nodelog"]
        fig = None
        if chart_type == 'box':
            fig = px.box(df, x=df.columns[0], y=df.columns[1], color=df.columns[2], log_x=log_x, log_y=log_y)
        elif chart_type == 'bar':
            fig = px.bar(df, x=df.columns[0], y=df.columns[1], color=df.columns[2], log_x=log_x, log_y=log_y)
        elif chart_type == 'scatter':
            fig = px.scatter(df, x=df.columns[0], y=df.columns[1], color=df.columns[2], log_x=log_x, log_y=log_y)
        elif chart_type == 'line':
            fig = px.line(df, x=df.columns[0], y=df.columns[1], color=df.columns[2], log_x=log_x, log_y=log_y)


        y2_columns_with_None = df.columns.to_list()
        y2_columns_with_None.insert(0, "None")

        return ([{'label': html.Span(i, style = {'font-size': 14 }),
                 'value': i} for i in df.columns],
                (("Time" if "Time" in df.columns else df.columns[0]) if data == "timelines" else "LogFilePath"),
                [{'label': html.Span(i, style = {'font-size': 14 }),
                 'value': i} for i in df.columns],
                (("Incumbent" if "Incumbent" in df.columns else df.columns[1]) if data == "timelines" else "Runtime"),
                [{'label': html.Span(i, style = {'font-size': 14 }),
                 'value': i} for i in df.columns],
                "LogFilePath" if "LogFilePath" in df.columns else df.columns[2],
                [{'label': html.Span(i, style = {'font-size': 14 }),
                    'value': i} for i in y2_columns_with_None],
                "None",
                fig if fig is not None else dash.no_update,
                # df.to_dict()
                # df
               )


    @app.callback(
        Output('graph', 'figure', allow_duplicate=True),
        Input('x-axis', 'value'),
        Input('y-axis', 'value'),
        Input('y-axis2', 'value'),
        Input('color', 'value'),
        Input('chart-type', 'value'),
        Input('log-x', 'value'),
        Input('log-y', 'value'),
        # State('df', 'data'),
        State('data', 'value'),
        prevent_initial_call=True
    )
    def update_graph(x, y, y2, color, chart_type, log_x, log_y, data):
        if data == "summary":
            df = summary
        else:
            df = timelines["nodelog"]
        # print(y)
        # print(df[y])

        # fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig = go.Figure()

        if chart_type == 'box':
            fig_trace_function = go.Box
        elif chart_type == 'bar':
            fig_trace_function = go.Bar
        elif chart_type == 'scatter':
            fig_trace_function = go.Scatter
        elif chart_type == 'line':
            fig_trace_function = go.Scatter

        x_data = df[x].values
        # replace 0s with 0.0001 to avoid log(0) error
        if data == "summary" and x == "Time":
            x_data = pd.Series(x_data).replace(0, 1).values
        else:
            x_data = pd.Series(x_data).replace(0, 0.0001).values
        y_data = df[y].values
        y_data = pd.Series(y_data).replace(0, 0.0001).values
        colors = px.colors.qualitative.Plotly
        double_df_color_unique = [str(i) + "_y1" for i in df[color].unique()]
        double_df_color_unique_y2 = [str(i) + "_y2" for i in df[color].unique()]
        double_df_color_unique.extend(double_df_color_unique_y2)
        print(len(double_df_color_unique))

        color_dict = {k: v for k, v in zip(double_df_color_unique, colors)}
        traces = []
        for i, c in color_dict.items():
            if i.endswith("_y1"):
                traces.append(fig_trace_function(x=x_data[df[color] == i.split("_y1")[0]],
                                                 y=y_data[df[color] == i.split("_y1")[0]],
                                                 name=i + "_y1",
                                                 marker=dict(color=c),
                                                 # log_x=log_x,
                                                 # log_y=log_y
                                                 ))
        if y2 is not None and y2 != "None":
            # color_dict = {k: v for k, v in zip(double_df_color_unique, colors)}
            y2_data = df[y2].values
            y2_data = pd.Series(y2_data).replace(0, 0.0001).values
            color_dict = {k: v for k, v in zip(double_df_color_unique, colors)}
            for i, c in color_dict.items():
                if i.endswith("y2"):
                    traces.append(fig_trace_function(x=x_data[df[color] == i.split("_y2")[0]],
                                                     y=y2_data[df[color] == i.split("_y2")[0]],
                                                     name=i + "_y2",
                                                     marker=dict(color=c),
                                                     # log_x=log_x,
                                                     # log_y=log_y
                                                     ))
        print(len(traces))
        to_add = f" and {y2}" if y2 is not None and y2 != "None" else ""
        layout = go.Layout(
            yaxis=dict(title=f"{y}{to_add}"),
            xaxis=dict(title=x),
        )
        fig = go.Figure(data=traces, layout=layout)
        if x == "LogFilePath":
            fig.update_xaxes(type="category")

            def split_string_into_chunks(s, chunk_size):
                return [s[i:i + chunk_size] for i in range(0, len(s), chunk_size)]

            fig.update_xaxes(tickmode='array',
                             tickangle=45,
                             tickvals=df[x].unique(),
                             ticktext=['<br>'.join(split_string_into_chunks(val, 20)) for val in df[x].unique()])
        elif log_x:
            fig.update_xaxes(type="log")
        else:
            fig.update_xaxes(type="linear")
        if log_y:
            fig.update_yaxes(type="log")
        else:
            fig.update_yaxes(type="linear")


        return fig

    app.run_server(debug=True,
                   dev_tools_hot_reload=False,
                   )

def run_dash_app(layout, options):
    app = dash.Dash()
    app.layout = layout
    app.run_server(debug=False)

if __name__ == '__main__':
    # main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_NT_LT_INIT_LINEAR_MIN_Tarray_3samples_8tasks"
    main_folder = r"E:\Git\Metabolic_Task_Score\Data\Pyomo_python_files\Jelle\test_runs\consensus_no_threshold_fractional_options_test_enum"
    an.check_all_json_models_in_dir_for_feasibility(main_folder)
    # if not os.getcwd().endswith("results"):
    #     os.chdir(os.path.join(os.getcwd(), "results"))
    run_plotly_dash(main_folder)

