import plotly.graph_objects as go
import plotly.express as px


def bar_plot(
    xs: list[int], 
    ys:list[int],
    title: str = None, 
):

    fig = go.Figure([go.Bar(x=xs, y=ys)])
    fig.show()