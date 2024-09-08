#!/usr/bin/env python

import pandas as pd 
import argparse 
import plotly.graph_objects as go
import plotly.express as px
  
def prepare_votes(votes_path:str, n:int=10) -> pd.DataFrame:
    votes = pd.read_csv(votes_path, sep=" ")
    votes = votes.sort_values("cal_Freq", ascending=False)
    
    votes = votes.head(n)
    votes = votes.reset_index()
    
    others = round(100 - sum(votes["cal_Freq"]), 2)
    votes = pd.concat([votes, pd.DataFrame({11:{"index":"others", "Freq":0, "cal_Freq":others}}).T])
    
    return votes

def match_colors(votes, path) -> list:
    cls = pd.read_excel(path, names=["class", "color"], converters={"color":str})
    cls["class"] = cls["class"].str.replace(", ", "_").str.replace(" ", "_")
    cls["color"] = "#"+cls["color"]
    cls = cls[cls["class"].isin(votes["index"].tolist())]
    return cls["color"].tolist()
    
def plot_methylation(votes:pd.DataFrame, output:str, colors:list) -> None:
    colors = colors + ['#f5f4f2']
    fig = go.Figure()
    
    fig.add_trace(
        go.Pie(
            labels=votes["index"],
            values=votes["cal_Freq"],
            texttemplate="%{label}<br>"
                        "%{percent:.1%}",
            showlegend=False,
            pull=[0.2, 0, 0, 0, 0, 0],
            hoverinfo=None,
            marker_colors=colors
        )
    )
    
    fig.update_traces(
        # Size of the hole in the middle [0-1]
        hole=0.3
    )
    
    fig.update_layout(
        # Add annotations in the center of the donut.
        annotations=[
            dict(
                text='methylation class', 
                x=0.5, y=0.5, 
                font_size=14,
                showarrow=False
            )
        ],
        width=750, height=750
    )
    fig.write_image(output)
    
    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("votes", help="votes file generates by rapid-cnc² classification")
    parser.add_argument("-o", "--output", default="predictions.png")
    parser.add_argument("-n", "--number", default=10, type=int, help="number of predictions displayed")

    args = parser.parse_args()
    
    votes = prepare_votes(args.votes, args.number)
    colors = match_colors(votes, "rdata/class_colour.xlsx")
    plot_methylation(votes, args.output, colors)