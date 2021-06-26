from mesa.visualization.ModularVisualization import ModularServer
from mesa.visualization.modules import CanvasGrid, ChartModule
from mesa.visualization.UserParam import UserSettableParameter

from TB.agents import Env, T, RestMP, InfectMP, ChronInfectMP, ActivatedMP, Source, Necrosis
from TB.model import TB

def TB_portrayal(agent):
    if agent is None:
        return
    
    portrayal = {}

    if (type(agent) is RestMP):
        portrayal["Color"] = ["#84e184", "#adebad", "#d6f5d6"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 1
        portrayal["w"] = 1
        portrayal["h"] = 1
    
    elif (type(agent) is ActivatedMP):
        portrayal["Color"] = ["#0000FF", "#0000FF", "#0000FF"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 1
        portrayal["w"] = 1
        portrayal["h"] = 1

    elif type(agent) is InfectMP:
        portrayal["Color"] = ["#FFA500", "#FFA500", "#FFA500"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 1
        portrayal["w"] = 1
        portrayal["h"] = 1

    elif type(agent) is ChronInfectMP:
        portrayal["Color"] = ["#FF0000", "#FFF000", "#FF0000"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 1
        portrayal["w"] = 1
        portrayal["h"] = 1  

    elif type(agent) is T:
        portrayal["Color"] = ["#FFC0CB", "#FFC0CB", "#FFC0CB"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 2
        portrayal["w"] = 1
        portrayal["h"] = 1   
    
    elif type(agent) is Source:
        portrayal["Color"] = ["#000000", "#000000", "#000000"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 0
        portrayal["w"] = 1
        portrayal["h"] = 1
    
    elif type(agent) is Necrosis:
        portrayal["Color"] = ["#CC7722", "#CC7722", "#CC7722"]
        portrayal["Shape"] = "rect"
        portrayal["Filled"] = 'true'
        portrayal["Layer"] = 3
        portrayal["w"] = 1
        portrayal["h"] = 1
    
    return portrayal

canvas_element = CanvasGrid(TB_portrayal, 100, 100, 500, 500)
server = ModularServer(
    TB, [canvas_element], "TB")
server.port = 8521
