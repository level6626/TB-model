"""
TB Model
"""


from mesa import Model
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
import TB.para as para

from TB.agents import Env, T, RestMP, InfectMP, ChronInfectMP, ActivatedMP, Source
from TB.schedule import RandomActivationByBreed

class TB(Model):

    height = 100
    width = 100

    verbose = True # Print-monitoring

    description = (
        "TB model"
    )

    def __init__(
        self,
        height = 100,
        width = 100
    ):
        super().__init__()
        self.height = height
        self.width = width
        self.env = Env(self.next_id(), self)
        self.schedule = RandomActivationByBreed(self)
        self.grid = MultiGrid(self.height, self.width, torus = False)
        self.datacollector = DataCollector(
            {
                "RestMP": lambda m: m.schedule.get_breed_count(RestMP),
                "InfectMP": lambda m: m.schedule.get_breed_count(InfectMP),
                "ChonInfectMP": lambda m: m.schedule.get_breed_count(ChronInfectMP),
                "ActivatedMP": lambda m: m.schedule.get_breed_count(ActivatedMP),
                "T": lambda m: m.schedule.get_breed_count(T),
            }
        )

        #Create Environment & Extracellular Bacteria
        for (x, y) in [(49, 49), (49, 50), (50, 49), (50, 50)]:
            self.env.BE[x, y] = para.BE_init / 4
        self.schedule.add(self.env)
        
        #Create resting macrophage
        for i in range(para.M_init):
            x = self.random.randrange(self.width)
            y = self.random.randrange(self.height)
            MP = RestMP(self.next_id(), (x, y), self, True)
            self.grid.place_agent(MP, (x, y))
            self.schedule.add(MP)
        MP_to_be_infected = RestMP(self.next_id(), (49, 49), self, True)
        self.grid.place_agent(MP_to_be_infected, (49, 49))
        self.schedule.add(MP_to_be_infected)  

        #Create blood vessel (Source)
        for pos in [(25, 25), (25, 75) ,(75, 25), (75, 75)]:
            src = Source(self.next_id(), pos, self)
            self.grid.place_agent(src, pos)
            self.schedule.add(src)
        
        
        self.running = True
        self.datacollector.collect(self)
    
    def step(self):
        self.schedule.step()

        # collect data for testing
        self.datacollector.collect(self)
        if self.verbose:
            print(
                [
                    self.schedule.time,
                    self.schedule.get_breed_count(RestMP),
                    self.schedule.get_breed_count(InfectMP),
                    self.schedule.get_breed_count(ChronInfectMP),
                    self.schedule.get_breed_count(ActivatedMP),
                    self.schedule.get_breed_count(T),
                ]
            )
    
    def run_model(self, step_count = para.t_total):
        if self.verbose:
            print("Initial number RestMPs: ", self.schedule.get_breed_count(RestMP))

        for i in range(step_count):
            self.step()

        if self.verbose:
            print("")
            print("Final number restMPs: ", self.schedule.get_breed_count(RestMP))
