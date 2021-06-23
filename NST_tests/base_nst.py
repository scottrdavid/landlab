"""
Example Code of initializing nst using a variety of options
"""
#%% --------------------------Import Modules-----------------------------------
import warnings

import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import pandas as pd

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph
from landlab.io import read_shapefile

from landlab.plot import plot_network_and_parcels

#%% ---------------------Temp Functions----------------------------------------
#will move these to their own files eventually

def calculate_flow_depth(Q,nst):# input discharge (Q) and nst model
    slope=nst._channel_slope # channel slope (m/m)
    Dg=nst._d_mean_active # mean grainsize of active layer  (m)
    B=nst._width # channel width (m)
    g=9.81 # gravity (m/s2)
    #Solve for depth
    H=(Q*((2*Dg)**(1/6))/(B*np.sqrt(g*slope*8.1)))**(3/5)
    return H

#%% -------------------------Model Parameters----------------------------------
# model time parameters
time = [0.0] # model start time
dt = 60  # timestep (seconds)=1 min timestep
timesteps=10 # number of timesteps

# Sediment Transport Equation
transport_eqn="WilcockCrowe"

#%% --------------------Set up model intitial conditions-----------------------

# Initial conditions 
# input the type of geometry data you will load in 
inputGeometryType='num'  #num = manually inputting numbers; shp=shapefiles

# slope threshold (min slope allowed in the domain)
slopethreshold=1e-2 


# fill in data if using shapefiles
if inputGeometryType=='shp':
    DATA_DIR1 = r"D:\Box Sync\Wildfires\PLI_Watersheds\PLI_Watersheds\001\delete"
    DATA_DIR = r"D:\Box Sync\Wildfires\PLI_Watersheds\PLI_Watersheds\001"
    file = os.path.join(DATA_DIR, "a001_network.shp")    
    points_shapefile = os.path.join(DATA_DIR1, "a001_nodes_att.shp")
    grid = read_shapefile(
        file,
        points_shapefile=points_shapefile,
        node_fields=["usarea_km2", "elev_m"],
        link_fields=["usarea_km2", "Length_m"],
        link_field_conversion={"usarea_km2": "drainage_area",
            "Slope":"channel_slope", "Length_m":"reach_length"},
        node_field_conversion={
            "usarea_km2": "drainage_area",
            "elev_m": "topographic__elevation"},
        threshold=slopethreshold,
        )
    print('add in a channel width component to pull from shapefile')
    
    
# fill in data if manually inputing data
if inputGeometryType=='num':
    y_of_node = (0, 0, 0, 0)
    x_of_node = (0, 100, 200, 300)
    nodes_at_link = ((0,1), (1,2), (2,3))
    grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    
    grid.add_field("bedrock__elevation", [3., 2., 1., 0.], at="node") # m # m
    grid.add_field("reach_length", [100., 100, 100.], at="link")  # m
    
    B=15#  Channel width
    grid.at_link["channel_width"] = B * np.ones(grid.number_of_links) # m


grid.add_field("topographic__elevation", np.copy(grid.at_node["bedrock__elevation"]), at="node")

#%% --------------------Set up model boundary conditions----------------------
# Boundary Conditions data

# Input flow or depth data----options for loading in raw discharge data
flowinput='Discharge' # Options-->flow_depth;RI
Qunits='cms'# units of discharg input either cfs or cms

# input discharge data can be constant value or xlsx with hydrographs 
# example using a constant Q
#Q_data=1000# Simulate constant Q

# Example using a .xlsx file with hydrographs
Q_data=r"D:\Box Sync\NST\landlab\NST_tests\Qdata.xlsx"


# if specifying a discharge
if flowinput=='Discharge':
    # figure out if input is a sting (usgs discharge input) or a constant Q 
    TFQstr=isinstance(Q_data,str) 
    if TFQstr==True:# if input data is a hydrograph 
        #need to update the line below to reflect the num of cols, rows in your data (exclude headers)
        Q_data = pd.read_excel (Q_data,sheet_name='Sheet1',usecols="B:D",skiprows=0)
        Q_data=Q_data.values
        
        if Qunits=="cfs":
            Q_cms=Q_data*(0.3048**3)
        elif Qunits=="cms":
            Q_cms=Q_data.copy()
        else:
            print('Units not specified!! Assuming m3/s')
        Q=Q_cms.copy()
    
    if TFQstr==False: # if input data a constant discharge--steady state simulation
       if Qunits=='cfs':
          Q_data=Q_data*(0.3048**3)# convert cfs to cms  
          
       # make a Q for all time steps
       Q_cms=np.ones(grid.number_of_links)*Q_data
       #make Q for all links
       Q=np.matlib.repmat(Q_cms,timesteps,1)
       
    # make a fake field to initialize  nst
    grid.at_link["flow_depth"] = np.nan * np.ones(grid.number_of_links)
   
 
# if specifing a flow depth
if flowinput=='flow_depth':
    H=2 # Channel Depth (m)
    grid.at_link["flow_depth"] = H * np.ones(grid.number_of_links) # m

#specify recurrence interval for which to calculate flow
if flowinput=='RI':
    print('need code to build out this section')
    
    # Q_data=r"D:\Box Sync\NST\landlab\NST_tests\[SOMEDAILYFLOWFILE].xlsx"
    # return_period_years = 2 #years
    
    # Qdata = os.path.join(DATA_DIR, "USGS 10172200 RED BUTTE CREEK AT FORT DOUGLAS, NEAR SLC, UT.xlsx")
    # data = pd.read_excel (
    #         file,
    #         sheet_name='Mean Daily Flow',
    #         header = 3)
    
    # #drop other columns if needed
    # #data.dropna(subset=['Stop Number','Lithology Category', 'Size (cm)'])
    # flow_at_gage = data['Discharge, ft3/s (mean)'].values
      
    # #convert Q in ft3/s to m3/s
    # if Qunits=='cfs':
    #     flow_at_gage=flow_at_gage*(0.3048**3)

        
    
    # #HELP HELP -- need to work on handling dates
    # #date = data['Date'].values
    
    # #plot hydrograph
    # #add date as x-axis
    # #plt.plot(flow_at_gage_cms)
    # #add ylabel('Daily streamflow, m^3/s')
    
    # #HELP HELP -- identify and remove bad values
    # #Q(Q==0)=NaN;
    
   
    # # calculate recurrence interval flow
    # recurrence_interval_flow = calc_recurrence_interval_flow(
    #     flow_at_gage,return_period_years)



#---------------------Set up Sediment Parcels----------------------------------


# element_id is the link on which the parcel begins. 
# put a parcel at all links
element_id = np.repeat(np.arange(grid.number_of_links), 1)
element_id = np.expand_dims(element_id, axis=1)

# make just a parcel at the starting link
#element_id= np.array([[0]])

volume = 1*np.ones(np.shape(element_id))  # (m3)
active_layer = np.ones(np.shape(element_id)) # 1= active, 0 = inactive
density = 2650 * np.ones(np.size(element_id))  # (kg/m3)
abrasion_rate = 0 * np.ones(np.size(element_id)) # (mass loss /m)

# # Lognormal GSD
# medianD = 0.05 # m
# mu = np.log(medianD)
# sigma = np.log(2) #assume that D84 = sigma*D50
# np.random.seed(0)
#   # (m) the diameter of grains in each parcel
# D = np.random.lognormal(mu, sigma,np.shape(element_id)) 

D = 0.05 * np.ones(np.shape(element_id)) # m

time_arrival_in_link = np.random.rand(np.size(element_id), 1) 
#location_in_link = np.random.rand(np.size(element_id), 1) 
location_in_link = 0 * np.ones(np.shape(element_id))

lithology = ["quartzite"] * np.size(element_id)

variables = {
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], density),
    "lithology": (["item_id"], lithology),
    "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
    "active_layer": (["item_id", "time"], active_layer),
    "location_in_link": (["item_id", "time"], location_in_link),
    "D": (["item_id", "time"], D),
    "volume": (["item_id", "time"], volume)}

items = {"grid_element": "link", "element_id": element_id}

_OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1

parcels = DataRecord(
    grid,
    items=items,
    time=[0.0],
    data_vars=variables,
    dummy_elements={"link":[NetworkSedimentTransporter.OUT_OF_NETWORK]},)


#---------------------Transport Sediment---------------------------------------
fd = FlowDirectorSteepest(grid, "topographic__elevation")
fd.run_one_step()

nst = NetworkSedimentTransporter(    
    grid,
    parcels,
    fd,
    bed_porosity=0.3,
    g=9.81,
    fluid_density=1000,
    transport_method="WilcockCrowe",)


if flowinput=='Discharge' or 'RI':
    for t in range(0, (timesteps * dt), dt):
        print(t)
        idx=int(t/dt)
        print (idx)
        H=calculate_flow_depth(Q[idx,:], nst)
        nst._grid.at_link.dataset.flow_depth.data=H
        nst.run_one_step(dt)  
        print("Model time: ", t/(60*60), "hrs passed")

if flowinput=='flow_depth':
    for t in range(0, (timesteps * dt), dt): 
        nst.run_one_step(dt)  
        print("Model time: ", t/(60*60), "hrs passed")
              