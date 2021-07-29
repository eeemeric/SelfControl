# Class Definitions
NeuroPhysObject, an object that contains all the neurophysiology data aquired using the TDT or Plexon system.

# NeuroPhysObject includes the following properties and methods. 
## Properties for class NeuroPhysObject:
- HeaderInformation
- Digital, spikes and events
- Analog, LFPs

## Methods for class NeuroPhysObject:
- NeuroPhysObject, constructor. returns an object of class NeuroPhysObject
- getSpikes, returns list of spike names (nSpikes x 1, cell vector) and the spike timestamps relative to trial start (nSpik         
- getAnalog, returns list of LFP names and data
- getWFs, returns list of spike names and the spike waveforms           
- getEvents, returns all of the behavioral events and timestamps        
- getSDFRaster, returns rasters and histograms aligned on an event of interest. 
---
For all of these methods, you can use the help function to get the documentation (F1 line) for the object class. 
You can do this 2 different ways from the matlab command line
1. You can type "help classDefinition.method" (e.g., help objNeuroPhys.getSpikes)
2. Once the object is loaded into the workspace, you can type "help objectName.method" (e.g., help objNeuroPhys.getSpikes)