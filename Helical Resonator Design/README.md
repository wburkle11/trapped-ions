# **Helical Resonator**

Paul ion traps require high radio-frequency voltages applied to particular electrodes in order to properly provide electromagnetic potentials for trapping. Helical resonators allow impedence matching between a radio frequency source and an ion trap enabling high voltages while reducing noise injected into the experimental system. These properties make resonators of this kind important to not only ion trapping but also in a wide range of physical sciences. 

## Helical Resonator Functions

The Helical_Resonator_Functions.py file contain all functions required for numerical calculations of expected resonant frequency, and corresponding Qfactor from a helical resonator. This code replicates the mathematics outlined in James D Siverns paper in guiding resonator design; *"On the application of radio frequency voltages to ion traps via helical resonators"* [[Siverns 2012](https://arxiv.org/abs/1106.5013)]. 


## Helical Resonator Design

We can reduce the number of input paramters to a select few for charectorization of resonator performace for sake of simplicity. Following the dependece tree below, we can compute what our expected resonant frequency will be as a function of trap capacitance, as well as its corresponding Qfactor. Subsequent plots, are used in optimizing the resonator dimensions neccessary for maximizing this Qfactor. 

![dep_tree](https://github.com/wburkle11/trapped-ions/assets/92954143/b6d27875-7c73-4a44-8f46-ee09c98547ac)

Using this code I have been able to rediscover some of the guidelines highlighted in [[Siverns 2012](https://arxiv.org/abs/1106.5013)]. However, there are some instances where these "rules" regarding resonator dimensions are violated. For that reason, optimization plots are made for every neccesarry parameter. 

**Example Input Code:** 
~~~~
# Helical Resonator Dimensions

# Sheild Height:
b = 0.178  
# Sheild Diameter
D = 0.0997  
# Main Coil Diameter
d = 0.055  
# Diameter of Main Coil Wire
d_0 = 0.003

#Trap Specific Parameters
R_t = 3 
C_t1 = 0.3e-10
~~~~
