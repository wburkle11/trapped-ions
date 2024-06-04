# **Helical Resonator**

Paul ion traps require high radio-frequency voltages applied to particular electrodes in order to properly provide electromagnetic potentials for trapping. Helical resonators allow impedence matching between a radio frequency source and an ion trap enabling high voltages while reducing noise injected into the experimental system. These properties make resonators of this kind important to not only ion trapping but also in a wide range of physical sciences. 

## Helical Resonator Functions

The Helical_Resonator_Functions.py file contain all functions required for numerical calculations of expected resonant frequency, and corresponding Qfactor from a helical resonator. This code replicates the mathematics outlined in James D Siverns paper in guiding resonator design; *"On the application of radio frequency voltages to ion traps via helical resonators"* [[Siverns 2012](https://arxiv.org/abs/1106.5013)]. 


## Helical Resonator Design

We can reduce the number of input paramters to a select few for charectorization of resonator performace for sake of simplicity. Following the dependece tree below, we can compute what our expected resonant frequency will be as a function of trap capacitance, as well as its corresponding Qfactor. Subsequent plots, are used in optimizing the resonator dimensions neccessary for maximizing this Qfactor. 

![dep_tree](https://github.com/wburkle11/trapped-ions/assets/92954143/b6d27875-7c73-4a44-8f46-ee09c98547ac)

Using this code I have been able to rediscover some of the guidelines highlighted in [[Siverns 2012](https://arxiv.org/abs/1106.5013)]. However, there are some instances where these "rules" regarding resonator dimensions are violated. For that reason, optimization plots are made for every neccesarry parameter. 

**Example Input:** 
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

Given these input parameters, "Helical_Resonator_Design.py" will output a collection of plots shown below: We can also collect values for resonant frequency and Q for a specific trap capacitance and resistance. Lastly, we return values for more optimal resonator dimensions for maximizing Q within a similar resonant frequency range.  

**Example Output:**

- Plot of Resonant Freq vs. Trap Capacitance

![res_fig_example](https://github.com/wburkle11/trapped-ions/assets/92954143/9b9a563a-b8e9-4759-b8ce-64a3366d94f5)

- Plot of Q vs. Trap Capacitance (depends on trap resistance)

![qfactor_example](https://github.com/wburkle11/trapped-ions/assets/92954143/bc4959d9-b5b8-476a-972f-6399a6f8d5ea)

Returned Values: The optimal parameters shown in this table are computed by finding the maximum Q in the functions below

![print](https://github.com/wburkle11/trapped-ions/assets/92954143/da39e499-7aec-44a8-9d1f-73ea70d5b286)

- Optimal Resonator Sheild Diameter (D)

![opt_D](https://github.com/wburkle11/trapped-ions/assets/92954143/6eee7178-6bb6-4659-bec8-5734f58331d5)

- Optimal Resonator Main Coil Diameter (d)

![opt_dd](https://github.com/wburkle11/trapped-ions/assets/92954143/56989dbb-18a4-4594-8706-0de8fda672dc)

- Optimal Resonator Main Coil Height (b)

![opt_b](https://github.com/wburkle11/trapped-ions/assets/92954143/c1223daf-a09b-4412-936b-278e479d7698)

- Optimal Resonator Main Coil Wire Diameter (d_0)

![opt_d_0](https://github.com/wburkle11/trapped-ions/assets/92954143/813b0b57-cbe6-49f5-b1a2-808fb204826a)

- Optimal Resonator Main Coil Winding Pitch (tau)

![opt_t](https://github.com/wburkle11/trapped-ions/assets/92954143/334a7c14-30d2-41b9-9699-93e29312ee02)

Using these outputs, by redefining the input parameters with the updated "Optimal" parameters you can begin to dial in on the changes to resonator dimensions that will provide a signifigant boost to the Q factor. **Updating these dimensions will cause a shift in resonant frequency.** However most of the time the resonant frequency curve will begin to "flatten out" taking similar values of omega_0 for a wide range of trap capacitances. 

**NOTE:** 

- Make sure to monitor optimal tau! All resonator functions assume the typical case where tau is roughly = 2*d_0 + 0.001: THERE ARE SOME CIRCUMSTANCES WHERE THIS IS NOT TRUE. Thus, you will need to update the definiton of tau throughout all of the functions
- The capacitance of the wire connections between resonator and trap takes on a defined value in my code of: C_w = 0.0001e-12 Farads : Depending on your expieremental apparatus, this capacitance can take on quite a large range of values, and should be updated accordingly.
- Plots of Q Factor "break down" or begin to diverge for trap resitances of < 0.5 Ohms, in comparison with plots made in [[Siverns 2012](https://arxiv.org/abs/1106.5013)].

[1] [Siverns, 2012: On the application of radio frequency voltages to ion traps via helical resonators](https://arxiv.org/abs/1106.5013) 







