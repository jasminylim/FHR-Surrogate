# FHR System Surrogate Model

## Models:  
* ```systemSurrogate``` : gFHR System Surrogate with detailed system state for reference governor.  
* ```systemSurrogatePump``` : gFHR System Surrogate with detailed system state for reference governor and pump degradation surrogate model.  

## Surrogate Model Assumptions:
1. Constant time-step of 5 seconds.
2. Two time-delay system state input (```xk```). One time-delay pump degradation state input (```xkPump```).
   * pump degradation model only depends on the previous loss coeffcient values in the ```xkPump``` vector.
3. System state values do not account for degradation. Degradation can be detected through the pump degradation state values (```resPump```).
4. Full power (100%) is $280\times10^6$ [W]. 
5. Maximum power change rate allowed is 20% of full power per minute. If the user input violates this constraint, the surrogate model will modify the input to the maximum change allowed. See the ```systemSurrogate.change_input()``` function for details on how the target power will be modified.
  
## How to run models
1. Ensure the following text files are accessible (found in ```../models```):
   1. ```varmax-surrogate-RGPump1.txt```
   2. ```varmax-surrogate-RGPump2.txt```
   3. ```varmax-surrogate-RGPump3.txt```
   4. ```varmax-surrogate-RGPump4.txt```
   5. ```varmax-surrogate-RGPump5.txt```

2. Ensure the following surrogate model files are accessible: 
   1. ```systemSurrogate.py``` : FHR System Surrogate Model  
   2. ```PP_SurrogateModel``` : Pump Degradation Surrogate Model

   
3. Refer to the test cases. This contains the information about each model, including their states and usage details:
   * ```test_LF.ipynb``` : Load Follows Case (w/ SAM data).  
   * ```test_LF_wPump.ipynb``` : Load Follows Case (w/ SAM data) using system surrogate model with pump degradation surrogate model integrated.

## System Surrogate Model Prediction
The ```systemSurrogate.predict()``` function runs the surrogate model for the user defined target power input. 

#### Function Call:

```
res = systemSurrogate.predict( 
         xk : float, 
         uP : float)
```

#### Inputs:
* ```xk``` [array (2,42)]: Current system state values, includes a 2 step time-delay.
* ```uP``` [array (n,2)]: Targer power distribution for ```n```-step time prediction.

#### Returns:
* ```res``` [array (n, 42)]: n-step system state prediction.
  
### ```xk``` : System States
The system surrogate state vector $\mathbf{x} \in \mathbb{R}^{42}$ is 42-dimensional and contains the states listed in the table below, where "Index" is reffering to the index of the array:

|Index|Variable|Units|
|----|------|----|
|0|Time|s|
|1|$N_I$ : Relative I-135 Concentration|$\Delta$ k/k|
|2|$N_{Xe}$ : Relative Xe-135 Concentration |$\Delta$ k/k|
|3|$T_{c,in}$ : Core Inlet Temperature|K|
|4|$T_{c,out}$ : Core Outlet Temperature|K|
|5|$P_{c,in}$ : Core Inlet Pressure|Pa|
|6|$P_{c,out}$ : Core Outlet Pressure|Pa|
|7|$T_{ihx,p,in}$ : IHX Primary Inlet Temperature|K|
|8|$T_{ihx,p,out}$ : IHX Primary Outlet Temperature|K|
|9|$P_{ihx,p,in}$ : IHX Primary Inlet Pressure|Pa|
|10|$P_{ihx,p,out}$ : IHX Primary Outlet Pressure|Pa|
|11|$T_{ihx,s,in}$ : IHX Secondary Inlet Temperature|K|
|12|$T_{ihx,s,out}$ : IHX Secondary Outlet Temperature|K|
|13|$P_{ihx,s,in}$ : IHX Secondary Inlet Pressure|Pa|
|14|$P_{ihx,s,out}$ : IHX Secondary Outlet Pressure|Pa|
|15|$T_{sg,in}$ : SG Inlet Temperature|K|
|16|$T_{sg,out}$ : SG Outlet Temperature|K|
|17|$P_{sg,in}$ : SG Inlet Pressure|Pa|
|18|$P_{sg,out}$ : SG Outlet Pressure|Pa|
|19|$\dot{m}_c$ : Core Flow |kg/s|
|20|$\dot{m}_s$ : Solar-salt Mass Flow Rate|kg/s|
|21|$\dot{m}_{sg}$ : SG Mass Flow Rate|kg/s|
|22|$\dot{Q}_{RX}$ : Reactor Core Power|W|
|23|$\dot{Q}_{HX}$ : IHX Heat Extraction Rate|W|
|24|$\dot{Q}_{SG}$ : SG Heat Extraction Rate|W|
|25|$C_1$ : Delayed Neutron Precursor Concentration - group 1|-|
|26|$C_2$ : Delayed Neutron Precursor Concentration - group 2|-|
|27|$C_3$ : Delayed Neutron Precursor Concentration - group 3|-|
|28|$C_4$ : Delayed Neutron Precursor Concentration - group 4|-|
|29|$C_5$ : Delayed Neutron Precursor Concentration - group 5|-|
|30|$C_6$ : Delayed Neutron Precursor Concentration - group 6|-|
|31|$\rho_{m}$ : Moderator Reactivity Feedback|$\Delta$ k/k|
|32|$\rho_c$ : Coolant Reactivity Feedback|$\Delta$ k/k|
|33|$\rho_{f}$ : Fuel Reactivity Feedback|$\Delta$ k/k|
|34|$\rho_{ext}$ : External Reactivity Insertion (control rod)|$\Delta$ k/k|
|35|$z_{cr}$ : Control Rod Position |m|
|36|$\dot{m}_{P,p}$ : Primary Pump Mass Flow Rate|kg/s|
|37|$\dot{m}_{P,s}$ : Secondary Pump Mass Flow Rate|kg/s|
|38|$n_p$ : Primary Pump Speed|RPM|
|39|$n_s$ : Secondary Pump Speed|RPM|
|40|$\Delta P_p$ : Primary Pump Head|Pa|
|41|$\Delta P_s$ : Secondary Pump Head|Pa|

The names of the system state values can be accessed with ```systemSurrogate.states``` and the state labels/units can be accessed with ```systemSurrogate.state_labels```. 

### ```uP``` : System Inputs

The system surrogate state vector $\mathbf{u} \in \mathbb{R}^{2}$ is 2-dimensional and contains the inputs listed in the table below, where "Index" is reffering to the index of the array:

|Index|Variable|Units|
|----|------|----|
|0|Time|s|
|1|$\dot{Q}_{RX,T}$ : Target Reactor Core Power|W|

The names of the input values can be accessed with ```systemSurrogate.inputs``` and the input labels/units can be accessed with ```systemSurrogate.input_labels```.
  
## System Surrogate w/ Pump Degradation Model Prediction
The ```systemSurrogatePump.predict()``` function runs the surrogate model for the user defined target power input. Inlcudes pump degradation surrogate model predictions.

#### Function Call:

```
res, pumpRes = systemSurrogatePump.predict( 
                xk : float, 
                xkPump : float,
                uP : float, 
                newPump : dict = None, 
                FullyDegradeTime : dict = {"Primary":10*365*24*3600,"Secondary":10*365*24*3600},
                FullyDegradeHeadLoss : dict = {"Primary": 0.5, "Secondary": 0.5}, 
                DegradeUncertainty : dict = {"Primary": 0.001, "Secondary": 0.001}, 
                DegradeInvScale : dict = {"Primary": 0., "Secondary": 0.}, 
                FlowRateLossCoeff : dict = {"Primary": 0., "Secondary": 0.},
                FullDegradeFlowRateInv : dict = {"Primary": 0., "Secondary": 0.},
                deltaP : int = 1,
                ShowNoDegradePower : bool = False)
```

#### Inputs:
* ```xk``` [array (2,42)]: Current system state values, includes a 2 step time-delay.
* ```xkPump``` [array (9,)] : Current pump degradation state values.
* ```uP``` [array (n,2)]: Targer power distribution for ```n/DeltaP```-step time prediction.
* ```newPump``` [dict["Primary","Secondary"]=(n,)]: Indicates at each time-step if the pump is changed. New pump denoted by ```1```, and no new pump with ```0```. Defined for each pump with dictionary keys ```Primary, Secondary```. The default value ```None``` indicates no pump change during prediction period. A new pump will reset the degradation values.
* ```FullyDegradeTime``` [dict["Primary","Secondary"]] : Time at which full head loss will occur. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```FulleDegradeHeadLoss``` [dict["Primary","Secondary"]]: Percent head loss when pump is fully degraded. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```DegradeUncertainty``` [dict["Primary","Secondary"]]: pump degradation uncertainity. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```DegradeInvScale``` : used to scale the dKdt_nominal. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```FlowRateLossCoeff``` [dict["Primary","Secondary"]]: Coefficient for flow rate degradation. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```FullDegradeFlowRateInv``` : Inverse of fully degrade flow rate change. Defined for each pump with dictionary keys ```Primary, Secondary```.
* ```deltaP``` [int]: Sets compute step for pump degradation surrogate model. Default value if 1.
* ```ShowNoDegradePower``` [bool]: Indicates if pump power without degradation is computed. Default is True. If set to False, the pump power without degradation will $\dot{Q}_{P,p} = \dot{Q}_{P,s} = -1$ and the health pump indices $\eta_p, \eta_s$ will no be computer correctly.


For more details about the pump degradation surrogate model, please refer to the [```PP_SurrogateModel``` file](PP_SurrogateModel.py).

#### Returns:
* ```res``` [array (n, 42)]: n-step system state prediction.
* ```pumpRes``` [array (n/deltaP, 9)] : n/deltaP pump degradation state prediction. If ```ShowNoDegradePower == False```, the non-degradation pump values will be -1 and pump health score will be inaccurate.
  
### ```xk``` : System States
The system surrogate state vector $\mathbf{x} \in \mathbb{R}^{42}$ is 42-dimensional and contains the states listed in the table below, where "Index" is reffering to the index of the array:

|Index|Variable|Units|
|----|------|----|
|0|Time|s|
|1|$N_I$ : Relative I-135 Concentration|$\Delta$ k/k|
|2|$N_{Xe}$ : Relative Xe-135 Concentration |$\Delta$ k/k|
|3|$T_{c,in}$ : Core Inlet Temperature|K|
|4|$T_{c,out}$ : Core Outlet Temperature|K|
|5|$P_{c,in}$ : Core Inlet Pressure|Pa|
|6|$P_{c,out}$ : Core Outlet Pressure|Pa|
|7|$T_{ihx,p,in}$ : IHX Primary Inlet Temperature|K|
|8|$T_{ihx,p,out}$ : IHX Primary Outlet Temperature|K|
|9|$P_{ihx,p,in}$ : IHX Primary Inlet Pressure|Pa|
|10|$P_{ihx,p,out}$ : IHX Primary Outlet Pressure|Pa|
|11|$T_{ihx,s,in}$ : IHX Secondary Inlet Temperature|K|
|12|$T_{ihx,s,out}$ : IHX Secondary Outlet Temperature|K|
|13|$P_{ihx,s,in}$ : IHX Secondary Inlet Pressure|Pa|
|14|$P_{ihx,s,out}$ : IHX Secondary Outlet Pressure|Pa|
|15|$T_{sg,in}$ : SG Inlet Temperature|K|
|16|$T_{sg,out}$ : SG Outlet Temperature|K|
|17|$P_{sg,in}$ : SG Inlet Pressure|Pa|
|18|$P_{sg,out}$ : SG Outlet Pressure|Pa|
|19|$\dot{m}_c$ : Core Flow |kg/s|
|20|$\dot{m}_s$ : Solar-salt Mass Flow Rate|kg/s|
|21|$\dot{m}_{sg}$ : SG Mass Flow Rate|kg/s|
|22|$\dot{Q}_{RX}$ : Reactor Core Power|W|
|23|$\dot{Q}_{HX}$ : IHX Heat Extraction Rate|W|
|24|$\dot{Q}_{SG}$ : SG Heat Extraction Rate|W|
|25|$C_1$ : Delayed Neutron Precursor Concentration - group 1|-|
|26|$C_2$ : Delayed Neutron Precursor Concentration - group 2|-|
|27|$C_3$ : Delayed Neutron Precursor Concentration - group 3|-|
|28|$C_4$ : Delayed Neutron Precursor Concentration - group 4|-|
|29|$C_5$ : Delayed Neutron Precursor Concentration - group 5|-|
|30|$C_6$ : Delayed Neutron Precursor Concentration - group 6|-|
|31|$\rho_{m}$ : Moderator Reactivity Feedback|$\Delta$ k/k|
|32|$\rho_c$ : Coolant Reactivity Feedback|$\Delta$ k/k|
|33|$\rho_{f}$ : Fuel Reactivity Feedback|$\Delta$ k/k|
|34|$\rho_{ext}$ : External Reactivity Insertion (control rod)|$\Delta$ k/k|
|35|$z_{cr}$ : Control Rod Position |m|
|36|$\dot{m}_{P,p}$ : Primary Pump Mass Flow Rate|kg/s|
|37|$\dot{m}_{P,s}$ : Secondary Pump Mass Flow Rate|kg/s|
|38|$n_p$ : Primary Pump Speed|RPM|
|39|$n_s$ : Secondary Pump Speed|RPM|
|40|$\Delta P_p$ : Primary Pump Head|Pa|
|41|$\Delta P_s$ : Secondary Pump Head|Pa|

The names of the system state values can be accessed with ```systemSurrogatePump.states``` and the state labels/units can be accessed with ```systemSurrogatePump.state_labels```. 

### ```xkPump``` : Pump Degradation States
The pump degradation state vector $\mathbf{x}_k \in \mathbb{R}^{9}$ is 9-dimensional and contains the states listed in the table below, where "Index" is reffering to the index of the array:

|Index|Variable|Units|
----|----|----|
|0|Time|s|
|1|$\dot{Q}_{P,D,p}$ : Primary Required Pump Power with Degradation|W|
|2|$\dot{Q}_{P,p}$ : Primary Required Pump Power|W|
|3|$K_p$ : Primary Pump Loss Coeffcient|-|
|4|$\eta_p$ : Primary Pump Health Index|-|
|5|$\dot{Q}_{P,D,s}$ : Secondary Required Pump Power with Degradation|W|
|6|$\dot{Q}_{P,s}$ : Secondary Required Pump Power|W|
|7|$K_s$ : Secondary Pump Loss Coefficient|-|
|8|$\eta_s$ : Secondary Pump Health Index|-|

The names of the pump degradation model state values can be accessed with ```systemSurrogatePump.pump_states``` and the pump degradation state labels/units can be accessed with ```systemSurrogatePump.pump_state_labels```.

### ```uP``` : System Inputs

The system surrogate state vector $\mathbf{u} \in \mathbb{R}^{2}$ is a 2-dimensional and contains the inputs listed in the table below, where "Index" is reffering to the index of the array:

|Index|Variable|Units|
|----|------|----|
|0|Time|s|
|1|$\dot{Q}_{RX,T}$ : Target Reactor Core Power|W|


The names of the input values can be accessed with ```systemSurrogatePump.inputs``` and the input labels/units can be accessed with ```systemSurrogatePump.input_labels```.


## Datasets
This repo uses two datasets:
1. ```../data/gFHR-LF.csv``` : Load Follows Case (100% full power -> 60% full power -> 100% full power)
2. ```../data/gFHR-LF-15.csv``` : Load Follows Profile Case (Random Power Level Changes between 55% and 100%)

## Last Updated
Last Updated: August 2, 2024  