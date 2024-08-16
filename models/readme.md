# VARMAX Models for systemSurrogateRGPump

This directory contains text files of the VARMAX models used in the ```systemSurrogatePump``` models. The VARMAX models are trained using the [statsmodels](https://www.statsmodels.org/dev/index.html) python package; the resulting [VARMAXResults](https://www.statsmodels.org/dev/generated/statsmodels.tsa.statespace.varmax.VARMAXResults.html) objects are then converted to text files.

### Models

|file|$\mathbf{s}$ : States|$\mathbf{u}$ : Inputs|$p$|$q$|
|------|------|------|--|--|
|```varmax-surrogate-RGPump1.txt```|['primary_pump_flow', 'solar_pump_flow']| ['target_power']|2|1|
|```varmax-surrogate-RGPump2.txt```|['Pump:pump_RPM', 'Solar-pump:pump_RPM', 'core_energy', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'Core_flow', 'secondary_flow', 'IHX_energy', 'SG_energy']|['target_power', 'primary_pump_flow', 'solar_pump_flow']|0|1|
|```varmax-surrogate-RGPump3.txt```|['core_in_P', 'IHX_primary_in_P', 'IHX_primary_out_P', 'IHX_secondary_out_t', 'IHX_secondary_in_P', 'core_out_P', 'SG_in_T', 'SG_out_T', 'SG_flow', 'CRD_pos:value'] |['primary_pump_flow', 'solar_pump_flow', 'core_energy', 'C_XE135']|0|1|
|```varmax-surrogate-RGPump4.txt```|['CRD_Reactivity', 'moderator_Reactivity', 'doppler_Reactivity', 'IHX_secondary_out_P'] |['target_power', 'CRD_pos:value', 'primary_pump_flow', 'solar_pump_flow']|0|1|
|```varmax-surrogate-RGPump5.txt```|['core_in_T', 'core_out_T', 'IHX_primary_in_T', 'IHX_primary_out_t', 'IHX_secondary_in_T', 'SG_out_P', 'SG_in_P', 'coolant_Reactivity']|['primary_pump_flow', 'solar_pump_flow', 'core_energy']|0|1|


### Versions
Last Updated : August 16, 2024
