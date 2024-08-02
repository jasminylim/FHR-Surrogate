def PP_DegradationSurrogate(
        WhichPump: str, 
        NextPumpFlowRate: float, 
        NextPumpHead: float, 
        NextPumpDensity: float, 
        PreviousK: float, 
        PreviousTime: float, 
        TimeStepSize: float, 
        NewPump: bool, 
        FullyDegradeTime: float, 
        FullyDegradeHeadLoss: float, 
        DegradeUncertainty: float, 
        DegradeInvScale: float, 
        FlowRateLossCoeff: float = 0.0, 
        PreviousPumpFlowRate: float = 0.0, 
        FullyDegradeFlowRateInv: float = 0.0,
        ShowNoDegradePower: bool = True
    ) -> tuple:
    """
    - Name: PP_DegradationSurrogate
    - Description: Homologous pump model w/ degradation using a loss coefficient
    - Author: Aaron Huxford
    - Created: 12-19-2023
    - Edits: 
        12-19-2023 by Jin Li: wrap the code into a function
        02-05-2024 by Aaron Huxford: updated to solve for pump power and more inputs
        04-02-2024 by Jin Li: updated to include degradation due to flow rate change
        05-23-2024 by Jin Li: define and calculate the pump health index

    - Normal Inputs:
        WhichPump:        pump being modeled ('Primary' or 'Secondary')   [m^3/s]
        NextPumpFlowRate: next state's pump flow rate                     [m^3/s]
        NextPumpHead:     next state's pump pressure head                 [Pa]
        NextPumpDensity:  next state's pump fluid density                 [kg/m^3]
        PreviousK:        previous state's degradation loss coefficient   [-]
        PreviousTime:     previous state's time                           [s]
        TimeStepSize:     time step from prev. state's time to next state [s]
        NewPump:          pump been overhauled since last model call?     [bool] 
        ShowNoDegradePower: return power w/o degradation?                 [bool]
        
    - Inputs for varying degradation rate:
        FullyDegradeTime: time when pump is fully degraded                [s]
        - real recommended value = 10 years in seconds (315,360,000)
        - can vary value to vary degradation rate versus time 
        
        FullyDegradeHeadLoss: fractional head loss when fully degraded    [-]
        - recommended value = 0.1
        - can vary value between 0 and 1 to vary degradation severity

        DegradeUncertainty: used to define the standard dev. in normal dist. [-]
        - recommended value = 0.5
        - standard deviation = value*dKdt_nominal
        - can vary value between 0 and 1 to change spread of random degradation
           
        DegradeInvScale: used to scale the dKdt_nominal [-]
        - recommended value = 1.0
        - can vary value between 0 and 1 to change spread of random degradation
    
    - Inputs for pump flow rate degradation loss:
        FlowRateLossCoeff: coefficient for flow rate degradation loss      [s/m^3]
    
    - Inputs for pump flow rate change degradation loss:
        PreviousPumpFlowRate: previous state's pump flow rate              [m^3/s]
        FullDegradeFlowRateInv: inverse of fully degrade flow rate change  [m^3/s]
        - recommended value = 1 / (2*0.22*3600*24*365)
        - can vary value to change degradation rate due to flow rate change

    - Outputs:
        PumpPowerDeg:    electrical power required for pump w/ degradation   [W]
        PumpPowerNoDeg:  electrical power required for pump w/o degradation  [W]
                         returned if ShowNoDegradePower is True, or else -1
        NextK:           degradation loss coefficient at the next state      [-]
        PumpHealthIndex: index indicating the pump health status             [-]
                         value is between 0 (end of life) to 1 (brand new)
                         assuming index = 0 when pump power = 2*no-degrade power
    """
    import numpy as np
    import sys
    from scipy.optimize import fsolve
    import scipy.stats as stats

    # pull inputs
    run_D = NextPumpDensity
    run_V = NextPumpFlowRate
    run_H = NextPumpHead
    dt    = TimeStepSize
    
    degrade_percent = FullyDegradeHeadLoss
    previousK       = PreviousK
    t_fulldegrade   = FullyDegradeTime

    spread    = DegradeUncertainty
    deviation = DegradeInvScale

    # if pump was overhauled since last timestep
    if NewPump:
        t_nodegrade = PreviousTime
        previousK = 0
    else:
        t_nodegrade = 0
    #--------------------------------------------------------------------------
    if WhichPump=='Primary':
        # primary pump parameters, FLiBe in gFHR model | units
        Area    = 0.47171642012  # flow area           | m^2
        rated_N = 1200           # rated speed         | rpm
        rated_V = 2*0.22         # rated flow rate     | m^3/s
        rated_H = 148188.23      # rated head          | Pa
        rated_D = 1900           # rated density       | kg/m^3
    elif WhichPump=='Secondary':
        # secondary pump parameters, Nitrate in gFHR model | units
        Area    = 0.51934454          # flow area          | m^2
        rated_N = 1200                # rated speed        | rpm
        rated_V = 2*0.3               # rated flow rate    | m^3/s
        rated_H = 74300.50            # rated head         | Pa
        rated_D = 1772                # rated density      | kg/m^3
    else:
        sys.exit('Invalid input for WhichPump, only allows Primary or Secondary as inputs')
    #--------------------------------------------------------------------------
    '''
    Functions
    '''
    def homol_head(x):
        '''
        Inputs
        - x: pi + arctan(Qbar/Nbar)
        
        Outputs
        - pump head (new pump w/o degradation)
        '''
    
        # polynomial coefficients in SAM
        ah = np.ones([3,7])
        # 0 to pi
        ah[0,0] =  0.63380980
        ah[0,1] =  0.46015764
        ah[0,2] = -2.4004049
        ah[0,3] =  3.17937240
        ah[0,4] = -1.7730449
        ah[0,5] =  0.46235776
        ah[0,6] = -0.04624640
        # pi to 3pi/2
        ah[1,0] =  431.96699
        ah[1,1] = -576.61438
        ah[1,2] =  301.00029
        ah[1,3] =  -75.465856
        ah[1,4] =    8.6754986
        ah[1,5] =   -0.26062352
        ah[1,6] =   -0.01596287
        # 3pi/2 to 2pi 
        ah[2,0] =  3240.358813
        ah[2,1] = -3096.775392
        ah[2,2] =  1184.122172
        ah[2,3] =  -226.2997
        ah[2,4] =    21.595177
        ah[2,5] =    -0.822322
        ah[2,6] =     0.0
        
        
        PI = np.pi 
        
        if 0 <= x < PI:
            i_h = 0
        elif PI <= x <= PI*3/2:
            i_h = 1
        else:
            i_h = 2
    
        return ah[i_h,0] + ah[i_h,1]*x + ah[i_h,2]*x**2 + ah[i_h,3]*x**3 \
                              + ah[i_h,4]*x**4 + ah[i_h,5]*x**5 + ah[i_h,6]*x**6 
    
    def homol_torque_hydraulic(x):
        '''
        Inputs
        - x: pi + arctan(Qbar/Nbar)
    
        Outputs
        - pump hydraulic torque (new pump w/o degradation)
        '''
        
        # polynomial coefficients used in SAM
        ah = np.ones([3,7])
        # 0 to pi
        ah[0,0] = -6.8436766E-01
        ah[0,1] =  2.7759909E+00
        ah[0,2] = -5.3988010E+00
        ah[0,3] =  6.8541205E+00
        ah[0,4] = -4.0757860E+00
        ah[0,5] =  1.0813311E+00
        ah[0,6] = -1.0475812E-01
        # pi to 3pi/2
        ah[1,0] = -1.1549471E+03
        ah[1,1] =  1.8584915E+03
        ah[1,2] = -1.2376683E+03
        ah[1,3] =  4.3601653E+02
        ah[1,4] = -8.5573772E+01
        ah[1,5] =  8.8627717E+00
        ah[1,6] = -3.7830487E-01
        # 3pi/2 to 2pi 
        ah[2,0] = -3.7981080E+02
        ah[2,1] =  7.2614914E+02
        ah[2,2] = -4.9625029E+02
        ah[2,3] =  1.6764136E+02
        ah[2,4] = -3.0366923E+01
        ah[2,5] =  2.8311896E+00
        ah[2,6] = -1.0681625E-01
        
        PI = np.pi
        
        if 0 <= x < PI:
            i_h = 0
        elif PI <= x <= PI*3/2:
            i_h = 1
        else:
            i_h = 2
    
        return ah[i_h,0] + ah[i_h,1]*x + ah[i_h,2]*x**2 + ah[i_h,3]*x**3 \
                              + ah[i_h,4]*x**4 + ah[i_h,5]*x**5 + ah[i_h,6]*x**6 
    
    
    def default_torque_friction(Nbar):
        '''
        Inputs
        - Nbar: normalized pump speed
    
        Outputs
        - pump frictional torque (new pump w/o degradation)
        '''
        
        if (Nbar >= 0 and Nbar < 0.01):
            val = 0.1 - 73.13 * Nbar**2
          
        elif (Nbar >= 0.01 and Nbar < 0.268):
            val = 0.00268 + 0.07 * Nbar**2;
    
        elif (Nbar >= 0.268):
            val = 0.00383 + 0.01071 * Nbar + 0.01406 * Nbar**2
    
        return val
    
    # Equation to solve for normalized pump speed
    def N_func(Nbar,Vbar,H):
        return rated_H * (Vbar**2 + Nbar**2) * homol_head(np.array(np.pi+np.arctan2(Vbar,Nbar))) - H 
    '''
    Calculate required pump electrical power
    '''
    bar_V = run_V/rated_V # normalized flow rate

    # define nominal degradation rate dK/dt
    dKdt_nom = deviation*(rated_H*degrade_percent/( (rated_D/2) * (rated_V/Area)**2))/(t_fulldegrade-t_nodegrade) # #/second
    fullyDegradeK = deviation*(rated_H*degrade_percent/( (rated_D/2) * (rated_V/Area)**2)) # fully degraded K

    NextK = -1
    dK_time = -1
    dK_mdot = -1
    dK_dmdot_dt = -1

    if spread > 0.0 and dKdt_nom > 0.0:
        # truncate normal distribution around lower and upper bounds for dK/dt
        lower = 0
        upper1 = 2*dKdt_nom
        mu1, sigma1 = dKdt_nom, dKdt_nom*spread
        X1 = stats.truncnorm((lower - mu1) / sigma1, (upper1 - mu1) / sigma1, loc=mu1, scale=sigma1)  

        while (dK_time < 0):
            # call truncated normal pdf
            dKdt_rand = X1.rvs(1)[0]
            # dK due to time
            dK_time = (1/deviation)*dKdt_rand*dt
    else:
        # dK due to time
        dK_time = dKdt_nom*dt

    if spread > 0.0 and FlowRateLossCoeff > 0.0:
        # truncate normal distribution around lower and upper bounds for dK/dmdot
        upper2 = 2*FlowRateLossCoeff
        mu2, sigma2 = FlowRateLossCoeff, FlowRateLossCoeff*spread
        X2 = stats.truncnorm((lower - mu2) / sigma2, (upper2 - mu2) / sigma2, loc=mu2, scale=sigma2)

        while (dK_mdot < 0):
            # dK due to flow rate
            dK_dmdot = X2.rvs(1)[0]
            dK_mdot = (1/deviation)*dK_time*dK_dmdot*(NextPumpFlowRate + PreviousPumpFlowRate)/2
    else:
        # dK due to flow rate
        dK_mdot = dK_time * FlowRateLossCoeff*(NextPumpFlowRate + PreviousPumpFlowRate)/2

    if spread > 0.0 and FullyDegradeFlowRateInv > 0.0:
        # truncate normal distribution around lower and upper bounds for dK/dmdot_dt
        upper3 = 2*FullyDegradeFlowRateInv
        mu3, sigma3 = FullyDegradeFlowRateInv, FullyDegradeFlowRateInv*spread
        X3 = stats.truncnorm((lower - mu3) / sigma3, (upper3 - mu3) / sigma3, loc=mu3, scale=sigma3) 
        
        while (dK_dmdot_dt < 0):
            # dK due to flow rate change
            dK_dmdot_dt = (1/deviation)*X3.rvs(1)[0]*fullyDegradeK*abs(NextPumpFlowRate - PreviousPumpFlowRate)
    else:
        # dK due to flow rate change
        dK_dmdot_dt = fullyDegradeK*FullyDegradeFlowRateInv*abs(NextPumpFlowRate - PreviousPumpFlowRate)

    #print(dK_mdot)
    #print(dK_time, dK_mdot, dK_dmdot_dt, NextPumpFlowRate, PreviousPumpFlowRate)
    NextK = previousK + dK_time + dK_mdot + dK_dmdot_dt

    # degradation is a function of K, D and V
    loss_H = 0.5*NextK*run_D*(run_V/Area)**2
    req_H  = run_H + loss_H
    
    bar_N = 1 # initial guess
    bar_N = fsolve(N_func,[bar_N], args=(bar_V,req_H))
    
    # get pump electrical power using pump speed
    f1 = 6  # SAM default for hydraulic torque
    f2 = 40 # SAM default for frictional torque
    
    # pump power
    PumpPowerDeg = bar_N*rated_N*( f1*(bar_V**2 + bar_N**2)*homol_torque_hydraulic(np.pi+np.arctan2(bar_V,bar_N)) + f2*default_torque_friction(bar_N) )

    PumpPowerNoDeg = -1
    if ShowNoDegradePower:
        # new pump values for reference
        bar_N_new = fsolve(N_func,[bar_N], args=(bar_V,run_H))
        PumpPowerNoDeg = bar_N_new*rated_N*( f1*(bar_V**2 + bar_N_new**2)*homol_torque_hydraulic(np.pi+np.arctan2(bar_V,bar_N_new)) + f2*default_torque_friction(bar_N_new) )
    
    # pump health index
    PumpHealthIndex = 2 - PumpPowerDeg/PumpPowerNoDeg

    return PumpPowerDeg, PumpPowerNoDeg, NextK, PumpHealthIndex