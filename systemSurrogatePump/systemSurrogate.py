# file: systemSurrogate.py
# description: systemSurrogate model for gFHR with pump degradation surrogate model
# author : Jasmin Lim
# date created: 2-15-2024

import numpy as np

# import pump degradation surrogate model
from PP_SurrogateModel import PP_DegradationSurrogate

class xenon_reactivity(object):
    '''
    - Name: xenon_reactivity
    - Description: analytical solver for the xenon concentration equations
    - Author: Jin Li
    - Created: 10-23-2023
    - Edits :
        02-15-2023 by Jasmin Lim : wrap code into object
    '''
    def __init__(self):
        # set Iodine/Xenon equation constants
        self.lambda_i = 2.9306e-05              # I-135 decay constant
        self.lambda_xe = 2.1066e-05             # Xe-135 decay constant

        self.sigma_xe1 = 1.58836E+05*1.0e-24    # Xe-135 microscopic absorption cross section, calculated by new methods (updated)

        self.phi0 = 7.36630E+14                 # total fluid cm^-2*s^-1 
        self.gamma_i = 6.38584e-02              # I-135 yield
        self.gamma_x = 4.44110e-03              # Xe-135 yield
        self.Sigma_f = 8.29802E-04              # Macroscopic fission cross section 

        self.xe_concentration_reactivity = -0.72360e-17 # Xenon concentration reactivity

        self.sf = self.sigma_xe1*self.phi0
        self.rt = self.gamma_i/(self.gamma_i+self.gamma_x)

    def xe_dydt(self, y, P):
        ''' ode function for the iodine/xenon concentration equations
        Equations 6.45 and 6.46 in SAM Theory Manual (ANL/NSE-17/4 Rev. 1)
        
        - Inputs
            y [(2,) array] : y[0] iodine relative concentration, y[1] xenon relative concentration
            P [float] : relative reactor power

        - Returns
            du/dt [(2,)] : computed ode values
        '''
        return np.array([self.lambda_i*(P-y[0]), \
                    self.lambda_xe*(P-y[1])+self.sf*P*(1-y[1])+self.rt*(self.lambda_xe+self.sf)*(y[0]-P)])
    

class homologous_pump(object):
    '''
    - Name: homologous_pump
    - Description: Homologous pump model 
    - Author: Aaron Huxford
    - Created: 12-19-2023
    - Edits:
        02-15-2024 by Jasmin Lim : wrap code into object
        06-11-2024 by Jasmin Lim : correct pump densities

    - Inputs
        type: pump being modeled ('p' for primary, 's' for secondary)
    '''
    def __init__(self, type="p"):
        self.type = type 
        
        if type not in ["p", "s"]:
            raise Exception(f"type must be either [\"p\", \"s\"] to specify primary or secondary pump respectivily")

        # assumed constants
        self.g = 9.81   # gravity [m/s^2]

        # set pump parameters 
        # primary pump
        if self.type == "p": 
            self.area    = 0.47171642012              # flow area [m^2]
            self.rated_N = 1200/1.595119e+03          # rated speed [rpm]
            self.rated_V = 2*0.22/1.147402e+03        # rated flow rate [m^3/s]
            self.rated_H = 148188.23/2.569589e+05     # rated head [Pa]
            self.rated_D = 1900                       # rated density [kg/m^3]
        # seconrady pump
        if self.type == "s":
            self.area    = 0.51934454                 # flow area                                  | m^2
            self.rated_N = 1200/3.386493e+03          # rated speed                                | rpm
            self.rated_V = 2*0.3/1.659512e+03         # rated flow rate                            | m^3/s
            self.rated_H = 74300.50/7.236768e+05      # rated head                                 | Pa
            self.rated_D = 1772                       # rated density                              | kg/m^3
        self.run_D   = self.rated_D # assume run density is approximately same as rated density


    def homologous_head(self, x):
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
    
        
        if 0 <= x < np.pi:
            i_h = 0
        elif np.pi <= x <= np.pi*3/2:
            i_h = 1
        else:
            i_h = 2

        return ah[i_h,0] + ah[i_h,1]*x + ah[i_h,2]*x**2 + ah[i_h,3]*x**3 \
                            + ah[i_h,4]*x**4 + ah[i_h,5]*x**5 + ah[i_h,6]*x**6 
    

    def compute_head(self, run_N, run_V):
        ''' compute head: computes pump head
        
        - Inputs:
            run_N [float] : pump speed [rpm]
            run_V [float] : pump mass flow rate [kg/s]
        
        - Returns
            Pump Head [Pa]'''

        run_V /= self.rated_D # convert pump flow rate from [kg/s] to [m^3/s]
        bar_N = run_N/self.rated_N # normalized speed
        bar_V = run_V/self.rated_V # normalized flow rate

        # call functions
        W_h   = self.homologous_head(np.pi + np.arctan2(bar_V,bar_N)) 

        bar_H = (bar_N**2 + bar_V**2)*W_h #normalized head

        return self.rated_H*bar_H # compute pump head

class varmax_surr_txt(object):
    def __init__(self, filename):

        '''
        - Name: varmax
        - Description: Import trained varmax model (.txt files) and use for prediction
        - Author: Jasmin Lim
        - Created: 05-12-2024
        - Edits: 
        
        - Inputs:
            filename [str] : varmax model text file name
        '''

        with open(filename, "r") as file:
            data = file.readlines()

        self.p = int(data[0].strip())
        self.q = int(data[1].strip())
        self.ns = int(data[2].strip())
        self.nu = int(data[3].strip())
        self.num_params = int(data[4].strip())

        self.states = data[5].strip().split(",")
        self.inputs = data[6].strip().split(",")
        self.param_names = data[7].strip().split(",")
        self.params = data[8].strip().split(",")
        self.params = [float(param) for param in self.params]
        self.error_cov_type = data[9].strip()
        
        # create coefficient matrices
        self.a = np.zeros(self.ns)
        for i in range(0, self.ns):
            idx = self.param_names.index('intercept.' + self.states[i])
            self.a[i] = self.params[idx]

        # AR coefficients matrix
        if self.p > 0:
            self.A = np.zeros((self.p, self.ns, self.ns))
            for i in range(0, self.p):  # axis=0 loop for lags
                for k in range(0, self.ns):  # axis=1 loop for columns
                    for j in range(0, self.ns):  # axis=2 loop for rows
                        idx = self.param_names.index('L' + str(i + 1) + '.' + self.states[k] + '.' + self.states[j])
                        self.A[i, j, k] = self.params[idx]
            self.A = np.hstack(self.A)
        else:
            self.A = np.zeros((self.ns, self.ns))

        # MA coefficient matrix
        self.M = np.zeros((self.q, self.ns, self.ns))
        for i in range(0, self.q):
            for k in range(0, self.ns):
                for j in range(0, self.ns):
                    idx = self.param_names.index('L' + str(i + 1) + '.e(' + self.states[k] + ').' + self.states[j])
                    self.M[i, j, k] = self.params[idx]
        self.M = np.hstack(self.M)

        # input coefficient matrix
        self.B = np.zeros((self.ns, self.nu))
        for i in range(self.ns):
            for j in range(self.nu):
                idx = self.param_names.index('beta.' + self.inputs[j] + "." + self.states[i])
                self.B[i,j] = self.params[idx]

        # noise covariance matrices
        self.Sigma = np.zeros((self.ns, self.ns))
        if self.error_cov_type == "diagonal":
            for i in range(self.ns):
                idx = self.param_names.index('sigma2.' + self.states[i])
                self.Sigma[i,i] = self.params[idx]
        elif self.error_cov_type == "unstructured":
            for i in range(self.ns):
                for j in range(self.ns):
                    if i==j:
                        idx = self.param_names.index('sqrt.var.' + self.states[i])
                    else:
                        idx = self.param_names.index('sqrt.cov.' + self.states[min(i,j)] + "." + self.states[min(i,j)])
                    
                    self.Sigma[i,j] = self.params.iloc[idx]
        self.Sigma = 0*self.Sigma

    def step(self, current_x, current_u, current_noise):
        '''One prediction step of VARMAX model
        
        Inputs
        ------
        current_x (ns, p) : current state values at time k
        current_u (nu, )  : current input values at time k 
        current_noise (ns, q) : current noise values at time k
        
        Returns
        --------
        res (ns) : states values at time k+1
        noise : noise values at time k+1 '''
        new_v = np.dot(self.Sigma, np.random.randn(self.ns))
        if self.p > 0:
            current_x = np.reshape(np.flip(current_x ,axis = 0), (self.p*self.ns))

        res = self.a + np.dot(self.A, current_x) + np.dot(self.B, current_u) + np.dot(self.M, current_noise) + new_v
        current_noise = np.hstack([new_v, current_noise[self.ns:]])

        return res, current_noise
    
class systemSurrogate(object):
    def __init__(self):
        '''
        - Name: systemSurrogate
        - Description: gFHR System Surrogate Model
        - Author: Jasmin Lim
        - Created: 02-13-2024
        - Edits: 
            03-03-2024 by Jasmin Lim : added states for full Task 3 requirements
            03-07-2024 by Jasmin Lim : re-ordered state vector variables, small bugs fixed
            03-20-2024 by Jasmin Lim : new surrogate models added with updated SAM data (control rod movement constraint)
        - Inputs:
            (none)

        '''

        self.V1 = varmax_surr_txt("../models/varmax-surrogate-RGPump1.txt") # varmax model #1
        self.V2 = varmax_surr_txt("../models/varmax-surrogate-RGPump2.txt") # varmax model #2
        self.V3 = varmax_surr_txt("../models/varmax-surrogate-RGPump3.txt") # varmax model #3
        self.V4 = varmax_surr_txt("../models/varmax-surrogate-RGPump4.txt") # varmax model #4
        self.V5 = varmax_surr_txt("../models/varmax-surrogate-RGPump5.txt") # varmax model #5
        self.xe = xenon_reactivity() # iodine/xenon physical model
        self.primary_pump = homologous_pump("p") # primary pump model
        self.secondary_pump = homologous_pump("s") # secondary pump model

        self.time_delay = 2
        
        # define state and inputs (SAM datafile names)
        self.states = ["time",
                        "C_I135",
                        "C_XE135",
                        "core_in_T",
                        "core_out_T",
                        "core_in_P",
                        "core_out_P",
                        "IHX_primary_in_T",
                        "IHX_primary_out_t",
                        "IHX_primary_in_P",
                        "IHX_primary_out_P",
                        "IHX_secondary_in_T",
                        "IHX_secondary_out_t",
                        "IHX_secondary_in_P",
                        "IHX_secondary_out_P",
                        "SG_in_T",
                        "SG_out_T",
                        "SG_in_P",
                        "SG_out_P",
                        "Core_flow",
                        "secondary_flow",
                        "SG_flow",
                        "core_energy",
                        "IHX_energy",
                        "SG_energy",
                        "C1",
                        "C2",
                        "C3",
                        "C4",
                        "C5",
                        "C6",
                        "moderator_Reactivity",
                        "coolant_Reactivity",
                        "doppler_Reactivity",
                        "CRD_Reactivity",
                        "CRD_pos:value",
                        "primary_pump_flow",
                        "solar_pump_flow",
                        "Pump:pump_RPM",
                        "Solar-pump:pump_RPM",
                        "Pump:pump_head",
                        "Solar-pump:pump_head"]
        self.inputs = ["time","target_power"]

        # define state and input labels (for plotting, includes units)
        self.state_labels = [r"Time [s]",
                            r"$N_I$ [$\Delta$k/k]",
                            r"$N_{Xe}$ [$\Delta$k/k]",
                            r"$T_{c,in}$ [K]",
                            r"$T_{c,out}$ [K]",
                            r"$P_{c,in}$ [Pa]",
                            r"$P_{c,out}$ [Pa]",
                            r"$T_{ihx,p,in}$ [K]",
                            r"$T_{ihx,p,out}$ [K]",
                            r"$P_{ihx,p,in}$ [Pa]",
                            r"$P_{ihx,p,out}$ [Pa]",
                            r"$T_{ihx,s,in}$ [K]",
                            r"$T_{ihx,s,out}$ [K]",
                            r"$P_{ihx,s,in}$ [Pa]",
                            r"$P_{ihx,s,out}$ [Pa]",
                            r"$T_{sg,in}$ [K]",
                            r"$T_{sg,out}$ [K]",
                            r"$P_{sg,in}$ [Pa]",
                            r"$P_{sg,out}$ [Pa]",
                            r"$\dot{m}_c$ [kg/s]",
                            r"$\dot{m}_s$ [kg/s]",
                            r"$\dot{m}_{sg}$ [kg/s]",
                            r"$\dot{Q}_{RX}$ [W]",
                            r"$\dot{Q}_{HX}$ [W]",
                            r"$\dot{Q}_{SG}$ [W]",
                            r"$C_1$",
                            r"$C_2$",
                            r"$C_3$",
                            r"$C_4$",
                            r"$C_5$",
                            r"$C_6$",
                            r"$\rho_{m}$ [$\Delta$k/k]",
                            r"$\rho_c$ [$\Delta$k/k]",
                            r"$\rho_{f}$ [$\Delta$k/k]",
                            r"$\rho_{ext}$ [$\Delta$k/k]",
                            r"$z_{cr}$ [m]",
                            r"$\dot{m}_{P,p}$ [kg/s]",
                            r"$\dot{m}_{P,s}$ [kg/s]",
                            r"$n_p$ [RPM]",
                            r"$n_s$ [RPM]",
                            r"$\Delta P_p$ [Pa]",
                            r"$\Delta P_s$ [Pa]"]
        self.input_labels = [r"Time [s]", r"Target Power [W]"]


        # surrogate model is normalized by the following total values
        self.full_power = 280e6 # full reactor power
        self.total_x = np.array([1.0,           # time
                        1331942173958500.0,     # C_I135
                        302373738645450.0,      # C_XE135
                        822.03406294783,        # core_in_T
                        923.00142227975,        # core_out_T
                        501530.6176105,         # core_in _P
                        152338.19995612,        # core_out_P
                        923.00380845784,        # IHX_primary_in_T
                        822.0016912683,         # IHX_primary_out_t
                        394936.00039368,        # IHX_primary_in_P
                        387701.7491615,         # IHX_primary_out_P
                        773.15132860419,        # IHX_secondary_in_T
                        882.79936658421,        # IHX_secondary_out_t
                        967993.37688452,        # IHX_secondary_in_P
                        1.,        # IHX_secondary_out_P
                        882.80010992894,        # SG_in_T
                        882.79982590876,        # SG_out_T
                        964777.58543592,        # SG_in_P
                        964749.99974328,        # SG_out_P
                        1147.4021545562,        # Core_flow
                        -1659.5118740692,       # secondary_flow
                        1659.5123258815,        # SG_flow
                        280000000.0,            # core_energy
                        112809775.47619,        # IHX_energy
                        279941560.79359,        # SG_energy
                        49.13847501317,         # C1
                        115.2505436281,         # C2
                        27.927296511265,        # C3
                        24.34472154563,         # C4
                        3.7252301624365,        # C5
                        0.44428547801139,       # C6
                        1.,                     # moderator_Reactivity
                        1.,                     # coolant_Reactivity
                        1.,                     # doppler_Reactivity
                        1.,                     # CRD_Reactivity
                        1.309676929559,         # CRD_pos:value
                        1147.4020330681,        # primary_pump_flow
                        1659.5118364404,        # solar_pump_flow
                        1595.1190867936,        # Pump:pump_RPM
                        3386.4927517663,        # Solar-pump:pump_RPM
                        256958.9350492,         # Pump:pump_head
                        723676.79132201])       # Solar-pump:pump_head                                
        self.total_u = np.array([1., 280000000.0]) # target power

        # initial SAM values
        self.states_fullpower = np.array([0.0,           # time
                        1331942173958500.0,     # C_I135
                        302373738645450.0,      # C_XE135
                        822.03406294783,        # core_in_T
                        923.00142227975,        # core_out_T
                        501530.6176105,         # core_in _P
                        152338.19995612,        # core_out_P
                        923.00380845784,        # IHX_primary_in_T
                        822.0016912683,         # IHX_primary_out_t
                        394936.00039368,        # IHX_primary_in_P
                        387701.7491615,         # IHX_primary_out_P
                        773.15132860419,        # IHX_secondary_in_T
                        882.79936658421,        # IHX_secondary_out_t
                        967993.37688452,        # IHX_secondary_in_P
                        965252.08799173,        # IHX_secondary_out_P
                        882.80010992894,        # SG_in_T
                        882.79982590876,        # SG_out_T
                        964777.58543592,        # SG_in_P
                        964749.99974328,        # SG_out_P
                        1147.4021545562,        # Core_flow
                        -1659.5118740692,       # secondary_flow
                        1659.5123258815,        # SG_flow
                        280000000.0,            # core_energy
                        112809775.47619,        # IHX_energy
                        279941560.79359,        # SG_energy
                        49.13847501317,         # C1
                        115.2505436281,         # C2
                        27.927296511265,        # C3
                        24.34472154563,         # C4
                        3.7252301624365,        # C5
                        0.44428547801139,       # C6
                        -2.5737495564977e-08,   # moderator_Reactivity
                        -1.7998903482294e-09,   # coolant_Reactivity
                        -6.3169877882546e-08,   # doppler_Reactivity
                        0.0,                    # CRD_Reactivity
                        1.309676929559,         # CRD_pos:value
                        1147.4020330681,        # primary_pump_flow
                        1659.5118364404,        # solar_pump_flow
                        1595.1190867936,        # Pump:pump_RPM
                        3386.4927517663,        # Solar-pump:pump_RPM
                        256958.9350492,         # Pump:pump_head
                        723676.79132201         # Solar-pump:pump_head
                        ])

    def check_inputs(self, xk, uP):
        '''
        - Name: check_inputs
        - Description: check inputs for prediction
        - Author: Jasmin Lim
        - Date Created: 04-22-2024
        - Edits:
            05-03-2024 : modify target power to set maximum power change rate to 20%/min
        - Inputs:
            xk [array: (2,42)] : time-delay previous state values (normalized)
            uP [array: (n,2)] : user input for prediction (normalized)
            
        - Returns:
            uP [(n,2)] : power input for prediction (modified if needed)'''
        
        # CHECK 1 : 2-time delays for the previous state values
        if xk.shape[0] != 2:
            raise ValueError("Input Error (xk) : must provide 2 time-delays for previous state values")
        
        # CHECK 2 : check constant time-step is 5 seconds
        if abs(5.-np.diff(xk[:,0])) > 1e-6:
            raise ValueError("Input Error (xk) : time vector must have a 5 second constant time-step")
        if abs(5.-(uP[0,0]-xk[-1,0])) > 1e-6:
            raise ValueError("InputError (xk,uP) : time step from xk to uP must be 5 seconds.")
        if abs(5.-np.diff(uP[:,0])).any() > 1e-6:
            raise ValueError("Input Error (uP) : time vector must have a 5 second constant time-step")
        
        # CHECK 3 : target power must be within [50%,100%] full power
        if (uP[:,1] < 0.5).any() or (uP[:,1] > 1.0).any():
            raise ValueError("Input Error (uP) : target power requested must be between [50%,100%] of full power, 260E+08 [W]")
        
        # CHECK 4 : maximum target power rate is 20% per minute
        if ((abs(uP[0,1]-xk[-1,22])-(.2/12)) > 1e-6):
            print("---WARNING---\n Change from current power level to requested target power (xk-->uP) exceeds the power change limit, modifying input\n-------------")
            uP = self.change_target_power(np.vstack([xk[-1,[0,22]], uP]))
            uP = uP[1:,:]
        if ((abs(np.diff(uP[:,1]))-(.2/12)) > 1e-6).any():
            print("---WARNING---\n Input target power ramp exceeds the power change limit, modifying input\n-------------")
            uP = self.change_target_power(uP)

        return uP
        
    def change_target_power(self, uP):
        '''
        - Name: change_target_power
        - Description: modifies the target power distribution to be in 20% full power/minute
        - Author: Jasmin Lim
        - Date Created: 05-01-2024
        - Edits:
        
        - Inputs:
            uP [array: (n,2)] : user input for prediction - (normalized) 
            
        - Returns:
            uP [array: (n,2)] : modified input for prediction'''
        
        # compute change in power
        powerStep = np.diff(uP[:,1])

        # determine the indicies where power is changing at a rate greater than .2/minute
        greaterThan20 = abs(powerStep)-(.2/12) > 1e-6 
        change = np.where(greaterThan20 == True)[0]
        
        # initialized modified input
        newuP = np.copy(uP)
        while(len(change) > 0):
            # determine
            startChange = change[0]

            # find next power hold 
            nextPowerChange = startChange + np.where(powerStep[startChange:] == 0.0)[0][0]

            # compute new power ramp
            power1 = newuP[startChange,1]
            power2 = uP[nextPowerChange,1]
            newRamp = np.arange(power1, power2, (power2-power1)/(abs(power2-power1))*(.2/12))

            if power2 not in newRamp: # add end power if np.arrange did not add it
                newRamp = np.append(newRamp, power2)

            # modify the power input
            if startChange+newRamp.shape[0] <= newuP.shape[0]: # if new ramp fits within user input
                newuP[startChange:startChange+newRamp.shape[0],1] = newRamp
            else: # if new ramp exceeds the user inputted time
                newuP[startChange:,1] = newRamp[:newuP[startChange:].shape[0]]

            # remove indicies from current power change segment
            change = change[change>nextPowerChange]
        
        return newuP

    def predict(self, 
                xk : float,
                uP : float):
        '''
        - Name: predict
        - Description: Prediction for gFHR SYSTEM surrogate
        - Author: Jasmin Lim
        - Date Created: 02-15-2024
        - Edits:
            03-25-2024 : (Jasmin Lim) time-delay, captures PID micro-dynamics
            04-23-2024 : (Jasmin Lim) added pump degradation model
            05-03-2024 : (Jasmin Lim) updated for most recently updated PP_DegradationSurrogate, added deltaP for option to not use PP_DegradationSurrogate each time step. Seperate inputs for each pump.
        
        - Inputs
            - xk [(2,42)]: current state vector
            - uP [(n,2)]: input for prediction
            
        - Returns
            res [(n,48)]: array with system surrogate prediction results
        '''
        num_steps = uP.shape[0]

        # normalize state and input values
        xk_ = xk/self.total_x
        uP_ = uP/self.total_u

        # check inputs
        uP_ = self.check_inputs(xk_, uP_)

        # set time-step
        dt = 5.
        
        # initialize system surrogate results
        res = np.zeros((num_steps, len(self.states)))
        res[:,0] = np.copy(uP_[:,0])
        current_x = np.copy(xk_[:,[36,37]])

        # initialize noise for VARMAX models
        noise1 = np.zeros(self.V1.ns)
        noise2 = np.zeros(self.V2.ns)
        noise3 = np.zeros(self.V3.ns)
        noise4 = np.zeros(self.V4.ns)
        noise5 = np.zeros(self.V5.ns)

        for i in range(num_steps):

            ### XENON PHYSICAL MODEL - xenon concentration
            if i==0:
                res[i,[1,2]] = xk_[-1,[1,2]] + dt*self.xe.xe_dydt(xk_[-1,[1,2]], xk_[-1,22])
            else:
                res[i,[1,2]] = res[i-1,[1,2]] + dt*self.xe.xe_dydt(res[i-1,[1,2]], res[i-1,22])

            ### VARMAX MODEL 1 - time delay, pump mass flow rate
            res[i,[36,37]], noise1 = self.V1.step(current_x, [uP_[i,1]], noise1)
            current_x = np.vstack([current_x[1:,:], res[i,[36,37]]])

            ### VARMAX MODEL 2 
            res[i,[38,39,22,25,26,27,28,29,30,19,20,23,24]], noise2 = self.V2.step(res[i,[38,39,22,25,26,27,28,29,30,19,20,23,24]], np.insert(res[i,[36,37]],0,uP_[i,1]), noise2)

            ### VARMAX MODEL 3 
            res[i,[5,9,10,12,13,6,15,16,21,35]], noise3 = self.V3.step(res[i,[5,9,10,12,13,6,15,16,21,35]], res[i,[36,37,22,2]], noise3)

            ## VARMAX MODEL 4 
            res[i,[34,31,33,14]], noise4 = self.V4.step(res[i,[34,31,33,14]], np.insert(res[i,[35,36,37]],0,uP_[i,1]), noise4)

            ### VARMAX MODEL 5 
            res[i,[3,4,7,8,11,18,17,32]], noise5 = self.V5.step(res[i,[3,4,7,8,11,18,17,32]], res[i,[36,37,22]], noise5)

            ### PUMP PHYSICAL MODEL - pump head
            res[i,40] = self.primary_pump.compute_head(res[i,38], res[i,36])
            res[i,41] = self.secondary_pump.compute_head(res[i,39], res[i,37])

        return res*self.total_x
    
class systemSurrogatePump(object):
    def __init__(self):
        '''
        - Name: systemSurrogatePump
        - Description: gFHR System Surrogate Model with Pump Degradation Model
        - Author: Jasmin Lim
        - Created: 02-13-2024
        - Edits: 
            03-03-2024 by Jasmin Lim : added states for full Task 3 requirements
            03-07-2024 by Jasmin Lim : re-ordered state vector variables, small bugs fixed
            03-20-2024 by Jasmin Lim : new surrogate models added with updated SAM data (control rod movement constraint)
            05-03-2024 by Jasmin Lim : seperated system states and pump degradation states, added change target power
        - Inputs:
            (none)

        '''

        self.V1 = varmax_surr_txt("../models/varmax-surrogate-RGPump1.txt") # varmax model #1
        self.V2 = varmax_surr_txt("../models/varmax-surrogate-RGPump2.txt") # varmax model #2
        self.V3 = varmax_surr_txt("../models/varmax-surrogate-RGPump3.txt") # varmax model #3
        self.V4 = varmax_surr_txt("../models/varmax-surrogate-RGPump4.txt") # varmax model #4
        self.V5 = varmax_surr_txt("../models/varmax-surrogate-RGPump5.txt") # varmax model #5
        self.xe = xenon_reactivity() # iodine/xenon physical model
        self.primary_pump = homologous_pump("p") # primary pump model
        self.secondary_pump = homologous_pump("s") # secondary pump model

        self.time_delay = 2
        
        # define state and inputs (SAM datafile names)
        self.states = ["time",
                        "C_I135",
                        "C_XE135",
                        "core_in_T",
                        "core_out_T",
                        "core_in_P",
                        "core_out_P",
                        "IHX_primary_in_T",
                        "IHX_primary_out_t",
                        "IHX_primary_in_P",
                        "IHX_primary_out_P",
                        "IHX_secondary_in_T",
                        "IHX_secondary_out_t",
                        "IHX_secondary_in_P",
                        "IHX_secondary_out_P",
                        "SG_in_T",
                        "SG_out_T",
                        "SG_in_P",
                        "SG_out_P",
                        "Core_flow",
                        "secondary_flow",
                        "SG_flow",
                        "core_energy",
                        "IHX_energy",
                        "SG_energy",
                        "C1",
                        "C2",
                        "C3",
                        "C4",
                        "C5",
                        "C6",
                        "moderator_Reactivity",
                        "coolant_Reactivity",
                        "doppler_Reactivity",
                        "CRD_Reactivity",
                        "CRD_pos:value",
                        "primary_pump_flow",
                        "solar_pump_flow",
                        "Pump:pump_RPM",
                        "Solar-pump:pump_RPM",
                        "Pump:pump_head",
                        "Solar-pump:pump_head"]
        self.pump_states = ["time",
                        "PrimaryPumpPowerDeg", 
                        "PrimaryPumpPowerNoDeg", 
                        "PrimaryK", 
                        "PrimaryHealthIndex",
                        "SecondaryPumpPowerDeg", 
                        "SecondaryPumpPowerNoDeg", 
                        "SecondaryK",
                        "Secondary Health Index"]
        self.inputs = ["time","target_power"]

        # define state and input labels (for plotting, includes units)
        self.state_labels = [r"t [s]",
                            r"$N_I$ [$\Delta$k/k]",
                            r"$N_{Xe}$ [$\Delta$k/k]",
                            r"$T_{c,in}$ [K]",
                            r"$T_{c,out}$ [K]",
                            r"$P_{c,in}$ [Pa]",
                            r"$P_{c,out}$ [Pa]",
                            r"$T_{ihx,p,in}$ [K]",
                            r"$T_{ihx,p,out}$ [K]",
                            r"$P_{ihx,p,in}$ [Pa]",
                            r"$P_{ihx,p,out}$ [Pa]",
                            r"$T_{ihx,s,in}$ [K]",
                            r"$T_{ihx,s,out}$ [K]",
                            r"$P_{ihx,s,in}$ [Pa]",
                            r"$P_{ihx,s,out}$ [Pa]",
                            r"$T_{sg,in}$ [K]",
                            r"$T_{sg,out}$ [K]",
                            r"$P_{sg,in}$ [Pa]",
                            r"$P_{sg,out}$ [Pa]",
                            r"$\dot{m}_c$ [kg/s]",
                            r"$\dot{m}_s$ [kg/s]",
                            r"$\dot{m}_{sg}$ [kg/s]",
                            r"$\dot{Q}_{RX}$ [W]",
                            r"$\dot{Q}_{HX}$ [W]",
                            r"$\dot{Q}_{SG}$ [W]",
                            r"$C_1$",
                            r"$C_2$",
                            r"$C_3$",
                            r"$C_4$",
                            r"$C_5$",
                            r"$C_6$",
                            r"$\rho_{m}$ [$\Delta$k/k]",
                            r"$\rho_c$ [$\Delta$k/k]",
                            r"$\rho_{f}$ [$\Delta$k/k]",
                            r"$\rho_{ext}$ [$\Delta$k/k]",
                            r"$z_{cr}$ [m]",
                            r"$\dot{m}_{P,p}$ [kg/s]",
                            r"$\dot{m}_{P,s}$ [kg/s]",
                            r"$n_p$ [RPM]",
                            r"$n_s$ [RPM]",
                            r"$\Delta P_p$ [Pa]",
                            r"$\Delta P_s$ [Pa]"]
        self.pump_state_labels = [r"t [s]",
                                r"$\dot{Q}_{P,D,p}$ [W]", 
                                r"$\dot{Q}_{P,p}$ [W]", 
                                r"$K_{p}$", 
                                r"$\eta_{p}$",
                                r"$\dot{Q}_{P,D,s}$ [W]", 
                                r"$\dot{Q}_{P,s}$ [W]", 
                                r"$K_{s}$",
                                r"$\eta_{s}$"]
        self.input_labels = [r"t [s]", r"Target Power [W]"]


        # surrogate model is normalized by the following total values
        self.full_power = 280e6 # full reactor power
        self.total_x = np.array([1.0,           # time
                        1331942173958500.0,     # C_I135
                        302373738645450.0,      # C_XE135
                        822.03406294783,        # core_in_T
                        923.00142227975,        # core_out_T
                        501530.6176105,         # core_in _P
                        152338.19995612,        # core_out_P
                        923.00380845784,        # IHX_primary_in_T
                        822.0016912683,         # IHX_primary_out_t
                        394936.00039368,        # IHX_primary_in_P
                        387701.7491615,         # IHX_primary_out_P
                        773.15132860419,        # IHX_secondary_in_T
                        882.79936658421,        # IHX_secondary_out_t
                        967993.37688452,        # IHX_secondary_in_P
                        1.,        # IHX_secondary_out_P
                        882.80010992894,        # SG_in_T
                        882.79982590876,        # SG_out_T
                        964777.58543592,        # SG_in_P
                        964749.99974328,        # SG_out_P
                        1147.4021545562,        # Core_flow
                        -1659.5118740692,       # secondary_flow
                        1659.5123258815,        # SG_flow
                        280000000.0,            # core_energy
                        112809775.47619,        # IHX_energy
                        279941560.79359,        # SG_energy
                        49.13847501317,         # C1
                        115.2505436281,         # C2
                        27.927296511265,        # C3
                        24.34472154563,         # C4
                        3.7252301624365,        # C5
                        0.44428547801139,       # C6
                        1.,                     # moderator_Reactivity
                        1.,                     # coolant_Reactivity
                        1.,                     # doppler_Reactivity
                        1.,                     # CRD_Reactivity
                        1.309676929559,         # CRD_pos:value
                        1147.4020330681,        # primary_pump_flow
                        1659.5118364404,        # solar_pump_flow
                        1595.1190867936,        # Pump:pump_RPM
                        3386.4927517663,        # Solar-pump:pump_RPM
                        256958.9350492,         # Pump:pump_head
                        723676.79132201])       # Solar-pump:pump_head                                
        self.total_u = np.array([1., 280000000.0]) # target power

        # initial SAM values
        self.states_fullpower = np.array([0.0,           # time
                        1331942173958500.0,     # C_I135
                        302373738645450.0,      # C_XE135
                        822.03406294783,        # core_in_T
                        923.00142227975,        # core_out_T
                        501530.6176105,         # core_in _P
                        152338.19995612,        # core_out_P
                        923.00380845784,        # IHX_primary_in_T
                        822.0016912683,         # IHX_primary_out_t
                        394936.00039368,        # IHX_primary_in_P
                        387701.7491615,         # IHX_primary_out_P
                        773.15132860419,        # IHX_secondary_in_T
                        882.79936658421,        # IHX_secondary_out_t
                        967993.37688452,        # IHX_secondary_in_P
                        965252.08799173,        # IHX_secondary_out_P
                        882.80010992894,        # SG_in_T
                        882.79982590876,        # SG_out_T
                        964777.58543592,        # SG_in_P
                        964749.99974328,        # SG_out_P
                        1147.4021545562,        # Core_flow
                        -1659.5118740692,       # secondary_flow
                        1659.5123258815,        # SG_flow
                        280000000.0,            # core_energy
                        112809775.47619,        # IHX_energy
                        279941560.79359,        # SG_energy
                        49.13847501317,         # C1
                        115.2505436281,         # C2
                        27.927296511265,        # C3
                        24.34472154563,         # C4
                        3.7252301624365,        # C5
                        0.44428547801139,       # C6
                        -2.5737495564977e-08,   # moderator_Reactivity
                        -1.7998903482294e-09,   # coolant_Reactivity
                        -6.3169877882546e-08,   # doppler_Reactivity
                        0.0,                    # CRD_Reactivity
                        1.309676929559,         # CRD_pos:value
                        1147.4020330681,        # primary_pump_flow
                        1659.5118364404,        # solar_pump_flow
                        1595.1190867936,        # Pump:pump_RPM
                        3386.4927517663,        # Solar-pump:pump_RPM
                        256958.9350492,         # Pump:pump_head
                        723676.79132201         # Solar-pump:pump_head
                        ])

    def check_inputs(self, xk, uP):
        '''
        - Name: check_inputs
        - Description: check inputs for prediction
        - Author: Jasmin Lim
        - Date Created: 04-22-2024
        - Edits:
            05-03-2024 : modify target power to set maximum power change rate to 20%/min
        - Inputs:
            xk [array: (2,42)] : time-delay previous state values (normalized)
            uP [array: (n,2)] : user input for prediction (normalized)
            
        - Returns:
            uP [(n,2)] : power input for prediction (modified if needed)'''
        
        # CHECK 1 : 2-time delays for the previous state values
        if xk.shape[0] != 2:
            raise ValueError("Input Error (xk) : must provide 2 time-delays for previous state values")
        
        # CHECK 2 : check constant time-step is 5 seconds
        if abs(5.-np.diff(xk[:,0])) > 1e-6:
            raise ValueError("Input Error (xk) : time vector must have a 5 second constant time-step")
        if abs(5.-(uP[0,0]-xk[-1,0])) > 1e-6:
            raise ValueError("InputError (xk,uP) : time step from xk to uP must be 5 seconds.")
        if abs(5.-np.diff(uP[:,0])).any() > 1e-6:
            raise ValueError("Input Error (uP) : time vector must have a 5 second constant time-step")
        
        # CHECK 3 : target power must be within [50%,100%] full power
        if (uP[:,1] < 0.5).any() or (uP[:,1] > 1.0).any():
            raise ValueError("Input Error (uP) : target power requested must be between [50%,100%] of full power, 260E+08 [W]")
        
        # CHECK 4 : maximum target power rate is 20% per minute
        if ((abs(uP[0,1]-xk[-1,22])-(.2/12)) > 1e-6):
            print("---WARNING---\n Change from current power level to requested target power (xk-->uP) exceeds the power change limit, modifying input\n-------------")
            uP = self.change_target_power(np.vstack([xk[-1,[0,22]], uP]))
            uP = uP[1:,:]
        if ((abs(np.diff(uP[:,1]))-(.2/12)) > 1e-6).any():
            print("---WARNING---\n Input target power ramp exceeds the power change limit, modifying input\n-------------")
            uP = self.change_target_power(uP)

        return uP
        
        
    def change_target_power(self, uP):
        '''
        - Name: change_target_power
        - Description: modifies the target power distribution to be in 20% full power/minute
        - Author: Jasmin Lim
        - Date Created: 05-01-2024
        - Edits:
        
        - Inputs:
            uP [array: (n,2)] : user input for prediction (normalized)
            
        - Returns:
            uP [array: (n,2)] : modified input for prediction'''
        
        # compute change in power
        powerStep = np.diff(uP[:,1])

        # determine the indicies where power is changing at a rate greater than .2/minute
        greaterThan20 = abs(powerStep)-(.2/12) > 1e-6 
        change = np.where(greaterThan20 == True)[0]
        
        # initialized modified input
        newuP = np.copy(uP)
        while(len(change) > 0):
            # determine
            startChange = change[0]

            # find next power hold 
            nextPowerChange = startChange + np.where(powerStep[startChange:] == 0.0)[0][0]

            # compute new power ramp
            power1 = newuP[startChange,1]
            power2 = uP[nextPowerChange,1]
            newRamp = np.arange(power1, power2, (power2-power1)/(abs(power2-power1))*(.2/12))

            if power2 not in newRamp: # add end power if np.arrange did not add it
                newRamp = np.append(newRamp, power2)

            # modify the power input
            if startChange+newRamp.shape[0] <= newuP.shape[0]: # if new ramp fits within user input
                newuP[startChange:startChange+newRamp.shape[0],1] = newRamp
            else: # if new ramp exceeds the user inputted time
                newuP[startChange:,1] = newRamp[:newuP[startChange:].shape[0]]

            # remove indicies from current power change segment
            change = change[change>nextPowerChange]
        
        return newuP

    def predict(self, 
                xk : float, 
                xkPump : float,
                uP : float, 
                newPump : dict = {"Primary":None, "Secondary":None}, 
                FullyDegradeTime : dict = {"Primary":10*365*24*3600,"Secondary":10*365*24*3600},
                FullyDegradeHeadLoss : dict = {"Primary": 0.5, "Secondary": 0.5}, 
                DegradeUncertainty : dict = {"Primary": 0.001, "Secondary": 0.001}, 
                DegradeInvScale : dict = {"Primary": 0., "Secondary": 0.}, 
                FlowRateLossCoeff : dict = {"Primary": 0., "Secondary": 0.},
                FullDegradeFlowRateInv : dict = {"Primary": 0., "Secondary": 0.},
                deltaP : int = 1,
                ShowNoDegradePower : bool = False):
        '''
        - Name: predict
        - Description: Prediction for gFHR SYSTEM surrogate
        - Author: Jasmin Lim
        - Date Created: 02-15-2024
        - Edits:
            03-25-2024 : (Jasmin Lim) time-delay, captures PID micro-dynamics
            04-23-2024 : (Jasmin Lim) added pump degradation model
            05-03-2024 : (Jasmin Lim) updated for most recently updated PP_DegradationSurrogate, added deltaP for option to not use PP_DegradationSurrogate each time step. Seperate inputs for each pump.
            05-23-2024 : (Jasmin Lim) update function for new pump surrogate output, health index score
        
        - Inputs
            - xk [(2,42)]: current state vector
            - xkPump [(9,)] : current pump degradation vector
            - uP [(n,2)]: input for prediction
            - (optional) newPump [dict, (n,)]: array indicating new pump. "1" for new pump, "0" for no new pump
            - (optional) FullyDegradeTime [dict]: full degrade time. ["Primary","Secondary"] keys define each pump respectivly. Default is 10 years.
            - (optional) FullyDegradeHeadLoss [dict]: full degraded head loss. ["Primary","Secondary"] keys define each pump respectivly. Default is 0.1
            - (optional) DegradeUncertainty [dict]: pump degradation uncertainity. ["Primary","Secondary"] keys define each pump respectivly. Default is 0.5
            - (optional) DegradeInvScale [dict]: pump power inverse problem uncertainty. ["Primary","Secondary"] keys define each pump respectivly. Default is 0.001
            - (optional) FlowRateLossCoeff [dict] : inverse of fully degraded flow rate change. ["Primary","Secondary"] keys define each pump respectivly. Default is 0
            - (optional) FullDegradeFlowRateInv [float] : inverse of fully degraded flow rate change. ["Primary","Secondary"] keys define each pump respectivly. Default is 0
            - (optional) deltaP [int]: sets compute step for pump degradation surrogate model. Default is 1
            - (optional) ShowNoDegradePower [bool]: If True, will include pump power without degradation computation pump degradation surrograte model. Default is False
            
        - Returns
            res [(n,42)]: array with system surrogate prediction results
            pump_res [(n//deltaP,9)] : array with pump degradation surrogate results
        '''
        num_steps = uP.shape[0]

        # normalize state and input values
        xk_ = xk/self.total_x
        uP_ = uP/self.total_u

        # check inputs
        uP_ = self.check_inputs(xk_, uP_)

        # set time-step
        dt = 5.
        
        # initialize system surrogate results
        res = np.zeros((num_steps, len(self.states)))
        res[:,0] = np.copy(uP_[:,0])
        current_x = np.copy(xk_[:,[36,37]])

        # initialize pump degradation surrogate results
        pumpRes = np.zeros((int(np.ceil(num_steps/deltaP)), 9))
        pumpRes[:,0] = uP[::deltaP,0]

        # initialize noise for VARMAX models
        noise1 = np.zeros(self.V1.ns)
        noise2 = np.zeros(self.V2.ns)
        noise3 = np.zeros(self.V3.ns)
        noise4 = np.zeros(self.V4.ns)
        noise5 = np.zeros(self.V5.ns)

        for i in range(num_steps):

            ### XENON PHYSICAL MODEL - xenon concentration
            if i==0:
                res[i,[1,2]] = xk_[-1,[1,2]] + dt*self.xe.xe_dydt(xk_[-1,[1,2]], xk_[-1,22])
            else:
                res[i,[1,2]] = res[i-1,[1,2]] + dt*self.xe.xe_dydt(res[i-1,[1,2]], res[i-1,22])

            ### VARMAX MODEL 1 - time delay, pump mass flow rate
            res[i,[36,37]], noise1 = self.V1.step(current_x, [uP_[i,1]], noise1)
            current_x = np.vstack([current_x[1:,:], res[i,[36,37]]])

            ### VARMAX MODEL 2 
            res[i,[38,39,22,25,26,27,28,29,30,19,20,23,24]], noise2 = self.V2.step(res[i,[38,39,22,25,26,27,28,29,30,19,20,23,24]], np.insert(res[i,[36,37]],0,uP_[i,1]), noise2)

            ### VARMAX MODEL 3 
            res[i,[5,9,10,12,13,6,15,16,21,35]], noise3 = self.V3.step(res[i,[5,9,10,12,13,6,15,16,21,35]], res[i,[36,37,22,2]], noise3)

            ## VARMAX MODEL 4 
            res[i,[34,31,33,14]], noise4 = self.V4.step(res[i,[34,31,33,14]], np.insert(res[i,[35,36,37]],0,uP_[i,1]), noise4)

            ### VARMAX MODEL 5 
            res[i,[3,4,7,8,11,18,17,32]], noise5 = self.V5.step(res[i,[3,4,7,8,11,18,17,32]], res[i,[36,37,22]], noise5)

            ### PUMP PHYSICAL MODEL - pump head
            res[i,40] = self.primary_pump.compute_head(res[i,38], res[i,36])
            res[i,41] = self.secondary_pump.compute_head(res[i,39], res[i,37])

            if (i%deltaP == 0):
                iP = i//deltaP

                ### PUMP PHYSICAL MODEL - pump power
                if i == 0: # for first step, use user input of last state values
                    PrimaryK = xkPump[3]; SecondaryK = xkPump[7]
                    PrimaryV = xk_[-1,36] * self.total_x[36] / self.primary_pump.rated_D 
                    SecondaryV = xk_[-1,37] * self.total_x[37] / self.secondary_pump.rated_D
                else:
                    PrimaryK = pumpRes[iP-1,3]; SecondaryK= pumpRes[iP-1,7] 
                    PrimaryV = res[i-deltaP,36] * self.total_x[36] / self.primary_pump.rated_D 
                    SecondaryV = res[i-deltaP,37] * self.total_x[37] / self.secondary_pump.rated_D

                pumpRes[iP,0] = res[i,0] # get current time
                pumpRes[iP,1], pumpRes[iP,2], pumpRes[iP,3], pumpRes[iP,4] = PP_DegradationSurrogate("Primary", 
                                            res[i,36] * self.total_x[36] / self.primary_pump.rated_D, # volumnetric flow rate
                                            res[i,40] * self.total_x[40], # pump head
                                            self.primary_pump.rated_D, # pump density
                                            PrimaryK, # previous degradation value
                                            res[i,0], # time
                                            deltaP*dt,
                                            0 if type(newPump["Primary"])==type(None) else newPump["Primary"][iP], 
                                            FullyDegradeTime["Primary"], 
                                            FullyDegradeHeadLoss["Primary"], 
                                            DegradeUncertainty["Primary"],
                                            DegradeInvScale["Primary"],
                                            FlowRateLossCoeff["Primary"],
                                            PrimaryV,
                                            FullDegradeFlowRateInv["Primary"],
                                            ShowNoDegradePower)
                pumpRes[iP,5], pumpRes[iP,6], pumpRes[iP,7], pumpRes[iP,8] = PP_DegradationSurrogate("Secondary", 
                                            res[i,37] * self.total_x[37] / self.secondary_pump.rated_D, # volumnetric flow rate
                                            res[i,41] * self.total_x[41], # pump head
                                            self.secondary_pump.rated_D, # pump density
                                            SecondaryK, # previous pump degradation value
                                            res[i,0], # time
                                            deltaP*dt,
                                            0 if type(newPump["Secondary"])==type(None) else newPump["Secondary"][iP], 
                                            FullyDegradeTime["Secondary"], 
                                            FullyDegradeHeadLoss["Secondary"], 
                                            DegradeUncertainty["Secondary"],
                                            DegradeInvScale["Secondary"],
                                            FlowRateLossCoeff["Secondary"],
                                            SecondaryV,
                                            FullDegradeFlowRateInv["Secondary"],
                                            ShowNoDegradePower)
                
        return res*self.total_x, pumpRes