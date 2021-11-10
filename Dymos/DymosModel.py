# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 16:28:42 2021

@author: renkert2
"""
import dymos as dm
import numpy as np
import openmdao.api as om
import importlib as impL
import json

def DymosPhase(mdl, tx, state_opts = {}, control_opts = {}, disturbance_opts = {}, parameter_opts = {}):
    # Instantiate a Dymos Trajectory
    meta = ImportMetadata(mdl)
    
    # Instantiate a Phase and add it to the Trajectory.
    # Here the transcription is necessary but not particularly relevant.
    phase = dm.Phase(ode_class=DymosModel, ode_init_kwargs={"Model":mdl}, transcription=tx)
    
    # Tell Dymos the states to be propagated using the given ODE.
    state_vars = [x["StateVariable"] for x in meta["StateTable"]]
    for var in state_vars:
        phase.add_state(var, rate_source=var+"_dot", **state_opts)
        
    # Tell Dymos the inputs to be propagated using the given ODE.
    input_vars = [x["InputVariable"] for x in meta["InputTable"]]
    for var in input_vars:
        phase.add_control(var, **control_opts)
        
    dist_vars = [x["DisturbanceVariable"] for x in meta["DisturbanceTable"]]
    for var in dist_vars:
        phase.add_control(var, **disturbance_opts)

    # Define constant parameters
    for param in meta["ParamTable"]:
        phase.add_parameter(param["SymID"], val=param["Value"], static_target=True, **parameter_opts)
        
    return phase

class DymosModel(om.Group):
    def initialize(self):
        self.options.declare('num_nodes', types=int) # Number of nodes property, required for Dymos models
        self.options.declare('Model', types=str, default = 'None') # "Model" property used to store name of folder containing Model .py functions and variable tables
        
    def setup(self):
        nn = self.options['num_nodes']
        mdl = self.options['Model']
        
        self.ImportModel(mdl)
        
        # Instantiate Dynamic Model Subsystem (f), promote x,u,d,theta
        f_subsys = _CalcF(num_nodes = nn, _Calc = self.Calc["f"], _CalcJ = self.CalcJ["f"], _Metadata=self.Metadata)
        self.add_subsystem(name="CalcF", subsys=f_subsys,
                           promotes = ["*"])
        
        # Instantiate Output Model Subsystem (g), promote x,u,d,theta
        g_subsys = _CalcG(num_nodes = nn, _Calc = self.Calc["g"], _CalcJ = self.CalcJ["g"], _Metadata=self.Metadata)
        self.add_subsystem(name="CalcG", subsys=g_subsys,
                           promotes = ["*"])
        
    def ImportModel(self,mdl):
        # sys.path.insert(0, mdl)
        meta_dict = ImportMetadata(mdl)

        Nx = meta_dict["Nx"]        
        Nu = meta_dict["Nu"]
        Nd = meta_dict["Nd"]
        Ny = meta_dict["Ny"]
        Ntheta = meta_dict["Ntheta"]
        
        func_list = ['f', 'g']
        calc = {}
        for func in func_list:
            mod = impL.import_module(mdl + ".Model_" + func)
            handle = getattr(mod, "Calc_" + func)
            calc[func] = handle
            
        var_list = ["x", "u", "d", "theta"]
        calcJ = {}
        for func in func_list:
            calcJ_inner = {}
            for var in var_list:
                partial = func + "_" + var
                mod = impL.import_module(mdl + ".ModelJ_" + partial)
                handle = getattr(mod, "CalcJ_" + partial)
                calcJ_inner[var] = handle
            calcJ[func] = calcJ_inner
                    
        self.Metadata = meta_dict
        self.Calc = calc
        self.CalcJ = calcJ

class _CalcF(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int) # Number of nodes property, required for Dymos models
        self.options.declare('_Metadata')
        self.options.declare('_Calc') # Number of nodes property, required for Dymos models
        self.options.declare('_CalcJ') # Number of nodes property, required for Dymos models
        
    def setup(self):
        nn = self.options['num_nodes']
        meta = self.options["_Metadata"]
        
        ### INPUTS ###
        # x
        state_table = meta["StateTable"]
        self.StateVars = [x["StateVariable"] for x in state_table]
        for var in state_table:
            self.add_input(var["StateVariable"], desc=var["Description"], shape=(nn,))
        # u
        input_table = meta["InputTable"]
        self.InputVars = [x["InputVariable"] for x in input_table]
        for var in input_table:
            self.add_input(var["InputVariable"], desc=var["Description"], shape=(nn,))
        # d
        disturbance_table = meta["DisturbanceTable"]
        self.DisturbanceVars = [x["DisturbanceVariable"] for x in disturbance_table]
        for var in disturbance_table:
            self.add_input(var["DisturbanceVariable"], desc=var["Description"], shape=(nn,))
            
        # theta
        param_table = meta["ParamTable"]
        self.ParamVars = [x["SymID"] for x in param_table]
        for param in param_table:
            self.add_input(param["SymID"], val=param["Value"], desc=var["Description"], tags=['dymos.static_target'])
        
        ### OUTPUTS ###
        #xdot
        self.StateVars_dot = [x["StateVariable"]+"_dot" for x in state_table]
        for var in state_table:
            self.add_output(var["StateVariable"]+"_dot", desc="Rate of Change: "+var["Description"], shape=(nn,))
            
        ### PARTIALS ###
        # Come Back to this - just get the model working
        self.declare_partials(["*"], ["*"], method="fd")
       
    def compute(self, inputs, outputs):
        nn = self.options['num_nodes']
        
        # Assemble Vectors
        X = []
        for i,x in enumerate(self.StateVars):
            X.append(inputs[x])
        
        U = []
        for i,u in enumerate(self.InputVars):
            U.append(inputs[u])

        D = []
        for i,d in enumerate(self.DisturbanceVars):
            D.append(inputs[d])
            
        Theta = []
        for i,theta in enumerate(self.ParamVars):
            Theta.append(inputs[theta])
        
        # Call Calc function
        f = self.options["_Calc"]
        X_dot = f(X,U,D,Theta,nn)
        
        # Process Outputs
        # - Output is list of arrays.  Outer level is the number of time segments, inner index corresponds to the outputs
        # - X_dot[i][j] corresponds to the j_th output of the i_th time segment
        
        # Assign to Outputs

        for i,x_dot in enumerate(self.StateVars_dot):
            # 
            out_temp = np.zeros(nn)
            for j in range(nn):
                out_temp[j] = X_dot[j][i]
            outputs[x_dot] = out_temp

class _CalcG(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int) # Number of nodes property, required for Dymos models
        self.options.declare('_Metadata')
        self.options.declare('_Calc') # Number of nodes property, required for Dymos models
        self.options.declare('_CalcJ') # Number of nodes property, required for Dymos models
        
    def setup(self):
        nn = self.options['num_nodes']
        meta = self.options["_Metadata"]
        
        ### INPUTS ###
        # x
        state_table = meta["StateTable"]
        self.StateVars = [x["StateVariable"] for x in state_table]
        for var in state_table:
            self.add_input(var["StateVariable"], desc=var["Description"], shape=(nn,))
        # u
        input_table = meta["InputTable"]
        self.InputVars = [x["InputVariable"] for x in input_table]
        for var in input_table:
            self.add_input(var["InputVariable"], desc=var["Description"], shape=(nn,))
        # d
        disturbance_table = meta["DisturbanceTable"]
        self.DisturbanceVars = [x["DisturbanceVariable"] for x in disturbance_table]
        for var in disturbance_table:
            self.add_input(var["DisturbanceVariable"], desc=var["Description"], shape=(nn,))
            
        # theta
        param_table = meta["ParamTable"]
        self.ParamVars = [x["SymID"] for x in param_table]
        for param in param_table:
            self.add_input(param["SymID"], val=param["Value"], desc=var["Description"], tags=['dymos.static_target'])
        
        ### OUTPUTS ###
        #y
        output_table = meta["OutputTable"]
        self.OutputVars = [x["OutputVariable"] for x in output_table]
        for var in output_table:
            self.add_output(var["OutputVariable"], desc=var["Description"], shape=(nn,))
            
        ### PARTIALS ###
        # Come Back to this - just get the model working
        self.declare_partials(["*"], ["*"], method="fd")
    
    def compute(self, inputs, outputs):
        nn = self.options['num_nodes']
        
        # Assemble Vectors
        X = []
        for i,x in enumerate(self.StateVars):
            X.append(inputs[x])
        
        U = []
        for i,u in enumerate(self.InputVars):
            U.append(inputs[u])

        D = []
        for i,d in enumerate(self.DisturbanceVars):
            D.append(inputs[d])
            
        Theta = []
        for i,theta in enumerate(self.ParamVars):
            Theta.append(inputs[theta])
        
        # Call Calc function
        g = self.options["_Calc"]
        Y = g(X,U,D,Theta,nn)
        
        # Process Outputs
        # - Output is list of arrays.  Outer level is the number of time segments, inner index corresponds to the outputs
        # - X_dot[i][j] corresponds to the j_th output of the i_th time segment
        
        # Assign to Outputs
        for i,y in enumerate(self.OutputVars):
            out_temp = np.zeros(nn)
            for j in range(nn):
                out_temp[j] = Y[j][i]
            outputs[y] = out_temp
        
        return 

### HELPER FUNCTIONS ###
def ImportMetadata(mdl):
    meta_file = open(mdl + '\ModelMetadata.json', 'r')
    meta_dict = json.load(meta_file)
    return meta_dict