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
import os
import sys

def DymosPhase(mdl, tx, include_disturbances = False, path = '', 
               model_opts = {}, state_opts = {}, control_opts = {}, disturbance_opts = {}, parameter_opts = {}):
    
    # Instantiate a Dymos Trajectory
    meta = ImportMetadata(mdl)
    
    # Instantiate a Phase and add it to the Trajectory.
    # Here the transcription is necessary but not particularly relevant.
    model_opts = {"Model":mdl, "Path":path, "include_disturbances":include_disturbances}
    phase = dm.Phase(ode_class=DymosModel, ode_init_kwargs=model_opts, transcription=tx)
    
    # Tell Dymos the states to be propagated using the given ODE.
    state_vars = [x["StateVariable"] for x in meta["StateTable"]]
    for var in state_vars:
        phase.add_state(var, rate_source=var+"_dot", **state_opts)
        
    # Tell Dymos the inputs to be propagated using the given ODE.
    input_vars = [x["InputVariable"] for x in meta["InputTable"]]
    for var in input_vars:
        phase.add_control(var, **control_opts)
        
    if include_disturbances:
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
        self.options.declare('include_disturbances', types=bool, default = False) # "Model" property used to store name of folder containing Model .py functions and variable tables
        self.options.declare('Path', types=str, default = '') # Path to directory containing model folder "Model", i.e. abs_path = Path/Model/
    
    def setup(self):
        nn = self.options['num_nodes']
        mdl = self.options['Model']
        
        self.ImportModel(mdl, self.options["Path"])
        shared_opts = {"num_nodes":nn, "_Metadata":self.Metadata, "include_disturbances":self.options["include_disturbances"]}
        # Instantiate Dynamic Model Subsystem (f), promote x,u,d,theta
        f_subsys = _CalcF(_Calc = self.Calc["f"], _CalcJ = self.CalcJ["f"], **shared_opts)
        self.add_subsystem(name="CalcF", subsys=f_subsys,
                           promotes = ["*"])
        
        # Instantiate Output Model Subsystem (g), promote x,u,d,theta
        g_subsys = _CalcG(_Calc = self.Calc["g"], _CalcJ = self.CalcJ["g"], **shared_opts)
        self.add_subsystem(name="CalcG", subsys=g_subsys,
                           promotes = ["*"])
        
    def ImportModel(self,mdl, mdl_path=''):
        # sys.path.insert(0, mdl)
        if mdl_path:
            sys.path.append(mdl_path)
            
        meta_dict = ImportMetadata(mdl)
        Nx = meta_dict["Nx"]        
        Nu = meta_dict["Nu"]
        Nd = meta_dict["Nd"]
        Ny = meta_dict["Ny"]
        Ntheta = meta_dict["Ntheta"]
        JacStruct = meta_dict["JacStruct"]
        
        # Functions
        func_list = ['f', 'g']
        calc = {}
        for func in func_list:
            mod = impL.import_module(mdl + ".Model_" + func)
            handle = getattr(mod, "Calc_" + func)
            calc[func] = handle
            
        # Jacobian
        # Append function handle to Jacobian metadata
        for f in JacStruct:
            for v in JacStruct[f]:
                if JacStruct[f][v]["NCalc"] > 0:
                    partial = JacStruct[f][v]["Name"]
                    mod = impL.import_module(mdl + ".ModelJ_" + partial)
                    handle = getattr(mod, "CalcJ_" + partial)
                    JacStruct[f][v]["Handle"] = handle
                
        self.Metadata = meta_dict
        self.Calc = calc
        self.CalcJ = JacStruct
        
        if mdl_path:
            sys.path.remove(mdl_path)

class _CalcF(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int) # Number of nodes property, required for Dymos models
        self.options.declare('_Metadata')
        self.options.declare('_Calc') # Number of nodes property, required for Dymos models
        self.options.declare('_CalcJ') # Number of nodes property, required for Dymos models
        self.options.declare('include_disturbances', types=bool, default = False) # "Model" property used to store name of folder containing Model .py functions and variable tables
        
    def setup(self):
        nn = self.options['num_nodes']
        meta = self.options["_Metadata"]
        calcJ = self.options["_CalcJ"]
        
        self.VarNames = {}
        
        ### INPUTS ###
        # x
        state_table = meta["StateTable"]
        self.VarNames["x"] = [x["StateVariable"] for x in state_table]
        for var in state_table:
            self.add_input(var["StateVariable"], desc=var["Description"], shape=(nn,))
        # u
        input_table = meta["InputTable"]
        self.VarNames["u"] = [x["InputVariable"] for x in input_table]
        for var in input_table:
            self.add_input(var["InputVariable"], desc=var["Description"], shape=(nn,))
        # d
        disturbance_table = meta["DisturbanceTable"]
        self.VarNames["d"] = [x["DisturbanceVariable"] for x in disturbance_table]
        if self.options["include_disturbances"]:
            for var in disturbance_table:
                self.add_input(var["DisturbanceVariable"], desc=var["Description"], shape=(nn,))
            
        # theta
        param_table = meta["ParamTable"]
        self.VarNames["theta"] = [x["SymID"] for x in param_table]
        for param in param_table:
            self.add_input(param["SymID"], val=param["Value"], desc=var["Description"], tags=['dymos.static_target'])
        
        ### OUTPUTS ###
        #xdot
        self.VarNames["x_dot"] = [x["StateVariable"]+"_dot" for x in state_table]
        for var in state_table:
            self.add_output(var["StateVariable"]+"_dot", desc="Rate of Change: "+var["Description"], shape=(nn,), units="1.0/s")
 
        ### PARTIALS ###
        # Come Back to this - just get the model working
        arange = np.arange(nn)
        c = np.zeros(nn)
        var_list = ["x", "u", "d", "theta"] if self.options["include_disturbances"] else ["x", "u", "theta"]
        self.VarList = var_list
        for v in var_list:
            J = calcJ[v] # Get the Jacobian information of f w.r.t v
            
            # Calculated Derivatives
            if J["NCalc"]:
                rCalc = np.array(J["rCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                cCalc = np.array(J["cCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NCalc"]):
                    of = self.VarNames["x_dot"][rCalc[i]]
                    wrt = self.VarNames[v][cCalc[i]]
                    rows = arange
                    cols = (c if v == "theta" else arange)
                    self.declare_partials(of=of, wrt=wrt, rows=rows, cols=cols)
            
            # Constant Derivatives
            if J["NConst"]:
                rConst = np.array(J["rConst"]) - 1 # Convert from 1 indexing to 0 indexing
                cConst = np.array(J["cConst"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NConst"]):
                    of = self.VarNames["x_dot"][rConst[i]]
                    wrt = self.VarNames[v][cConst[i]]
                    rows = arange
                    cols = (c if v == "theta" else arange)
                    val = J["valConst"][i]*np.ones(nn)
                    self.declare_partials(of=of, wrt=wrt, val=val, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        nn = self.options['num_nodes']
        arg_list = AssembleVectors(self, inputs, nn)
           
        # Call Calc function
        f = self.options["_Calc"]
        X_dot = f(*arg_list) # Tuple of outputs, X_dot[i][j] corresponds to jth time step of ith output
        
        # Assign to Outputs
        for i,x_dot in enumerate(self.VarNames["x_dot"]): 
            outputs[x_dot] = X_dot[i]
    
    def compute_partials(self, inputs, partials):
        # Assemble Input Variables
        nn = self.options['num_nodes']
        calcJ = self.options["_CalcJ"]
        
        # Assemble Vectors
        arg_list = AssembleVectors(self, inputs, nn)
        
        # Calculate and assign computed derivatives to Partials
        for v in self.VarList:
            J = calcJ[v]
            if J["NCalc"]:
                J_out = J["Handle"](*arg_list)
                rCalc = np.array(J["rCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                cCalc = np.array(J["cCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NCalc"]):
                    of = self.VarNames["x_dot"][rCalc[i]]
                    wrt = self.VarNames[v][cCalc[i]]
                    partials[of, wrt] = J_out[i]


class _CalcG(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int) # Number of nodes property, required for Dymos models
        self.options.declare('_Metadata')
        self.options.declare('_Calc') # Number of nodes property, required for Dymos models
        self.options.declare('_CalcJ') # Number of nodes property, required for Dymos models
        self.options.declare('include_disturbances', types=bool, default = False) # "Model" property used to store name of folder containing Model .py functions and variable tables
        
    def setup(self):
        nn = self.options['num_nodes']
        meta = self.options["_Metadata"]
        calcJ = self.options["_CalcJ"]
        
        self.VarNames = {}
        
        ### INPUTS ###
        # x
        state_table = meta["StateTable"]
        self.VarNames["x"] = [x["StateVariable"] for x in state_table]
        for var in state_table:
            self.add_input(var["StateVariable"], desc=var["Description"], shape=(nn,))
        # u
        input_table = meta["InputTable"]
        self.VarNames["u"] = [x["InputVariable"] for x in input_table]
        for var in input_table:
            self.add_input(var["InputVariable"], desc=var["Description"], shape=(nn,))
        # d
        disturbance_table = meta["DisturbanceTable"]
        self.VarNames["d"] = [x["DisturbanceVariable"] for x in disturbance_table]
        if self.options["include_disturbances"]:
            for var in disturbance_table:
                self.add_input(var["DisturbanceVariable"], desc=var["Description"], shape=(nn,))
            
        # theta
        param_table = meta["ParamTable"]
        self.VarNames["theta"] = [x["SymID"] for x in param_table]
        for param in param_table:
            self.add_input(param["SymID"], val=param["Value"], desc=var["Description"], tags=['dymos.static_target'])
        
        ### OUTPUTS ###
        #y
        output_table = meta["OutputTable"]
        self.VarNames["y"] = [x["OutputVariable"] for x in output_table]
        for var in output_table:
            self.add_output(var["OutputVariable"], desc=var["Description"], shape=(nn,))
            
        ### PARTIALS ###
        # Come Back to this - just get the model working
        arange = np.arange(nn)
        c = np.zeros(nn)
        var_list = ["x", "u", "d", "theta"] if self.options["include_disturbances"] else ["x", "u", "theta"]
        self.VarList = var_list
        for v in var_list:
            J = calcJ[v] # Get the Jacobian information of f w.r.t v
            
            # Calculated Derivatives
            if J["NCalc"]:
                rCalc = np.array(J["rCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                cCalc = np.array(J["cCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NCalc"]):
                    of = self.VarNames["y"][rCalc[i]]
                    wrt = self.VarNames[v][cCalc[i]]
                    rows = arange
                    cols = (c if v == "theta" else arange)
                    self.declare_partials(of=of, wrt=wrt, rows=rows, cols=cols)
            
            # Constant Derivatives
            if J["NConst"]:
                rConst = np.array(J["rConst"]) - 1 # Convert from 1 indexing to 0 indexing
                cConst = np.array(J["cConst"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NConst"]):
                    of = self.VarNames["y"][rConst[i]]
                    wrt = self.VarNames[v][cConst[i]]
                    rows = arange
                    cols = (c if v == "theta" else arange)
                    val = J["valConst"][i]*np.ones(nn)
                    self.declare_partials(of=of, wrt=wrt, val=val, rows=rows, cols=cols)
    
    def compute(self, inputs, outputs):
        nn = self.options['num_nodes']
        arg_list = AssembleVectors(self, inputs, nn)
           
        # Call Calc function
        g = self.options["_Calc"]
        Y = g(*arg_list) # Tuple of outputs, X_dot[i][j] corresponds to jth time step of ith output
        
        # Assign to Outputs
        for i,y in enumerate(self.VarNames["y"]): 
            outputs[y] = Y[i]
            
    def compute_partials(self, inputs, partials):
        # Assemble Input Variables
        nn = self.options['num_nodes']
        calcJ = self.options["_CalcJ"]
        
        # Assemble Vectors
        arg_list = AssembleVectors(self, inputs, nn)
        
        # Calculate and assign computed derivatives to Partials
        for v in self.VarList:
            J = calcJ[v]
            if J["NCalc"]:
                J_out = J["Handle"](*arg_list)
                rCalc = np.array(J["rCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                cCalc = np.array(J["cCalc"]) - 1 # Convert from 1 indexing to 0 indexing
                for i in range(J["NCalc"]):
                    of = self.VarNames["y"][rCalc[i]]
                    wrt = self.VarNames[v][cCalc[i]]
                    partials[of, wrt] = J_out[i]

### HELPER FUNCTIONS ###
def ImportMetadata(mdl):
    meta_file = open(mdl + '\ModelMetadata.json', 'r')
    meta_dict = json.load(meta_file)
    return meta_dict

def AssembleVectors(obj, inputs, nn):
    arg_list = []
    for v in ["x", "u", "d", "theta"]:
        arg_list_inner = []
        for x in obj.VarNames[v]:
            if v != "d":
                arg_list_inner.append(inputs[x])
            else: 
                if obj.options["include_disturbances"]:
                    arg_list_inner.append(inputs[x])
                else: 
                    arg_list_inner.append(np.zeros(nn))
        arg_list.append(arg_list_inner)
    return arg_list