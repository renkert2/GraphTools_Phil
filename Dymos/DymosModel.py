# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 16:28:42 2021

@author: renkert2
"""

import numpy as np
import openmdao.api as om

class DymosModel(om.ExplicitComponent()):

    def initialize(self):
        self.options.declare('num_nodes', types=int)
        
    def setup(self):
        nn = self.options['num_nodes']
        
        ### DYNAMIC STATES ###
        for i in range(9):
            self.add_input(name=('x'+str(i)), shape=(nn,))
        
        ### INPUTS ###

        ### DISTURBANCES ###
        
        ### PARAMETERS ###
        
        
    def compute(self, inputs, outputs):
        
    
    def compute_partials(self, inputs, partials):
        