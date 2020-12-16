#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# python > 3.6

#from math import pi
import numpy as np
#from matplotlib import pyplot


###########################################################################################
#                functions
###########################################################################################


#---------------------------------------------------------------------------------- 

def two_port_metworks_type_T(Z1, Z2, Z3):

        #[[A, B],[C, D]]

        #------------------------------- A
        A = Z1/Z3 + 1.
        #------------------------------- B
        B = -1.*( (Z3+Z1)*(-1.*(Z2+Z3)/Z3) + Z3)  
        #------------------------------- C
        C = 1./Z3
        #------------------------------- D
        D = (Z2+Z3)/Z3
        #-------------------------------
        
        return np.matrix([[A,B],[C,D]])



def calculation_parallel_R(R1, R2):
        return 1./(1./R1 + 1./R2)


def calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1, Z2, Z3):
        return Z1 + calculation_parallel_R(Z2, Z3)
        
                                                                                                                                             
#----------------------------------------------------------------------------------


if __name__ == "__main__":

        z = np.around( two_port_metworks_type_T(2., 4., 10.), decimals=2)

        assert np.all( z == np.around( np.matrix([[1.2,6.8],[0.1,1.4]]), decimals=2) )
        
