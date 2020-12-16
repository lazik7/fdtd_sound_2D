# -*- coding: utf-8 -*-
# python > 3.6


#---------------------------------------------------------- damping_factor

# https://en.wikipedia.org/wiki/Attenuation
# http://www.kayelaby.npl.co.uk/general_physics/2_4/2_4_1.html
# http://www.kayelaby.npl.co.uk/general_physics/2_4/2_4_6.html
air_damping_factor = 1.64 # 20 °C, dB/(MHz*cm)
water_damping_factor = 0.0022 #


#----------------------------------------------------------

class liquid():

    def wavelength(self, f):
        return self.c/f

    def mass_stream_for_heat(self, P, dT):
        """
        P - powet, W
        dt - temperature difference, C

        return mass stream, kg/s
        """
        return P/(self.cp*dT) #  strumien masy, kg/s


    def mass_stream_for_heat__kg_per_min(self, P, dT):
        """
        P - powet, W
        dt - temperature difference, C

        return mass stream, kg/s
        """
        return P/(self.cp*dT)*60 #  strumien masy, kg/min


    def mass_stream_for_heat__kg_per_h(self, P, dT):
        """
        P - powet, W
        dt - temperature difference, C

        return mass stream, kg/s
        """
        return P/(self.cp*dT)*3600 #  strumien masy, kg/s

    


###########################################################################################
#                functions
###########################################################################################


#---------------------------------------------------------------------------------- ocean water


def ocean_water_density(S):
    """
    Acustical oceanography principies and applications
    Clarence S. Clay
    Herman Medwin
    str 90

    S - salinity, ppt
    """
    return (1. + S*1e-3)*1e3 # kg/m3


def gauge_pressure_due_to_water_column(q,z,g=9.81):
    """
    Acustical oceanography principies and applications
    Clarence S. Clay
    Herman Medwin
    str 90

    q - density, kg/m3
    z - depth, m    
    
    return P - gauge pressure due to water column, Pa
    """

    return q*g*z # Pa


def SpeedOfSound_OceanWater(T=20., S=0., P=0):
    """
    Acustical oceanography principies and applications
    Clarence S. Clay
    Herman Medwin
    str 88

    T - temperature, *C
    S - salinity, ppt
    P - gauge pressure due to water column, Pa
    
    return c, sound speed, m/s
    """
    return 1449.2 + 4.6*T - 0.055*T**2 + 0.00029*T**3 + (1.34-0.010*T)*(S-35.0) + 1.58e-6*P # m/s



        
                                                                                                                                             
#---------------------------------------------------------------------------------- air

def SpeedOfSound_air(T):
        """
        T, C
        """
        
        return 331.5*(1+T/273.15)**0.5 # m/s

def density_air(T):
	"""
        T, [C]

        return
        density, [kg/m3]

        approximation range: -50 - 1200 C
	"""
	a6 = 3.51339e-18	     #  3.356e-19    (9.553%)
	a5 = -1.49119e-14	     #  1.176e-15    (7.887%)
	a4 = 2.55675e-11	     #  1.547e-12    (6.052%)
	a3 = -2.31452e-08	     #  9.456e-10    (4.085%)
	a2 = 1.24366e-05	     #  2.683e-07    (2.157%)
	a1 = -0.00447181	     #  3.088e-05    (0.6906%)
	a0 = 1.28977		     #  0.001049     (0.08134%)

	return a6*T**6+a5*T**5+a4*T**4+a3*T**3+a2*T**2+a1*T+a0


class Air(liquid):
    def __init__(self, T=20.):
        """
        T - temperature, *C
        """

        self.calculate_the_parameters(T)

    def calculate_the_parameters(self, T=20):
        """
        T - temperature, *C
        """
        
        self.T = T
        self.ρ = density_air(T)
        self.c = SpeedOfSound_air(T)
        self.Z = self.ρ*self.c

    def calculate_Z(self, T=20.):
        return density_air(T) * SpeedOfSound_air(T)


air = Air()



def density_water(T):
# T - temperatura, C
	return 999.732+0.07935*T-0.00857*T**2+5.83e-5*T**3-2.677e-7*T**4+4.843e-10*T**5 # kg/m3

def SpeedOfSound_water(T):
	return 1402.736 + 5.03358*T-0.0579506*T**2+3.31636e-4*T**3-1.45262e-6*T**4+3.0449e-9*T**5


class Water(liquid):
    def __init__(self, T=20.):
        """
        T - temperature, *C
        """

        self.calculate_the_parameters(T)

    def calculate_the_parameters(self, T=20):
        """
        T - temperature, *C
        """
        
        self.T = T
        self.ρ = density_water(T)
        self.c = SpeedOfSound_water(T)
        self.Z = self.ρ*self.c
        self.cp = self.heat_capacity(T)

    def calculate_Z(self, T=20.):
        return density_water(T) * SpeedOfSound_water(T)


    def heat_capacity(self, T):
        """
        T - temperature, *C
        """
        return (4.214-0.00220*T+4.21e-5*T**2-2.817e-7*T**3+8.4525e-10*T**4)*1000. # J/(kg*K)


    def heat(self, T_start, T_stop, m):
        T_average = (T_start + T_stop)/2.
        dT = T_stop-T_start
        cp = self.heat_capacity(T_average)
        
        return cp*dT*m # J

    def time_to_heat(self, T_start, T_stop, m, P):
        """

        return time, s
        """
                
        W = self.heat(T_start, T_stop, m)

        return W/P # s

    def heat_power(self, T_start, T_stop, m, t):
        """
        T_start - temp on start, *C
        T_stop - temp stop, *C
        m - mass, kg
        t - time, s
        """

        return self.heat(T_start,T_stop,m)/t

        
        
        


water = Water()
        
