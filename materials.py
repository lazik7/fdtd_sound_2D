# -*- coding: utf-8 -*-
# python > 3.6

from math import pi
import numpy as np
from matplotlib import pyplot

import scipy.optimize as so
from scipy.special import jv, j0, j1, y0, y1 #, jn, yn

import electrical
import liquids



material_list = []


class Material():
        def __init__(self, name="", ρ=0, c=0, λ=0, E=0, ν=0, damping_factor_dB_na_mm=0, damping_coefficient=0, cp=0, α=0):
                self.name = name
                self.ρ = ρ # density, kg/m3
                self.c = c # sound speed, m/s
                self.λ = λ # Thermal conductivity, W/(m*K)
                self.E = E # Young, Pa
                self.ν = ν # Poisson
                self.damping_factor_dB_na_mm = damping_factor_dB_na_mm
                self.damping_coefficient = damping_coefficient
                self.cp = cp # Heat capacity, J/kgK
                self.α = α # thermal expansion, 1/K
                self.useMori = False

                if c == 0 and E != 0:
                        self.sound_speed()

                if c != 0 and E == 0 and ρ!=0:
                        self.young()
                        
                self.impedance()


        #---------------------------------------- temp

        def time_to_have_temp(self, Tstart, Tstop, m, P):
                Q = self.cp*m*(Tstop-Tstart)
                return Q/P
                
        #----------------------------------------

        #---------------------------------------- sound

        def time_echo(self, L):
                return 2*L/self.c

        #----------------------------------------
                

        def impedance(self):
                self.Z = self.ρ * self.c

        def sound_speed(self):
                self.c = (self.E/self.ρ)**0.5

        def sound_speed_poisson(self):
                self.c = (self.E/(self.ρ*(1.-self.ν**2)))**0.5

        def sound_speed_shaft(self, l, D, k=1.):
                r = D/2.
                self.c = (self.E/self.ρ)**0.5 / (1. + (k**2 * self.ν**2 * pi**2 * r**2 / (2.*l**2) ) )**0.5

        def sound_speed_shaft_Mori(self, f=21e3, D=50e-3):

                def myfunc(o):
                        return lambda z: z*jv(0,z)-(1-o)*jv(1,z)

                def bessel_equation(o):
                        alpha = 1.84+0.68*o  #approximate solution
                        z = np.arange(0., np.pi*2, 0.001)
                        equation = myfunc(o)
                        beta = so.bisect(equation, alpha-1., alpha+1., maxiter=1000000)
                        return beta
                
                #c = (self.E/(self.ρ*(1-self.ν**2)))**0.5
                c = (self.E/self.ρ)**0.5
                x = bessel_equation(self.ν)
                k_0 = 2*np.pi*f/c
                self.c = c*np.sqrt((1 - (1-self.ν**2)*(k_0*D/(2*x))**2)/(1 - (1-3*self.ν**2-2*self.ν**3)*(k_0*D/(2*x))**2))
                return self.c


        def sound_speed_lamba(self):
                """
                Awram Levi - „Ultradźwiękowe Badania Nieniszczące Własności Mechanicznych Cienkich Elementów Konstrukcyjnych”. - Prace IPPT – IFTR Reports 1/2010, Warszawa 2010, s.104

                E - nmodul Younga, Pa, domyslnie dla brazu
                ρ - gestosc, kg/m3, domyslnie dla brazu
                ν=0.33 - liczba Poissona, domyslnie dla brazu
                """

                return ( self.E/( self.ρ*(1.-self.ν**2) ) )**0.5
                

        def young(self):
                self.E = self.ρ*self.c**2


        ###########################################################################################
        # resonances
        ###########################################################################################


        def sonotroda_schodkowa(self, f=np.arange(15e3,50e3,1), D1=32.02e-3, D2=30.04e-3, l1=30.55e-3, l2=35.25e-3, useMori=False):

                
                #self.c = self.sound_speed_shaft_Mori(f, D=(D2+D1)/2)
                
                f, Z1, Z2, Z3 = self.Z1_Z2_Z3_shaft_Two_port_network(f=f, D=D1, d=0, l=l1, useMori=False)
                Z_c1 = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1, Z2, Z3)

                f, Z1, Z2, Z3 = self.Z1_Z2_Z3_shaft_Two_port_network(f=f, D=D2, d=0, l=l2, useMori=False)
                Z = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1+Z_c1, Z2, Z3)

                
                #pyplot.semilogy(f*1e-3, np.absolute(Z_c1))
                pyplot.semilogy(f*1e-3, np.absolute(Z))
                pyplot.show()


        #-------------- shaft

        def Z_acustic_shaft_one_frequency(self, f=1e6, D=51e-3, d=0, l=121e-3):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                D = srednica, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                l - dlugosc, m
                """

                q=self.ρ
                vz=self.c

                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2 # m2, A1=A2=A

                if self.useMori:
                        vz = self.sound_speed_shaft_Mori(f, D)
                        
                h = vz/f
                kz = 2.*np.pi/h
        
                return 1j*A*q*vz*np.tan(kz*l)
        

        def Z_shaft(self, f=np.arange(15e3,60e3,1), D=51e-3, d=0, l=121e-3, useMori=True):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                f - tablica czestotliwosci, Hz
                D - srednica, m
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                self.useMori = useMori
                        
                impedancja = np.vectorize( self.Z_acustic_shaft_one_frequency, excluded=['D', 'd', 'l'])

                Z = impedancja(f, D, d, l)
                #O = np.degrees(np.angle(Z))
                #ZdB = np.real(Zul_dB(Z)) # 20.*np.log10(zul/50.+1.)

                #Z = ((10.**ZdB/20.)-1)*50.

                return f, Z

        def Z_shaft_plot(self, f=np.arange(15e3,66e3,1), D=51e-3, d=0, l=121e-3, useMori=True):
                f, Z = self.Z_shaft(f, D, d, l, useMori)
                print(f[np.log(Z).argmin()])
                pyplot.semilogy(f*1e-3, np.absolute(Z))
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("|Z|")
                pyplot.show()

                return f, Z



        def __Z1_Z2_Z3_acustic_shaft_Two_port_network(self, f=1e6, D=51e-3, d=0, l=121e-3, q=2810., vz=7100):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 69

                f - czestotliwosc, Hz
                D - srednica, m
                d - srednica, ms
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                if self.useMori:
                        vz = self.sound_speed_shaft_Mori(f, D)
                else:
                        vz = self.c

                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2 # m2, A1=A2=A

                h = vz/f
                kz = 2.*np.pi/h

                Zc = q*vz

                Z1 = 1j * A * Zc * np.tan(0.5*kz*l)
                Z2 = Z1
                Z3 = (A*Zc/1j) * 1./(np.sin(kz*l))

                return Z1, Z2, Z3


        def Z1_Z2_Z3_shaft_Two_port_network(self, f=np.arange(15e3,66e3,1), D=51e-3, d=0, l=121e-3, useMori=True, use_speed_of_sound_approximation=False):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                f - tablica czestotliwosci, Hz
                D - srednica zewnetrzna, m
                d - srednica otworu, mm
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                if use_speed_of_sound_approximation:
                        self.c = self.speed_of_sound_approximation(l)

                
                self.useMori = useMori
        
                impedancja = np.vectorize( self.__Z1_Z2_Z3_acustic_shaft_Two_port_network, excluded=['D', 'd', 'l', 'q', 'vz'])

                Z1, Z2, Z3 = impedancja(f, D, d, l, self.ρ, self.c)

                return f, Z1, Z2, Z3

        def Z_shaft_Two_port_network(self, f=np.arange(15e3,66e3,1), D=51e-3, d=0, l=121e-3, useMori=True, use_speed_of_sound_approximation=False, Z1_load=0, Z2_load=0):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                f - tablica czestotliwosci, Hz
                D - srednica zewnetrzna, m
                d - srednica otworu, mm
                l - dlugosc, m
                vz - prędkość fazowa fali osiowej, m/s
                """
                
                self.useMori = useMori

                f, Z1, Z2, Z3 = self.Z1_Z2_Z3_shaft_Two_port_network(f, D, d, l, useMori, use_speed_of_sound_approximation)

                A = np.pi*(D/2.)**2
                
                Z = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1, Z2, Z3)
                
                return f, Z, Z1, Z2, Z3


        def Z_shaft_Two_port_network_plot(self, f=np.arange(15e3,80e3,1), D=51e-3, d=0, l=121e-3, useMori=True, use_speed_of_sound_approximation=False, Z1_load=0, Z2_load=0):
                
                f, Z, Z1, Z2, Z3 = self.Z_shaft_Two_port_network(f, D, d, l, useMori, use_speed_of_sound_approximation, Z1_load, Z2_load)
                print("%.3f kHz" % (1e-3*f[np.log(Z).argmin()]))
                pyplot.cla()
                pyplot.semilogy(f*1e-3, np.absolute(Z1), label="Z1")
                #pyplot.semilogy(f*1e-3, np.absolute(Z2))
                pyplot.semilogy(f*1e-3, np.absolute(Z3), label="Z3")
                pyplot.semilogy(f*1e-3, np.absolute(Z), color='red', label="Z")
                pyplot.legend()
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("|Z|")
                pyplot.show()

                return f, Z, Z1, Z2, Z3



        def resonant_frequency_shaft(self, D=51e-3, l=121.1e-3, skok=1):
                """
                l - length, m
                D - diameter, m
                """

                def find_the_first_minimum(Z):

                        i = Z.argmax()

                        if i < len(Z):
                                while Z [i] > Z[i+1]:
                                        if i >= len(Z)-2:
                                                return i
                                        i = i + 1

                        return i
                        

                c = (self.E/self.ρ)**0.5
                f = c/(2*l)

                f, Z = duralumin_PA6.Z_shaft(f=np.arange(100,f+1e3,skok), D=D, l=l, useMori=True)

                Z = np.log10(np.absolute(Z))

                # by uniknac znajdowania maksimow tam dzie jest "nan", zamienia sie "nan" na 0
                Z = np.nan_to_num(Z)
                
                i = find_the_first_minimum(Z)

                # moze sie znalesc odnalezienie minimum na kocu danych
                while i >= len(Z)-2:
                        Z = Z[:int(-10*skok)]
                        i = find_the_first_minimum(Z)
                                                
                
                return f[i]
                #return 0
                

        #def shatf_find_lenth_for_freq(self, f, D=51e-3):
                #"""
                #f - freq to find
                #"""

                #pass


        def Zac2_shaft(self, f=20e3, D=51e-3, d=0e-3, l=100e-3):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                f - czestotliwosc, Hz
                D - srednica zewnetrzna, m
                d - srednica otworu, mm
                l - dlugosc, m
                """

                A = np.pi*(D/2)**2 - np.pi*(d/2)**2

                w = 2*np.pi*f
                kz = w/self.c

                return 1j*A*self.ρ*self.c*np.tan(kz*l)


        def fr_shaft(self, l=100e-3, n=1):
                """
                n = 1,2,3...
                """
                return n*self.c/(2.*l)

        def fa_shaft(self, l=100e-3, m=1):
                """
                m = 1,3,5...
                """
                return m*self.c/(4.*l)               


        #-------------- cone

        def Z_acustic_cone_one_frequency(self, f=1e6, D1=30e-3, D2=20e-3, l=65e-3, d=0):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                D1 = srednica, m
                D2 = srednica, m
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                
                q=self.ρ
                vz=self.c

                r1 = D1/2.
                r2 = D2/2.
                y = -l*r1/(r1-r2)

                A1 = np.pi*(D1/2.)**2 - np.pi*(d/2.)**2 # m2

                if self.useMori:
                        D = np.absolute(D1-D2)/2**0.5
                        vz = self.sound_speed_shaft_Mori(f, D1+D)

                h = vz/f
                kz = 2.*np.pi/h

                Z = ( 1j*A1*q*vz*(l+y)/(y**2) ) * ( (kz*y*(l+y)+1./kz)*np.tan(kz*l)-l ) / (np.tan(kz*l) + kz*y)
        
                return Z

        def Z_cone(self, f=np.arange(15e3,55e3,1), D1=51e-3, D2=40e-3, l=46.9e-3, useMori=True):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                D1 = srednica, m
                D2 = srednica, m
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                self.useMori = useMori
        
                impedancja_cone = np.vectorize( self.Z_acustic_cone_one_frequency, excluded=['D1', 'D2', 'l', 'd'])

                Z = impedancja_cone(f, D1, D2, l)
                #O = np.degrees(np.angle(Z))
                #ZdB = np.real(Zul_dB(Z)) # 20.*np.log10(zul/50.+1.)

                #Z = ((10.**ZdB/20.)-1)*50.

                return f, Z


        def Z_cone_plot(self, f=np.arange(15e3,55e3,1), D1=51e-3, D2=40e-3, l=46.9e-3, useMori=True):
                f, Z = self.Z_cone(f, D1, D2, l, useMori=useMori)
                pyplot.semilogy(f, np.absolute(Z))
                pyplot.show()


        def Z_acustic_cone_Two_port_networke_one_frequency(self, f=1e6, D1=30e-3, D2=20e-3, l=65e-3, d=0):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 69
                """

                vz=self.c

                if self.useMori:
                        #D = np.absolute(D1-D2)/2**0.5
                        #vz = self.sound_speed_shaft_Mori(f, D1+D)
                        vz = self.sound_speed_shaft_Mori(f, D1)
                        

                A1 = np.pi*(D1/2.)**2 - np.pi*(d/2.)**2 # m2
                #A2 = np.pi*(D2/2.)**2 - np.pi*(d/2.)**2 # m2

                r1 = D1/2.
                r2 = D2/2.
                #r1 = (A1/np.pi)**0.5
                #r2 = (A2/np.pi)**0.5
                y = -l*r1/(r1-r2)

                

                h = vz/f
                kz = 2.*np.pi/h

                Zc = self.ρ*vz

                Z0 = A1*Zc/(1j)
                y0 = (l+y)/y
                tankzl = np.tan(kz*l)
                sinkzl = np.sin(kz*l)

                Z1 = Z0 * ( 1./tankzl + 1./(kz*y) - y0/sinkzl )
                Z2 = Z0 * y0 * ( y0*( 1./tankzl - 1./(kz*y) ) - 1./sinkzl )
                Z3 = Z0 * y0/sinkzl

                return Z1, Z2, Z3

        def Z_cone_Two_port_networke(self, f=np.arange(15e3,55e3,1), D1=51e-3, D2=40e-3, d=0, l=46.9e-3, useMori=True, Z1_load=0, Z2_load=0):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                D1 = srednica, m
                D2 = srednica, m
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                self.useMori = useMori
        
                impedancja_cone = np.vectorize( self.Z_acustic_cone_Two_port_networke_one_frequency, excluded=['D1', 'D2', 'l', 'd'])

                Z1, Z2, Z3 = impedancja_cone(f, D1, D2, l, d)
                #O = np.degrees(np.angle(Z))
                #ZdB = np.real(Zul_dB(Z)) # 20.*np.log10(zul/50.+1.)

                #Z = ((10.**ZdB/20.)-1)*50.

                #A1 = np.pi*(D1/2.)**2
                #A2 = np.pi*(D2/2.)**2
                
                #Z = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1+A1*Z1_load, Z2+A2*Z2_load, Z3)
                Z = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1, Z2, Z3)

                return f, Z, Z1, Z2, Z3
                

        def Z_cone_Two_port_networke_plot(self, f=np.arange(15e3,55e3,1), D1=51e-3, D2=40e-3, l=46.9e-3, useMori=True, Z1_load=0, Z2_load=0):
                f, Z, Z1, Z2, Z3 = self.Z_cone_Two_port_networke(f, D1, D2, l, useMori=useMori, Z1_load=Z1_load, Z2_load=Z2_load)
                pyplot.semilogy(f, np.absolute(Z))
                pyplot.show()


        #-------------- shaft plus cone


        def Z_shaft_plus_cone_Two_port_networke(self, f=np.arange(15e3,50e3,1), D=51e-3, d=0, l=20e-3, D1=51e-3, D2=69.7e-3,L=47.1e-3, useMori=True, use_speed_of_sound_approximation=False, Z1_load=0, Z2_load=0):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 67

                f - tablica czestotliwosci, Hz
                D - srednica zewnetrzna, m
                d - srednica otworu, mm
                l - dlugosc, m
                q - gestosc, kg/m3
                vz - prędkość fazowa fali osiowej, m/s
                """

                self.useMori = useMori

                f, Z1, Z2, Z3 = self.Z1_Z2_Z3_shaft_Two_port_network(f, D, d, l, useMori, use_speed_of_sound_approximation)

                #A1 = np.pi*(D1/2.)**2
                A = np.pi*(D2/2.)**2

                f_c, Z_c, Z1_c, Z2_c, Z3_c = self.Z_cone_Two_port_networke(f, D1, D2, L-l, useMori=useMori, Z1_load=Z1_load, Z2_load=0)
                
                Z = electrical.calculation_Z_two_port_metworks_type_T_closed_on_the_right(Z1+Z_c, Z2+A*Z2_load, Z3)

                return f, Z, Z1, Z2, Z3


        def Z_shaft_plus_cone_Two_port_network_plot(self, f=np.arange(15e3,50e3,1), D=51e-3, d=12e-3, l=20e-3, D1=51e-3, D2=70e-3, L=47.1e-3, useMori=True, use_speed_of_sound_approximation=False, Z1_load=0, Z2_load=0):
                
                f, Z, Z1, Z2, Z3 = self.Z_shaft_plus_cone_Two_port_networke(f, D, d, l, D1, D2, L, useMori, use_speed_of_sound_approximation, Z1_load, Z2_load)
                print("%.3f kHz" % (1e-3*f[np.log(Z).argmin()]))
                pyplot.cla()
                #pyplot.semilogy(f*1e-3, np.absolute(Z1), label="Z1")
                #pyplot.semilogy(f*1e-3, np.absolute(Z2))
                #pyplot.semilogy(f*1e-3, np.absolute(Z3), label="Z3")
                pyplot.semilogy(f*1e-3, np.absolute(Z), color='red', label="Z")
                pyplot.legend()
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("|Z|")
                pyplot.show()

                return f, Z



        #-------------- Cylindrical shell

        def cylindrical_shell__radial_vibration_with_axial_symetry(self, D=30e-3, h=3e-3, m=0, L=300e-3):
                """
                E.P. Mechel pt.”Formulas of Acoustics” str 1189
                see also: Dym (1973)

                self.ρ = ρ # density, kg/m3
                self.E = E # Young, Pa
                self.ν = ν # Poisson

                D - outer diameter, m, default D=30e-3
                h - thickness of the wall, m, default h=3e-3
                m - vibration mod, default, m=0
                L - the length of the pipe, m

                return f, Hz
                """

                R = D/2.

                w = ( self.E/(self.ρ*(1.-self.ν**2)) * ( 1./R**2 + ((h**2)/12.)*(m*np.pi/L)**4 ) )**0.5

                return w/(2.*np.pi)

        def cylindrical_shell__vibration(self, D=30e-3, h=3e-3, m=0, n=0, L=300e-3, x0=(1,)):
                """
                E.P. Mechel pt.”Formulas of Acoustics” str 1189
                see also: Dym (1973)

                self.ρ = ρ # density, kg/m3
                self.E = E # Young, Pa
                self.ν = ν # Poisson

                D - outer diameter, m, default D=30e-3
                h - thickness of the wall, m, default h=3e-3
                m - vibration mod, default, m=0
                L - the length of the pipe, m

                return f, Hz
                """

                R = D/2.

                a = 0.5*(1.-self.ν)
                H = ((h/R)**2)/12.
                Λ = m*np.pi*R/L

                Q0 = (1-self.ν**2)*Λ**4
                Q0 += H*( (Λ**2 + n**2)**4 + (9./4.)*(1.-self.ν**2)*Λ**4 + 4.*Λ**2*n**2 + n**4 + 6.*Λ**4*n**2 - 8*Λ**2*n**4 -2.*n**6 )
                Q0 += H**2 * ( 0.25*n**2 - (3./2.)*self.ν**2*Λ**4*n**2 - 0.5*n**6 + (9./4.)*Λ**8 + 4.*Λ**6*n**2 + (3./2.)*Λ**4*n**4 + 0.25*n**8 )
                Q0 += H**3*(a*(1.-a)*Λ**4*n**4)
                Q0 = a*Q0

                Q1 = (5.-4.*a)*Λ**2 + n**2
                Q1 += H*( (9./4.)*Λ**2 + (1./a + 1./4.)*n**2 - (2./a + 4.*a)*Λ**2*n**2 - (2./a)*n**4 + ((1.+a)/a)*(Λ**2+n**2)**3 )
                Q1 += H**2*( (9./4.)*Λ**6 + (11./4.-a)*Λ**4*n**2 + (3./4.-a)*Λ**2*n**4 + 0.25*n**6 )
                Q1 = a*Q1

                Q2 = (Λ**2+n**2)**2
                Q2 += H*( (9./4.)*Λ**4 + (3./2. + a + 1/a)*Λ**2*n**2 + (5./4.)*n**4 )
                Q2 += 0.25*H**2*n**4
                Q2 = a*Q2

                Q3 = 1. + H*(Λ**2 + n**2)
                
                Q4 = (1.+a)*(Λ**2 + n**2) + H*( (9./4.)*Λ**2 + (1.+a*4.)*n**2 )

                fun = lambda K: (K**6 - (Q3+Q4)*K**4 + (Q1+Q2)*K**1 - Q0)

                sol = so.root(fun, x0, method='hybr')
                                
                w = (sol.x/R)*( self.E/(self.ρ*(1.-self.ν**2)) )**0.5

                return w/(2.*np.pi)



        #---------------------------------------------------------

        
        def Obl_isotropowe_C11(self):
                return self.E*(1.-self.ν)/( (1.+self.ν)*(1-2*self.ν) )

        def Obl_isotropowe_C22(self):
                return self.Obl_isotropowe_C11()

        def Obl_isotropowe_C33(self):
                return self.Obl_isotropowe_C11()

        def Obl_isotropowe_C44(self):
                return self.Obl_isotropowe_C11() * (1.-2*self.ν)/(2.*(1.-self.ν))

        def Obl_isotropowe_C12(self):
                return self.Obl_isotropowe_C11() * (self.ν/(1.-self.ν))

        def Obl_isotropowe_C66(self):
                return self.Obl_isotropowe_C11() * (1.-2*self.ν)/(2.*(1.-self.ν))
        
        
        def Obl_C11_C22_C12_C66(self):
                print( self.Obl_isotropowe_C11() )
                print( self.Obl_isotropowe_C22() )
                print( self.Obl_isotropowe_C12() )
                print( self.Obl_isotropowe_C66() )

        def Obl_C11_C22_C12_C66_GPa(self):
                print( "C11 = ",round(self.Obl_isotropowe_C11()/1.e9, 1), "GPa" )
                print( "C22 = ",round(self.Obl_isotropowe_C22()/1.e9, 1), "GPa" )
                print( "C12 = ",round(self.Obl_isotropowe_C12()/1.e9, 1), "GPa" )
                print( "C66 = ",round(self.Obl_isotropowe_C66()/1.e9, 1), "GPa" )



        # C11 C22 C33 C12 C23 C31 C44 C55 C66
        def Obl_C11_C22_C33_C12_C23_C31_C44_C55_C66_GPa(self):
                print( "C11 = ",round(self.Obl_isotropowe_C11()/1.e9, 1), "GPa" )
                print( "C22 = ",round(self.Obl_isotropowe_C22()/1.e9, 1), "GPa" )
                print( "C33 = ",round(self.Obl_isotropowe_C33()/1.e9, 1), "GPa" )
                print( "C12 = ",round(self.Obl_isotropowe_C12()/1.e9, 1), "GPa" )
                print( "C23 = ",round(self.Obl_isotropowe_C12()/1.e9, 1), "GPa" )
                print( "C31 = ",round(self.Obl_isotropowe_C12()/1.e9, 1), "GPa" )
                print( "C44 = ",round(self.Obl_isotropowe_C44()/1.e9, 1), "GPa" )
                print( "C55 = ",round(self.Obl_isotropowe_C44()/1.e9, 1), "GPa" )
                print( "C66 = ",round(self.Obl_isotropowe_C66()/1.e9, 1), "GPa" )
        
                

        ################################################################################### resonances
                



###########################################################################################
#                physics constants
###########################################################################################

e0 = 8.854187817e-12 # Vacuum permittivity, F/m


# plexi
plexi = Material(
                        name = "plexi",
                        ρ = 1180.,      # density, kg/m3
                        c = 2700.       # sound speed, m/s
                )

material_list.append(plexi)

# brass
brass = Material(
                        name = "brass",
                        E = 113.5e9,     # GPa, Wiki mean
                        ν = 0.33,       # Wiki for copper
                        ρ = 8400.,      # density, kg/m3
                        λ = 110.        # W/(m*K)
                )

material_list.append(brass)


# milk
milk = Material(
                        name = "milk",
                        ρ = 1033.,      # density, kg/m3
                        c = 1630.       # sound speed, m/s
                )




# hard rubber
hard_rubber = Material(
                        name = "hard_rubber",
                        ρ = 1200.,      # density, kg/m3
                        c = 1570.       # sound speed, m/s
                )

material_list.append(hard_rubber)


# polystyrene
 
polystyrene = Material(
                        name = "polystyrene",
                        ρ = 1020.,              # density, kg/m3
                        c = 2300.,              # sound speed, m/s
                        ν = 0.34,               # https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
                        E = (2.28+3.28)*1e9,    # http://www.dielectriccorp.com/downloads/thermoplastics/polystyrene.pdf
                        cp = (1.2e3 + 1.3e3)/2,
                        λ = 0.03        # W/(m*K)
                )

material_list.append(polystyrene)



polystyrene_crystalline = Material(
                        name = "polystyrene_crystalline",
                        ρ = 1020.,              # density, kg/m3
                        c = 2300.,              # sound speed, m/s
                        ν = 0.34,               # https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
                        E = (2.28+3.28)*1e9,    # Pa, http://www.dielectriccorp.com/downloads/thermoplastics/polystyrene.pdf
                        cp = (1.2e3 + 1.3e3)/2,
                        λ = 0.316               # W/(m*K) http://tuprints.ulb.tu-darmstadt.de/2145/1/Algaer_Dissertation.pdf
                )

material_list.append(polystyrene_crystalline)

"""crystalline_polystyrene_c = Material(
                        ρ = 1020.,      # density, kg/m3
                        c = 2259.4        # sound speed, m/s
                )"""



########################
# Superalloy
########################

def Inconel_718_h_E(t):
        """
        t - temperature, *C

        reurn E - Young's modulus
        """

        a0=209.291489971234
        a1=-0.0388105997015295
        a2=-3.04677573841227e-05

        return (a0+a1*t+a2*t**2)*1.e9


# http://www.bibusmetals.pl/fileadmin/editors/countries/bmpl/Data_sheets/Inconel_718_karta_katalogowa.pdf
Inconel_718 = Material(
                        name = "Inconel_718",
                        ρ = 8190,       # density, kg/m3
                        λ = 11.4,        # W/mK
                        E = Inconel_718_h_E(20)
                )

material_list.append(Inconel_718)




########################
# resins
########################

# resins 509FM
resins_509FM = Material(
                        name = "resins_509FM",
                        E = 4.92e9,
                        ν = 0.37,
                        ρ = 1125.,      # density, kg/m3
                        c = 2471.4      # sound speed, m/s
                )

material_list.append(resins_509FM)



# resins Epo_Tek_353ND
resins_Epo_Tek_353ND_epoxy = Material(
                        name = "resins_Epo_Tek_353ND_epoxy",
                        E = 4.92e9,
                        ν = 0.37,
                        ρ = 2600., #ρ = 2800.,      # density, kg/m3
                        c = 2471.4,      # adopted sound speed, m/s
                        damping_factor_dB_na_mm = 31.
                )

material_list.append(resins_Epo_Tek_353ND_epoxy)



# resins Epo_Tek_301
resins_Epo_Tek_301_epoxy = Material(
                        name = "resins_Epo_Tek_301_epoxy",
                        E = 4.92e9,
                        ρ = 2800.,      # adopted, density, kg/m3
                        c = 2700.,      # sound speed, m/s
                        damping_factor_dB_na_mm = 13.2
                )

material_list.append(resins_Epo_Tek_301_epoxy)


# resins LM285
# resins Epo_Tek_301
resins_LM285 = Material(
                        name = "resins_LM285",
                        ρ = 1166.,      # density, kg/m3
                        c = 2750.,      # sound speed, m/s
                )

material_list.append(resins_LM285)


resins_LM285_z_corundum = Material(
                        name = "resins_LM285_z_corundum",
                        ρ = 2130.,      # density, kg/m3
                        c = 3103.,      # sound speed, m/s
                )

material_list.append(resins_LM285_z_corundum)


# microballoons 20%
resins_LM285_z_microballoons_0_2 = Material(
                        name = "resins_LM285_z_microballoons_0_2",
                        ρ = 403.,      # density, kg/m3
                        c = 2101.,     # m/s,  # not pass above 2MHz
                )

material_list.append(resins_LM285_z_microballoons_0_2)


########################
# cast iron
########################


cast_iron = Material(
                        name = "cast_iron",
                        ρ = 7250.,              # density, kg/m3
                        E = 1e9*(80+150)/2.,    # Young module, Pa
                )

material_list.append(cast_iron)



########################
# steel
########################

class Stell(Material):
        def E_carbon_steel_C_03(self,T):
                list_T = [-200,-129,-73,21,93,149,204,260,316,371,427,482,538,593] # *C
                list_E = np.array( [31.4,30.8,30.2,29.5,28.8,28.3,27.7,27.3,26.7,25.5,24.2,22.4,20.4,18.0] ) * 6.8948 *1e9 # from psi to Pa

                i = 0
                while list_T[i]<T:
                        i = i+1

                a, b, c = np.polyfit( [list_T[i-1],list_T[i], list_T[i+1]] , [list_E[i-1],list_E[i], list_E[i+1]], 2)

                return a*T**2+b*T+c # Pa

        def E_carbon_steel_C_03_sound_speed(self, T):
                self.E = self.E_carbon_steel_C_03(T)
                self.sound_speed()


        def E_carbon_steel_C_03_sound_speed_shaft(self, T, l, D, k=1.):
                self.E = self.E_carbon_steel_C_03(T)
                self.sound_speed_shaft(l, D, k)


        def E_carbon_steel_C_mor_03(self,T):
                list_T = [-200,-129,-73,21,93,149,204,260,316,371,427,482,538,593,649] # *C
                list_E = np.array( [31.2,30.6,30.0,29.3,28.6,28.1,27.5,27.1,26.5,25.3,24.0,22.2,20.2,17.9,15.4] ) * 6.8948 *1e9 # from psi to Pa

                i = 0
                while list_T[i]<T:
                        i = i+1

                a, b, c = np.polyfit( [list_T[i-1],list_T[i], list_T[i+1]] , [list_E[i-1],list_E[i], list_E[i+1]], 2)

                return a*T**2+b*T+c # Pa

        def E_carbon_steel_C_mor_03_sound_speed(self, T):
                self.E = self.E_carbon_steel_C_mor_03(T)
                self.sound_speed()


        def E_carbon_steel_C_mor_03_sound_speed_shaft(self, T, l, D, k=1.):
                self.E = self.E_carbon_steel_C_mor_03(T)
                self.sound_speed_shaft(l, D, k)
                

                

# "A new topological structure for the Langevin-type ultrasonic transducer"; Xiaolong Lu, Junhui Hu, Hanmin Peng, Yuan Wang
steel = Stell(
                        name = "steel",
                        ρ = 7840.,      # density, kg/m3
                        E = 19.86e10,   # adopted, Young module, Pa
                        ν = 0.29,       # Poisson
                        λ = 83.5,       # Thermal conductivity, W / (m * K)
                        damping_coefficient = 1.09e-4,
                        cp = 640.       # cieplo walasciwe, J/kgK                        
                )


material_list.append(steel)





steel_40HM = Stell(
                        name = "steel_40HM",
                        ρ = 7850.,      # density, kg/m3
                        E = 19.86e10,   # adopted, Young module, Pa
                )

material_list.append(steel_40HM)



#----------------------------------------------------------- 1.2379

steel_1_2379 = Stell(
                        name = "steel_1_2379",
                        ρ = 7699,       # +-8 density, kg/m3 - measured
                        E = 209.2e9,      # Young module, Pa
                        ν = 0.27,       # adopted, Poisson
                )

material_list.append(steel_1_2379)

#----------------------------------------------------------- St52_3

steel_St52_3 = Stell(
                        name = "steel_St52_3",
                        ρ = 7685.6,      # density, kg/m3
                        E = 19.86e10,   # adopted, Young module, Pa
                        ν = 0.27,       # adopted, Poisson
                        λ = 83.5        # adopted, Thermal conductivity, W/(m*K)
                )

material_list.append(steel_St52_3)

#----------------------------------------------------------- S355

# http://www.meadinfo.org/2015/08/s355-steel-properties.html
# http://www.steelss.com/Carbon-steel/s355.html

steel_S355 = Stell(
                        name = "steel_S355",
                        ρ = 7850.,      # density, kg/m3
                        E = 205.3e9,     # wartosc zmierzona 190e9,      #0.95*190e9, # Young module, Pa 97.8e9
                        ν = 0.27,       # adopted, Poisson
                        λ = 54.         # Thermal conductivity, W/(m*K), http://www.s-k-h.com/media/de/Service/Werkstoffblaetter_englisch/Dickwand__Hohlprofile/Hollow_section_acc._to_10210.pdf
                )


material_list.append(steel_S355)

#----------------------------------------------------------- steel_tool

# Designing and modelling of the power ultrasonic transtucers
# Milan D. Radamanivić
# str str 1822


steel_tool = Stell(
                        name = "steel_tool",
                        ρ = 7850.,      # density, kg/m3
                        E = 218e9,      # Young module, Pa
                        ν = 0.29,       # adopted, Poisson
                        c = 5250,       # m.s
                )

material_list.append(steel_tool)

#----------------------------------------------------------- Stainless Steel 304

steel_304 = Stell(
                        name = "steel_316",
                        ρ = 8000.0,      # density, kg/m3
                        E = 193e9,      # Young module, Pa
                        ν = 0.29,       # adopted, Poisson
                )

material_list.append(steel_304)

#----------------------------------------------------------- Stainless Steel 316

steel_316 = Stell(
                        name = "steel_316",
                        ρ = 7970.0,      # density, kg/m3
                        E = 197.5e9,      # Young module, Pa
                        ν = 0.27,       # adopted, Poisson
                )

material_list.append(steel_316)

#----------------------------------------------------------- steel 50Hs
"""
ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
str 148
"""

steel_50Hs = Stell(
                        name = "steel_50Hs",
                        ρ = 7850.0,      # density, kg/m3
                        E = 199.0e9,      # Young module, Pa
                        ν = 0.29,       # adopted, Poisson
                )

material_list.append(steel_50Hs)

#*****************************************************************************
# carbon steel
#*****************************************************************************

def Thermal_conductivity__carbon_steel(T):
        """
        approx:
        http://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
        range od 25 do 225
        R = 0.993

        T -temperature, *C
        """
        return -0.035*T+55.0416666667 # W/(m*K)

#*****************************************************************************
# stainless steel
#*****************************************************************************

def Thermal_conductivity__stainless_steel(T):
        """
        approx:
        http://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
        range od 25 do 225
        T - temperature, *C
        """
        return 15.90625+0.00249999999999996*T+5.00000000000001e-05*T**2 # W/(m*K)



########################
# aluminum
########################

duralumin_AZW_7022 = Material(
                        name = "duralumin_AZW_7022",
                        ρ = 2820.,      # density, kg/m3
                        E = 72.6e9,     # Young module, Pa
                        ν = 0.34,       # Poisson
                )

material_list.append(duralumin_AZW_7022)


# http://www.dostal.com.pl/metale-kolorowe-aluminium.html
# http://www.zod.tu.koszalin.pl/dokumenty_pdf/materialy_uzupelniajace/materialoznawstwo/6_Stopy_CuAlMgLiTi.pdf

class ClassDuraluminPA6(Material):

        def speed_of_sound_approximation(self,l):
                """
                l - length, m
                """

                l = l*1e3

                a0=-111.933744254896
                a1=203.586581224638
                a2=-3.03243919134729
                a3=0.0200909714454779
                a4=-4.94663300246339e-05

                return a0+a1*l+a2*l**2+a3*l**3+a4*l**4

        def resonant_frequency_shaft_fi50(self, l):
                """
                l - length, m
                """
                l *= 1e3
                return 91.7138195346413-1.10384776001922*l+0.00431726494543734*l**2

        

        


duralumin_PA6 = ClassDuraluminPA6(
                        name = "duralumin_PA6",
                        ρ = 2811.,      # density, kg/m3 - 2811 +- 3 measured 2018.02.27, 2790 on page 
                        E = 76.2e9,      #75.9e9,            #74.35e9,    #1.3*72.5e9, # Young module, Pa
                        ν = 0.33,       # Poisson
                        cp = 873.,      # Heat capacity, J/kgK
                        α = 22.9e-9,    # thermal expansion, 1/K
                        λ = 134.,       # Thermal conductivity, W/(m*K)
                        c = 0.
                )


material_list.append(duralumin_PA6)



aluminum = Material(
                        name = "aluminum",
                        ρ = 2720.,      # density, kg/m3
                        E = 69.e9,      # Young module, Pa
                )

material_list.append(aluminum)


#----------------------------------------------------------- duraluminium

# Designing and modelling of the power ultrasonic transtucers
# Milan D. Radamanivić
# str str 1822


duraluminium = Material(
                        name = "duraluminium",
                        ρ = 2790.,      # density, kg/m3
                        E = 740e9,      # Young module, Pa
                        ν = 0.34,       # adopted, Poisson
                        c = 5250,       # m.s
                )

material_list.append(duraluminium)


#--------------------------------------------------
# http://www.ceramit.pl/wlasnosci.html
# Thermal conductivity, W / (m * K)

Al2O3_λ = 25.   # Thermal conductivity, W / (m * K)
ZrO2_λ = 2.     # Thermal conductivity, W / (m * K)
Si3N4_λ = 18.   # Thermal conductivity, W / (m * K)


#----------------------------------------------------------

########################
# titanium
########################

titanium = Material(
                        name = "titanium",
                        ρ = 4506.,      # density, kg/m3 - 2811 +- 3 measured 2018.02.27, 2790 on page 
                        E = 116e9,            #74.35e9,    #1.3*72.5e9, # Young module, Pa
                        ν = 0.32,       # Poisson
                        λ = 21.9,       # Thermal conductivity, W/(m*K)
                )

material_list.append(titanium)




#---------------------------------------------------------- Tytan titanium grade 2

titanium_grade_2 = Material(
                        name = "titanium",
                        ρ = 4510.,      # density, kg/m3 - 2811 +- 3 measured 2018.02.27, 2790 on page 
                        E = 105e9,            #74.35e9,    #1.3*72.5e9, # Young module, Pa
                        ν = 0.37,       # Poisson
                        λ = 16.4,       # Thermal conductivity, W/(m*K)
                )

material_list.append(titanium)

#----------------------------------------------------------


########################
# glass
########################

#---------------------------------------------------------- Young's modulus


# glass
glass = Material(
                        name = "glass",
                        ρ = 2600.,      # adopted density, kg/m3
                        E = 72.e9       # Young module, Pa
                )

material_list.append(glass)


def Thermal_conductivity__fiberglass_econo_1260(T):
        """
        T - temperature, *C
        """
        return (0.2025*T-10.)*0.001 # W/(m*K)


fiberglass_composite = Material(
                        name = "fiberglass_composite",
                        ρ = 2600.,      # adopted density, kg/m3
                        E = 150.e9,     # Young module, Pa
                        λ = Thermal_conductivity__fiberglass_econo_1260(20)
                )

material_list.append(fiberglass_composite)


# https://www.pgo-online.com/intl/macor-machinable-glass-ceramic.html
macor = Material(
                        name = "macor",
                        ρ = 2520.,      # density, kg/m3
                        E = 66.9e9,     # Young module, Pa
                        ν = 0.29,       # Poisson
                        λ = 1.46,       # Thermal conductivity, W/(m*K)
                )

material_list.append(macor)


# https://www.pgo-online.com/intl/BK7.html
BK7Schott = Material(
                        name = "BK7Schott",
                        ρ = 2510.,      # density, kg/m3
                        E = 82.e9,      # Young module, Pa
                        ν = 0.206,      # Poisson
                        λ = 1.114,      # Thermal conductivity, W/(m*K)
                )

material_list.append(BK7Schott)


########################
# piezoceramic
########################





# PIC18http://fisica.cab.cnea.gov.ar/bt/images/d/d3/PICat.pdf1 - ceramika twarda
# https://books.google.pl/books?id=Q4z4DAAAQBAJ&pg=PA127&lpg=PA127&dq=PIC181+material+properties&source=bl&ots=EqhaY2Y-c9&sig=8HAvebwnhBO_kefoG8CF-EqVTQ4&hl=pl&sa=X&ved=0ahUKEwjq_OWP7_jSAhWBHywKHTVgAxEQ6AEIQDAF#v=onepage&q=PIC181%20material%20properties&f=false
# Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
# http://piezomat.org/materials/273

# https://www.google.pl/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0ahUKEwiR0uHfu7XTAhXhK5oKHf3CCYoQFgg5MAE&url=https%3A%2F%2Fwww.piceramic.com%2Fen%2Fproducts%2Fpiezoelectric-materials%2F%3Ftype%3D5600%26downloadUid%3D1111%26downloadFileUid%3D2281%26cHash%3D62c0ee110a0c7980e1c29d99e080f216&usg=AFQjCNE4J1np_3aHNkbtfHUAkVy5JFYyDA

PIC181_c = 4629.63
PIC181_cp = 350. # Pojemnosc cieplna, J/(kg*K)
PIC181_λ = 1.1 # Thermal conductivity, W / (m * K)
PIC181_ν = 0.34 # liczba Poissona
PIC181_α3 = -5e-6 # rozszezalnoc cieplna w kierunku polaryzacji, aprox. -4 to 6 e-6, 1/K
PIC181_α1 = 6e-6 # rozszezalnoc cieplna w kierunku prostopadlym polaryzacji, aprox. 4 to 8 e-6, 1/K

PIC181_ρ = 7856.7 #7850. # kg/m3 od producenta 7800., wartosc z Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł 
PIC181_Tc = 330. # *C
PIC181_eT33pe0 = 1200.
PIC181_eT33 = PIC181_eT33pe0*e0
PIC181_eT11pe0 = 1500. # http://piezomat.org/materials/273
PIC181_eT11 = PIC181_eT11pe0*e0 # http://piezomat.org/materials/273
PIC181_e_e = 1500.
PIC181_tan_delta = 3.e-3

PIC181_kp = 0.56
PIC181_kt = 0.46
PIC181_k31 = 0.32
PIC181_k33 = 0.66
PIC181_k15 = 0.63

PIC181_d31 = -1.08e-10 # C/N, m*V-1 od producenta -120.*10**-12, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł 
PIC181_d33 = 2.53e-10 # C/N, m*V-1 od producenta 265.*10**-12, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł 
PIC181_d15 = 3.89e-10 # C/N, m*V-1 od producenta 475.*10**-12, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_e31 = -4.5 # N*V-1*m^-1, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_e33 = 14.7 # N*V-1*m^-1, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_e15 = 11. # N*V-1*m^-1, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_g31 = -11.2*10**-3 # Vm/N
PIC181_g33 = 25.*10**-3 # Vm/N

PIC181_Np = 2270. # Hzm
PIC181_N1 = 1640. # Hzm
PIC181_N3 = 2010. # Hzm
PIC181_Nt = 2110. # Hzm

PIC181_SE11 = 1.175e-11 # m2/N od producenta 11.8e-12, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł 
PIC181_SE12 = -4.07e-12 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SE13 = -4.996e-12 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł 
PIC181_SE33 = 1.411e-11 # m2/N od producenta 14.2e-12, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SE44 = 3.533e-11 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SE66 = 2.*(PIC181_SE11-PIC181_SE12) # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_YE11 = 1./PIC181_SE11 # Young module, Pa
PIC181_YE12 = 1./PIC181_SE12 # Young module, Pa
PIC181_YE13 = 1./PIC181_SE13 # Young module, Pa
PIC181_YE33 = 1./PIC181_SE33 # Young module, Pa
PIC181_YE44 = 1./PIC181_SE44 # Young module, Pa
PIC181_YE66 = 1./PIC181_SE66 # Young module, Pa

PIC181_SD11 = 1.058e-11 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SD12 = -5.235e-12 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SD13 = -2.268e-12 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SD33 = 1.137e-11 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SD44 = 2.134e-11 # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_SD66 = 2.*(PIC181_SD11-PIC181_SD12) # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_YD11 = 1./PIC181_SD11 # Young module, Pa
PIC181_YD12 = 1./PIC181_SD12 # Young module, Pa
PIC181_YD13 = 1./PIC181_SD13 # Young module, Pa
PIC181_YD33 = 1./PIC181_SD33 # Young module, Pa
PIC181_YD44 = 1./PIC181_SD44 # Young module, Pa
PIC181_YD66 = 1./PIC181_SD66 # Young module, Pa

PIC181_cE11 = 152.3e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cE12 = 89.1e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cE13 = 85.5e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cE33 = 134.e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cE44 = 28.3e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cE66 =  0.5*(PIC181_cE11-PIC181_cE12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł



PIC181_cD11 = 155.e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cD12 = 91.8e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cD13 = 70.6e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cD33 = 16.6e10#1.02*166.4e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cD44 = 46.9e9 # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_cD66 = 0.5*(PIC181_cD11-PIC181_cD12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł


PIC181_Qm = 2000.
PIC181_TK_S33 = 3.e-3 # 1/K

PIC181_eT11 = 12241.*e0 # F*m^-1 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_eT33 = 1135.*e0 # F*m^-1 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_eS11 = 740*e0 # F*m^-1, od producenta 1200*e0, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł
PIC181_eS33 = 634*e0 # F*m^-1, Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC181_h33 = 26.8e8 # PIC181_e33/PIC181_eS33 #2.5e8 * PIC181_kt * (PIC181_cD33 * PIC181_eS33)**0.5  # 26.8e8 # jako PZT8- 





class Piezoceramic(Material):

        def __init__(self, name="",cp=PIC181_cp, λ=PIC181_λ, ν=PIC181_ν, α3=PIC181_α3, α1=PIC181_α1, ρ=PIC181_ρ, Tc=PIC181_Tc, eT33pe0=PIC181_eT33pe0, eT33=PIC181_eT33, eT11pe0=PIC181_eT11pe0, e_e=PIC181_e_e, tan_delta=PIC181_tan_delta, kp=PIC181_kp, kt=PIC181_kt, k31=PIC181_k31, k33=PIC181_k33, k15=PIC181_k15,d31=PIC181_d31, d33=PIC181_d33, d15=PIC181_d15,e31=PIC181_e31, e33=PIC181_e33, e15=PIC181_e15,g31=PIC181_g31, g33=PIC181_g33,Np=PIC181_Np, N1=PIC181_N1, N3=PIC181_N3, Nt=PIC181_Nt,SE11=PIC181_SE11, SE12=PIC181_SE12, SE13=PIC181_SE13, SE33=PIC181_SE33, SE44=PIC181_SE44, SE66=PIC181_SE66,SD11=PIC181_SD11, SD12=PIC181_SD12, SD13=PIC181_SD13, SD33=PIC181_SD33, SD44=PIC181_SD44, SD66=PIC181_SD66,cE11=PIC181_cE11, cE12=PIC181_cE12, cE13=PIC181_cE13, cE33=PIC181_cE33, cE44=PIC181_cE44, cE66=PIC181_cE66,cD11=PIC181_cD11, cD12=PIC181_cD12, cD13=PIC181_cD13, cD33=PIC181_cD33, cD44=PIC181_cD44, cD66=PIC181_cD66,Qm=PIC181_Qm, TK_S33=PIC181_TK_S33,eT11=PIC181_eT11, eS11=PIC181_eS11, eS33=PIC181_eS33, vz=0., C0=0, N=0, vd=0, gamma_p_w=0, Mw=0, damping_coefficient=0, cE55=0, YE11=PIC181_YE11, YE33=PIC181_YE33, c=PIC181_c, Z0=0, CT0=0, h31=0, h33=PIC181_h33):
                self.name = name
                self.cp=cp
                self.λ=λ
                self.ν=ν
                self.α3=α3
                self.α1=α1
                self.ρ=ρ
                self.Tc=Tc
                self.eT33pe0=eT33pe0
                self.eT33=eT33
                self.eT11pe0=eT11pe0
                self.e_e=e_e
                self.tan_delta=tan_delta
                self.kp=kp
                self.kt=kt
                self.k31=k31
                self.k33=k33
                self.k15=k15
                self.d31=d31
                self.d33=d33
                self.d15=d15
                self.e31=e31
                self.e33=e33
                self.e15=e15
                self.g31=g31
                self.g33=g33
                self.Np=Np
                self.N1=N1
                self.N3=N3
                self.Nt=Nt
                
                self.SE11=SE11
                self.SE12=SE12
                self.SE13=SE13
                self.SE33=SE33
                self.SE44=SE44
                self.SE66=SE66

                self.YE11 = 1./self.SE11
                self.YE12 = 1./self.SE12
                self.YE13 = 1./self.SE13
                self.YE33 = 1./self.SE33
                self.YE44 = 1./self.SE44
                self.YE66 = 1./self.SE66
                
                self.SD11=SD11
                self.SD12=SD12
                self.SD13=SD13
                self.SD33=SD33
                self.SD44=SD44
                self.SD66=SD66

                self.YD11 = 1./self.SD11
                self.YD12 = 1./self.SD12
                self.YD13 = 1./self.SD13
                self.YD33 = 1./self.SD33
                self.YD44 = 1./self.SD44
                self.YD66 = 1./self.SD66
                
                self.cE11=cE11
                self.cE12=cE12
                self.cE13=cE13
                self.cE33=cE33
                self.cE44=cE44
                self.cE66=cE66
                self.cD11=cD11
                self.cD12=cD12
                self.cD13=cD13
                self.cD33=cD33
                self.cD44=cD44
                self.cD66=cD66
                self.Qm=Qm
                self.TK_S33=TK_S33
                self.eT11=eT11
                self.eS11=eS11
                self.eS33=eS33
                self.vz=vz
                self.C0=C0
                self.N=N
                self.vd=vd
                self.gamma_p_w=gamma_p_w
                self.Mw=Mw
                self.damping_coefficient=damping_coefficient
                self.cE55=cE55
                self.YE11=YE11
                self.YE33=YE33
                self.Z0=Z0
                self.CT0=CT0

                self.E = YE33

                if h31 != 0:
                        self.h31 = h31
                else:
                        self.h31 = self.g31/self.SD11
                        
                self.h33 = h33

                if c==0:
                        self.c = ( 1/(self.SE33*self.ρ) )**0.5
                else:
                        self.c = c

                self.impedance()

        def piezoelectric_effect(self, T, E=0, d=0, eT=0):
                """
                D = d*T + eT*E

                T - stress, Pa
                d - piezoelectric module, C/N or V/m
                E - intensity of the electric field, V/m
                eT - electrical permeability for permanent stress, F/m

                return D
                """
                if d == 0:
                        d = self.d33

                if eT == 0:
                        eT = self.eT33

                return d*T + E*eT

        def piezoelectric_effect_reversed(self, E, T=0, d=0, sE=0):
                """
                S = d*E + sE*T

                T - stress, Pa
                d - piezoelectric module, C/N or V/m
                E - intensity of the electric field, V/m
                sE - mechanical susceptibility with constant electric field strength, Pa^-1

                return D
                """
                if d == 0:
                        d = self.d33

                if sE == 0:
                        sE = self.SE33

                return d*E + sE*T

        def phase_velocity(self):
                """
                A. Milewski, P. Kogut, W. Kardyś, P. Kluk
                „Eksperymentalna walidacja elektromechanicznych modeli ultradźwiękowych przetworników piezoceramicznych”
                Elektronika Nr 8/2012.
                """
                
                return ( (self.cE33 + (self.e33**2)/(self.eS33))/self.ρ )**0.5

        def external_force(self, f, A, v1, v2, l, I):
                """
                A. Milewski, P. Kogut, W. Kardyś, P. Kluk
                „Eksperymentalna walidacja elektromechanicznych modeli ultradźwiękowych przetworników piezoceramicznych”
                Elektronika Nr 8/2012.
                """

                h = l/2.

                ω = 2.*np.pi*f
                vz = self.phase_velocity()
                k3 = ω/vz                
                
                F1 = -1j*self.ρ*vz*A*( (1./np.tan(2.*h*k3))*v1 + (1./np.sin(2.*h*k3))*v2 ) + self.e33/(1j*ω*self.eS33)*I
                
                F2 = -1j*self.ρ*vz*A*( (1./np.sin(2.*h*k3))*v1 + (1./np.tan(2.*h*k3))*v2 ) + self.e33/(1j*ω*self.eS33)*I
        
                U = 1./(1j*ω) * ( (self.e33/self.eS33)*v1 + (self.e33/self.eS33)*v2 + 2.*h/(self.eS33*A)*I)

                return F1, F2, U


        def electromechanical_coefficient_axial_vibrations(self):
                return (1./(PIC181.eS33*PIC181.cE33/PIC181.e33**2 + 1))**0.5

        def calculation_of_ceramic_parameters(self, fs, fp, D=50e-3, d=19.6e-3, l=5e-3, C0=3.6e-9, m=65.253e-3):

                """
                A. Milewski, P. Kogut, W. Kardyś, P. Kluk
                „Eksperymentalna walidacja elektromechanicznych modeli ultradźwiękowych przetworników piezoceramicznych”
                Elektronika Nr 8/2012.

                fs - resonance, Hz
                fp - antyresonance, Hz
                D - diameter, m
                d - hole diameter, m
                l - length, m
                C0 - capacity, F
                m - mass, kg                
                """
                
                h = l/2.
                
                A = np.pi*(D/2)**2 - np.pi*(d/2)**2

                self.ρ = m/(A*l)

                kt2 = (np.pi/2.) * (fs/fp) * np.tan(np.pi/2.*(fp-fs)/fp)

                self.kt = kt2**0.5
                
                self.cE33 = (1.-kt2)*(4*h*fp)**2 * self.ρ
                self.cD33 = (4*h*fp)**2 * self.ρ

                self.eT33 = ( (2*h*C0)/A )

                self.eS33 = self.eT33*(1.-kt2)*(1.-0.55**2)

                self.e33 = ( kt2*self.cE33*self.eS33/(1.-kt2) )**0.5

                self.h33 = self.e33/self.eS33

                self.vz = ( (self.cE33 + (self.e33**2)/(self.eS33))/self.ρ )**0.5

                return self.ρ, self.eT33, self.eS33, self.cE33, self.cD33, self.e33, self.h33, self.vz, self.kt


        def __electrical_impedance_axial_vibrations_Kogut_round_ceramics(self, f, D=50e-3, d=19e-3, l=5e-3, Z1=0., Z2=0.):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                """
                        
                l = 0.5*l # strona 29 rysunek 5.3, grubosc to 2h - tu stosowane l
                ω = 2.*np.pi*f
        
                self.vz = self.phase_velocity()
                kz = ω/self.vz

                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2

                CS0 = self.eS33*A/(2.*l)

                self.kt2 = (self.e33**2)/(self.eS33*self.cD33)

                Zel = (1./(1j*ω*CS0)) * ( 1. - self.kt2/(l*kz) * ( (1./2j) * (Z1+Z2) + self.ρ*self.vz*np.tan(l*kz) )/( (Z1*Z2)/(self.ρ*self.vz) - 1j*(Z1+Z2)/(np.tan(2.*l*kz)) + self.ρ*self.vz ) )

                return Zel


        def __electrical_impedance_axial_vibrations_Kogut_round_ceramics_z1_and_z1_array(self, i, D=50e-3, d=19e-3, l=5e-3):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                """
                        
                l = 0.5*l # strona 29 rysunek 5.3, grubosc to 2h - tu stosowane l
                ω = 2.*np.pi*self.f[i]
        
                self.vz = self.phase_velocity()
                kz = ω/self.vz

                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2

                CS0 = self.eS33*A/(2.*l)

                self.kt2 = (self.e33**2)/(self.eS33*self.cD33)

                Zel = (1./(1j*ω*CS0)) * ( 1. - self.kt2/(l*kz) * ( (1./2j) * (self.Z1[i]+self.Z2[i]) + self.ρ*self.vz*np.tan(l*kz) )/( (self.Z1[i]*self.Z2[i])/(self.ρ*self.vz) - 1j*(self.Z1[i]+self.Z2[i])/(np.tan(2.*l*kz)) + self.ρ*self.vz ) )

                return Zel
                

        def electrical_impedance_axial_vibrations_Kogut_round_ceramics(self, f=np.arange(1e3, 600e3, 100), D=50e-3, d=19.6e-3, l=5e-3, Z1=429., Z2=429.):
                """
                A. Milewski, P. Kogut, W. Kardyś, P. Kluk
                „Eksperymentalna walidacja elektromechanicznych modeli ultradźwiękowych przetworników piezoceramicznych”
                Elektronika Nr 8/2012.

                f - array frequence, Hz
                D - diameter, m
                d - hole diameter, m
                l - length, m
                Z1 - impedance, Rayl, default 429 for air
                Z2 - impedance, Rayl, default 429 for air
                """
                
                impedancja = np.vectorize( self.__electrical_impedance_axial_vibrations_Kogut_round_ceramics, excluded=["D", "d", "l", "Z1", "Z2"])

                Z = impedancja(f, D, d, l, Z1, Z2)

                #return f, np.real(Zul_dB(Z)), np.real(Z), np.degrees(np.angle(Z)), Z
                return f, Z


        def electrical_impedance_axial_vibrations_Kogut_round_ceramics_z1_and_z1_array(self, f=np.arange(1e3, 600e3, 100), D=50e-3, d=19.6e-3, l=5e-3, Z1=429.*np.ones(len(np.arange(1e3, 600e3, 100))), Z2=429.*np.ones(len(np.arange(1e3, 600e3, 100)))):
                """
                A. Milewski, P. Kogut, W. Kardyś, P. Kluk
                „Eksperymentalna walidacja elektromechanicznych modeli ultradźwiękowych przetworników piezoceramicznych”
                Elektronika Nr 8/2012.

                f - array frequence, Hz
                D - diameter, m
                d - hole diameter, m
                l - length, m
                Z1 - impedance, Rayl, default 429 for air
                Z2 - impedance, Rayl, default 429 for air
                """
                
                impedancja = np.vectorize( self.__electrical_impedance_axial_vibrations_Kogut_round_ceramics_z1_and_z1_array, excluded=["D", "d", "l"])

                self.f = f
                self.Z1 = Z1
                self.Z2 = Z2
                i = range(len(f))
                
                Z = impedancja(i, D, d, l)

                #return f, np.real(Zul_dB(Z)), np.real(Z), np.degrees(np.angle(Z)), Z
                return f, Z

        def electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot(self, f=np.arange(1e3, 600e3, 100), D=50e-3, d=19.6e-3, l=5e-3, Z1=429., Z2=429.):

                f, Z = self.electrical_impedance_axial_vibrations_Kogut_round_ceramics(f, D, d, l, Z1, Z2)

                print(f[np.absolute(Z).argmin()]*1e-3)

                pyplot.semilogy(f*1e-3, np.absolute(Z))
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("|Z|, $\Omega$")
                pyplot.show()


        def electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot_z1_and_z1_array(self, f=np.arange(1e3, 600e3, 100), D=50e-3, d=19.6e-3, l=5e-3, Z1=429.*np.ones(len(np.arange(1e3, 600e3, 100))), Z2=429.*np.ones(len(np.arange(1e3, 600e3, 100)))):

                f, Z = self.electrical_impedance_axial_vibrations_Kogut_round_ceramics_z1_and_z1_array(f, D, d, l, Z1, Z2)

                pyplot.semilogy(f*1e-3, np.absolute(Z))
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("|Z|, $\Omega$")
                pyplot.show()


        def __electrical_impedance_radial_vibrations_Kogut_round_ceramics(self, f, D=50e-3, d=19e-3, l=5e-3):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                """
        
                h = 0.5*l
                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2

                #------------------------------------------------------ 
                # effective material constants
                # "ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających"
                # page number 32, equations (5.27)
                cP11 = self.cE11-(self.cE13**2)/self.cE33 
                cP12 = self.cE12-(self.cE13**2)/self.cE33

                eP31 = self.e31 - self.e33*(self.cE13**2)/self.cE33
                #eP31 *= 1e-10
                εP33 = self.eS33 + (self.e33**2)/self.cE33
                
                #------------------------------------------------------

                w = 2.*np.pi*f

                vr = (cP11/self.ρ)**0.5
                self.vr = vr
                kr = w/vr

                oP = cP12/cP11 # -SE12/SE11
                kp2 = (eP31**2)/(εP33*cP11)
                Cp0 = εP33*A/(2.*h)

                
                def jr(r, kr, oP):
                        return kr*j0(kr*r) + (oP-1.)*j1(kr*r)/r # scipy.special.j0, scipy.special.j1
                        #return kr*jn(0, kr*r) + (oP-1.)*jn(1, kr*r)/r # scipy.special.j0, scipy.special.j1

                def yr(r, kr, oP):
                        return kr*y0(kr*r) + (oP-1.)*y1(kr*r)/r # scipy.special.y0, scipy.special.y1
                        #return kr*yn(0, kr*r) + (oP-1.)*yn(1, kr*r)/r # scipy.special.y0, scipy.special.y1

                a = D/2
                b = d/2

                ja = jr(a, kr, oP)
                jb = jr(b, kr, oP)

                ya = yr(a, kr, oP)
                yb = yr(b, kr, oP)

                A1 = (ya-yb)/(jb*ya-ja*yb)
                A2 = -(ja-jb)/(jb*ya-ja*yb)

                Delta_J1 = a*j1(kr*a)-b*j1(kr*b)
                Delta_Y1 = a*y1(kr*a)-b*y1(kr*b)

                Zel = 1./( 1j*w*Cp0 * (2.*np.pi*kp2*(A1*Delta_J1+A2*Delta_Y1)/A + 1.) )

                return Zel



        def electrical_impedance_radial_vibrations_Kogut_round_ceramics(self, f=np.arange(15e3,0.6e6,100), D=50e-3, d=19e-3, l=5e-3):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                """
        
                impedancja = np.vectorize( self.__electrical_impedance_radial_vibrations_Kogut_round_ceramics, excluded=["D", "d", "l"])

                Z = impedancja(f, D, d, l)

                #return f, np.real(Zul_dB(Z)), np.real(Z), np.degrees(np.angle(Z)), Z
                return f, Z



        #------------------------------------------- Mason model


        def __model_Mason(self, f=1e6, D=50e-3, d=19.6e-3, l=5e-3, Z1_load=liquids.air.Z, Z2_load=liquids.air.Z, p=1):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 39
                """

                cD33=self.cD33
                q=self.ρ
                eS33=self.eS33
                h33=self.h33
                                                
                
                
                vz = (cD33/q)**0.5
                
                λ = vz/f
                w = 2.*np.pi*f
                kz = w/vz

                A = np.pi*(D/2.)**2 - np.pi*(d/2.)**2
                
                
                if p == 1:
                        l = l/2.
                        
                        Zc = q*vz*A
                        
                        CS0 = eS33*A/(2*l)
        
                        n = h33*CS0
        
                        Zx = -1./(1j*w*CS0)

                

                        Z1a = Z1_load * A #* np.tan(w/Z1_load.c * 10)
                        Z1 = 1j*Zc*np.tan(p*kz*l) + Z1a
                
                        Z2a = Z2_load * A #* np.tan(w/Z1_load.c * 10)
                        Z2 = -1j*Zc/np.sin(2*p*kz*l) + Z2a

                else:
                        h = l/2
                        L = p*l
                        CS0 = eS33*A/L
                        n = h33*CS0/p
                        Zc = q*vz
                        y = 2.*np.arcsin( h*(-0.5*np.tan(L*kz/(2*p))*np.sin(L*kz/p)) )

                        Z1a = Z1_load * A
                        Z1 = Zc*np.tanh(p*h*y) + Z1a

                        Z2a = Z2_load * A #* np.tan(w/Z1_load.c * 10)
                        Z2 = Zc/(np.sinh(2*p*y)) + Z2a
                
                        
                        Zx = -1./(1j*w*CS0)
                        
                                
                return Z1, Z2, CS0, Zx, n
                

        def model_Mason(self, f=np.arange(15e3,0.6e6,100), D=50e-3, d=19e-3, l=5e-3, p=1, Z1_load=liquids.air.Z, Z2_load=liquids.air.Z):
                """
                ROZPRAWA DOKTORSKA mgr inż. Paweł Kogut Metody Modelowania i Projektowania Ultradźwiękowych Układów Drgających
                str 39

                f - array frequence, Hz
                D - diameter, m
                d - hole diameter, m
                l - length, m
                p - amount of ceramics
                """
                
                impedancja = np.vectorize( self.__model_Mason, excluded=['D', 'd', 'l', 'Z1_load', 'Z2_load', 'p'])

                Z1, Z2, CS0, Zx, n = impedancja(f, D, d, l, Z1_load, Z2_load, 1)

                Z1 = Z1/p
                Z2 = Z2/p
                Zx = Zx/p

                Z1 = 1./(1/Z1+1/Z1)

                Z = Z2 + Z1
                #n = h33*CS0
                Z = Z/n**2 # (Zp/Zs) = Np/Ns
                Z = Zx + Z
                Zc = 1./(1j*2.*np.pi*f*CS0)
                Z = 1./(1./Z + p/Zc) # p/Zc

                return Z, Z1, Z2, CS0, Zx, n, Zc


        def model_Mason_plot(self, f=np.arange(15e3,0.6e6,100), D=50e-3, d=19e-3, l=5e-3, p=1, Z1_load=liquids.air.Z, Z2_load=liquids.air.Z):

                Z, Z1, Z2, CS0, Zx, n, Zc = self.model_Mason(f=f, D=D, d=d, l=l, p=p, Z1_load=Z1_load, Z2_load=Z2_load)

                #c1 = np.loadtxt("D:/Tomasz/dokumentacje/ceramiki/PIC181/pic181 - nr4 - D50 d19.6 t5/ceramika 1/skan_Hz.dat", unpack=True)

                #pyplot.semilogy(c1[0]*1e-3, c1[1], color="blue", label="ceramika 1")

                print(f[np.absolute(Z).argmin()]*1e-3)

                pyplot.semilogy(f*1e-3, np.absolute(Z))
                pyplot.xlabel("f, kHz")
                pyplot.ylabel("Z, $\Omega$")
                pyplot.show()


        #------------------------------------------------------ three-dimensional Model of Piezoceramic Ring, Milan D. Radanovic "Designing and modeling of the power ultrasonic transducers" page 46

        def __three_dimensional_Model_of_Piezoceramic_Ring(self, f, a, b, h, v1, v2, v3, v4):
                """
                three-dimensional Model of Piezoceramic Ring, Milan D. Radanovic "Designing and modeling of the power ultrasonic transducers" page 46
                """

                w = 2.*np.pi*f

                vr = (self.cD11/self.ρ)**0.5
                vz = (self.cD33/self.ρ)**0.5

                kz = w/vz
                kr = w/vr
                
                P = np.pi*(a**2 - b**2)
                C0 = self.eS33*P/(2*h)
                
                z35 = self.h33/(1j*w)
                z55 = 1./(1j*w*C0)

                z33 = self.cD33*kz*P/(1j*w*np.tan(2.*kz*h))
                z34 = self.cD33*kz*P/(1j*w*np.sin(2.*kz*h))

                z23 = 2.*np.pi*a*self.cD13/(1j*w)
                z25 = 4.*np.pi*a*h*self.h31/(1j*w*P)

                z13 = 2.*np.pi*b*self.cD13/(1j*w)
                z15 = 4.*np.pi*b*h*self.h31/(1j*w*P)

                A0 = j1(kr*b)*y1(kr*a)-j1(kr*a)*y1(kr*b)
                A1 = y1(kr*a)/A0
                A2 = y1(kr*b)/A0

                B0 = j1(kr*a)*y1(kr*b)-j1(kr*b)*y1(kr*a)
                B1 = j1(kr*a)/B0
                B2 = j1(kr*b)/B0

                z21 = ( -4.*np.pi*kr*b*h*self.cD11/(1j*w) ) * ( A1*j0(kr*a) + B1*y0(kr*a) )
                
                z12 = ( -4.*np.pi*kr*b*h*self.cD11/(1j*w) ) * ( A2*j0(kr*b) + B2*y0(kr*b) )

                z22 = ( 4.*np.pi*h/(1j*w) ) * ( self.cD12-self.cD11*(1+kr*a*( A2*j0(kr*a) + B2*y0(kr*a) )) )

                z22 = ( -4.*np.pi*h/(1j*w) ) * ( self.cD12-self.cD11*(1-kr*b*( A1*j0(kr*b) + B1*y0(kr*b) )) )

                Z = np.matrix('z11 z12 z13 z13 z15; z21 z22 z23 z23 z25; z13 z23 z33 z34 z35; z13 z23 z34 z33 z35; z15 z25 z35 z35 z55')

                v = np.matrix('v1; v2; v3; v4; 1')

                F = Z*v

                return F[-1]

                
                

        def three_dimensional_Model_of_Piezoceramic_Ring(self, f=np.arange(15e3,600e3,10), D=50e-3, d=19.6e-3, l=5e-3):

                #impedancja = np.vectorize( self.__three_dimensional_Model_of_Piezoceramic_Ring, excluded=['a', 'b', 'h', 'v1', 'v2', 'v3', 'v4'])

                #Z = impedancja(f=f, a=D/2., b=d/2., h=l/2., v1=liquids.air.c, v2=liquids.air.c, v3=liquids.air.c, v4=liquids.air.c)

                print(self.__three_dimensional_Model_of_Piezoceramic_Ring(15e3, a=16.5e-3,b=7.5e-3,h=2.5e-3,v1=340.,v2=340., v3=340., v4=340.))
                #print(self.__three_dimensional_Model_of_Piezoceramic_Ring(430e3, a=16.5e-3,b=7.5e-3,h=2.5e-3,v1=340.,v2=340., v3=340., v4=340.))

                #pyplot.plot(f, Z)
                #pyplot.show()

        
                
                
PIC181 = Piezoceramic()



# SonoxP8 - Kogut
SonoxP8 = Piezoceramic(name = "",ρ=7513., eT33=9.1e-9, eS33=5.25e-9,
           cE33=12.8e10, cD33=15.49e10,
           e33=11.87, h33=22.6e8,
           vz=4540., kt=0.416
)




# Motorola 3203HD -
Motorola_3203HD = Piezoceramic(
                        name = "",
                        ρ = 7800., # kg/m3
                        cD33 = 1.77*(1.+0.023j)*10**11, # N/m2
                        eS33 = 1.06*(1.-0.053j)*10**-8, # F/m
                        h33 = 2.19*(1.+0.029j)*10**9, # V/m
                        kt = 0.536*(1.-0.005j),
                        C0 = 1.87*(1.-0.035j)*1e-9,# F
                        N = 4.11*(1.-0.024j), # C/m
                        vd = 4674.*(1.+0.012j), # m/s
                        gamma_p_w = 2.1*(1.-0.012j)*10**-4, # s/m
                        Mw = 3.33*(1+0.017j)*10**5 # Vs/mkg
)




# PZT8
# Morgan Electro Ceramics Web Site
# www.morgan-electroceramics.com
# PROPERTIES OF MORGAN
# ELECTRO CERAMIC CERAMICS
# D. Berlincourt and H. H. A. Krueger revised by
# C. Near

PZT8_ρ = 7600. # kg/m3
PZT8_Tc = 300. # C
PZT8_eT33pe0 = 1000.
PZT8_eT33 = PZT8_eT33pe0*e0
PZT8_eT11pe0 = 1290.
PZT8_eT11 = PZT8_eT11pe0*e0
PZT8_eS33pe0 = 582.
PZT8_eS33 = PZT8_eS33pe0*e0
PZT8_eS11pe0 = 900.
PZT8_eS11 = PZT8_eS11pe0*e0
PZT8_tan_delta = 4.e-3

PZT8_kp = -0.51
PZT8_kt = 0.48
PZT8_k31 = -0.30
PZT8_k33 = 0.64
PZT8_k15 = 0.55
PZT8_d31 = -37.e-10 # C/N
PZT8_d33 = 225.e-10 # C/N
PZT8_d15 = 330.e-10 # C/N
PZT8_dh = 31.e-10 # C/N
PZT8_g31 = -10.9e-3 # Vm/N
PZT8_g33 = 25.4e-3 # Vm/N
PZT8_g15 = 28.9e-3 # Vm/N
PZT8_gh = 3.6e-3 # Vm/N

PZT8_h31 = -7.8e8 # -7.7e-10 # V/m
PZT8_h33 = 26.9e8 # 26.4e-10 # V/m
PZT8_his = 12.9e-10 # V/m

PZT8_e31 = -4. # -4.1 # C/m^2
PZT8_e33 = 13.8 # 14.0 # C/m^2
PZT8_e15 = 10.3 # C/m^2

PZT8_SE11 = 11.5e-12 # m2/N
PZT8_SE33 = 13.5e-12 # m2/N
PZT8_SE44 = 31.9e-12 # m2/N
PZT8_SE66 = 30.4e-12 # m2/N
PZT8_SE12 = -3.38e-12 # -3.7e-12 # m2/N
PZT8_SE13 = -4.69e-12 #-4.8e-12 # m2/N
#PZT8_SE66 = 2.*(PZT8_SE11-PZT8_SE12) # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PZT8_SD11 = 10.1e-12 # m2/N
PZT8_SD33 = 8.5e-12 # m2/N
PZT8_SD44 = 22.6e-12 # m2/N
PZT8_SD12 = -4.5e-12 # m2/N
PZT8_SD13 = -2.5e-12 # m2/N
PZT8_SD66 = 2.*(PZT8_SD11-PZT8_SD12) # m2/N , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PZT8_cE11 = 13.7e10 #14.9e10 # N*m^-2
PZT8_cE33 = 12.4e10 # 13.2e10 # N*m^-2
PZT8_cE44 = 3.13e10 # N*m^-2
PZT8_cE66 = 3.4e10 # N*m^-2
PZT8_cE12 = 6.97e10 # 8.11e10 # N*m^-2
PZT8_cE13 = 7.16e10 # 8.11e10 # N*m^-2
#PZT8_cE66 =  0.5*(PIC181_cE11-PIC181_cE12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PZT8_cD11 = 14.e10 # 15.2e10 # N*m^-2
PZT8_cD33 = 16.1e10 # 16.9e10 # N*m^-2
PZT8_cD44 = 4.46e10 # N*m^-2
PZT8_cD12 = 7.28e10 # 8.41e10 # N*m^-2
PZT8_cD13 = 6.08e10 # 7.03e10 # N*m^-2
PZT8_cD66 = 0.5*(PIC181_cD11-PIC181_cD12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PZT8_Qm = 1000.
PZT8_Qe = 250.
PZT8_N1 = 1700.
PZT8_N3t = 2070.
PZT8_N3a = 2000.

PZT8_VD3 = 4720.
PZT8_VD4 = 4720.
PZT8_VE4 = 4720.


PZT8 = Piezoceramic(ρ=PZT8_ρ, Tc=PZT8_Tc, eT33pe0=PZT8_eT33pe0, eT33=PZT8_eT33, eT11pe0=PZT8_eT11pe0, eT11=PZT8_eT11, tan_delta= PZT8_tan_delta,
          kp=PZT8_kp, kt=PZT8_kt, k31=PZT8_k31, k33=PZT8_k33, k15=PZT8_k15,
          d31=PZT8_d31, d33=PZT8_d33, d15=PZT8_d15,
          e31=PZT8_e31, e33=PZT8_e33, e15=PZT8_e15,
          g31=PZT8_g31, g33=PZT8_g33,
          N1=PZT8_N1, N3=PZT8_N3a, Nt=PZT8_N3t,
          SE11=PZT8_SE11, SE12=PZT8_SE12, SE13=PZT8_SE13, SE33=PZT8_SE33, SE44=PZT8_SE44, SE66=PZT8_SE66,
          SD11=PZT8_SD11, SD12=PZT8_SD12, SD13=PZT8_SD13, SD33=PZT8_SD33, SD44=PZT8_SD44, SD66=PZT8_SD66,
          cE11=PZT8_cE11, cE12=PZT8_cE12, cE13=PZT8_cE13, cE33=PZT8_cE33, cE44=PZT8_cE44, cE66=PZT8_cE66,
          cD11=PZT8_cD11, cD12=PZT8_cD12, cD13=PZT8_cD13, cD33=PZT8_cD33, cD44=PZT8_cD44, cD66=PZT8_cD66,
          Qm=PZT8_Qm,
          eS11=PZT8_eS11, eS33=PZT8_eS33,
          h33=PZT8_h33
)



# PIC255 - ceramika mieka
# http://fisica.cab.cnea.gov.ar/bt/images/d/d3/PICat.pdf
PIC255_ρ = 7800. # kg/m3
PIC255_Tc = 350. # C
PIC255_eT33pe0 = 1750.
PIC255_eT33 = PIC255_eT33pe0*e0
PIC255_eT11pe0 = 1650.
PIC255_eT11 = PIC255_eT11pe0*e0
PIC255_e_e = 1650.
PIC255_tan_delta = 20.e-3

PIC255_kp = 0.62
PIC255_kt = 0.47
PIC255_k31 = 0.35
PIC255_k33 = 0.69
PIC255_k15 = 0.66
PIC255_d31 = -180.e-12 # C/N
PIC255_d33 = 400.e-12 # C/N
PIC255_d15 = 550.e-12 # C/N
PIC255_g31 = -11.3e-3 # Vm/N
PIC255_g33 = 25.e-3 # Vm/N

PIC255_Np = 2000. # Hzm
PIC255_N1 = 1420. # Hzm
#PIC255_N3 = 0. # Hzm
PIC255_Nt = 2000. # Hzm
PIC255_SE11 = 16.1e-12 # m2/N
PIC255_SE33 = 20.7e-12 # m2/N
PIC255_cD33 = 10.e10 # N/m2 # przyjecte z PIC151, dla PIC255 brak danych
PIC255_Qm = 80.
PIC255_TK_S33 = 4.e-3 # 1/K

PIC255_C = -1.0 # %
PIC255_CK = -1.0 # %

# http://www.sea-acustica.es/fileadmin/Oporto16/189.pdf
PIC255_cE11 = 1.23e11 # N*m^-2
PIC255_cE12 = 2.226e10 # N*m^-2
PIC255_cE13 = 7.67e10 # N*m^-2
PIC255_cE33 = 9.71e10 # N*m^-2
PIC255_cE44 = 7.025 # N*m^-2
PIC255_cE66 =  0.5*(PIC181_cE11-PIC181_cE12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł

PIC255_cD33 = PIC255_cE33/(1.-PIC255_kt**2) #http://ira.lib.polyu.edu.hk/bitstream/10397/566/1/material-parameters_97.pdf

# http://www.sea-acustica.es/fileadmin/Oporto16/189.pdf
PIC255_e31 = -7.15 # C/m^2
PIC255_e33 = 13.7 # C/m^2
PIC255_e15 = 11.9 # C/m^2

PIC255_c = 4485. # m/s, http://www.sea-acustica.es/fileadmin/Oporto16/189.pdf
PIC255_Z = PIC255_ρ*PIC255_c # Rayl




PIC255 = Piezoceramic(name = "",ρ=PIC255_ρ, Tc=PIC255_Tc, eT33pe0=PIC255_eT33pe0, eT33=PIC255_eT33, eT11pe0=PIC255_eT11pe0, eT11=PIC255_eT11, e_e=PIC255_e_e, tan_delta= PIC255_tan_delta,
          kp=PIC255_kp, kt=PIC255_kt, k31=PIC255_k31, k33=PIC255_k33, k15=PIC255_k15,
          d31=PIC255_d31, d33=PIC255_d33, d15=PIC255_d15,
          e31=PIC255_e31, e33=PIC255_e33, e15=PIC255_e15,
          g31=PIC255_g31, g33=PIC255_g33,
          Np=PIC255_Np, N1=PIC255_N1, Nt=PIC255_Nt,
          SE11=PIC255_SE11, SE33=PIC255_SE33,
          cE11=PIC255_cE11, cE12=PIC255_cE12, cE13=PIC255_cE13, cE33=PIC255_cE33, cE44=PIC255_cE44, cE66=PIC255_cE66,
          cD33=PIC255_cD33,
          Qm=PIC255_Qm, TK_S33=PIC255_TK_S33
)



# NCE80 - ceramika twarda
# http://piezomat.org/materials/86

NCE80_ρ = 7800. # kg/m3
NCE80_Tc = 305. # *C
NCE80_d15 = 380.0e-12 # C/N
NCE80_d31 = -100.0e-12 # C/N
NCE80_d33 = 240.0e-12 # C/N
NCE80_g31 = -11.0e-3 # Vm/N
NCE80_g33 = 27.0e-3 # Vm/N
NCE80_tan_delta = 0.020
NCE80_eT33 = 1000.
NCE80_eT11 = 1050.
NCE80_kt = 0.48
NCE80_kp = 0.55
NCE80_k31 = 0.30
NCE80_k33 = 0.68
NCE80_SE11 = 11.e-12 # m2/N
NCE80_SE12 = -3.7e-12 # m2/N
NCE80_SE13 = -4.8e-12 # m2/N
NCE80_SE33 = 14.e-12 # m2/N
NCE80_SE44 = 31.9e-12 # m2/N
NCE80_SE66 = 30.4e-12 # m2/N
NCE80_Qm = 1000.
NCE80_N1 = 1610. # Hzm
NCE80_N3 = 1500. # Hzm
NCE80_Np = 2270. # Hzm
NCE80_Nt = 2050. # Hzm

NCE80 = Piezoceramic(ρ= NCE80_ρ, Tc= NCE80_Tc, eT33= NCE80_eT33, eT11= NCE80_eT11, tan_delta= NCE80_tan_delta,
          kp=NCE80_kp, kt=NCE80_kt, k31=NCE80_k31, k33=NCE80_k33, 
          d15=NCE80_d15, d31=NCE80_d31, d33=NCE80_d33, 
          g31=NCE80_g31, g33=NCE80_g33,
          Np=NCE80_Np, N1=NCE80_N1, N3=NCE80_N3, Nt=NCE80_Nt,
          SE11=NCE80_SE11, SE12=NCE80_SE12, SE13=NCE80_SE13, SE33=NCE80_SE33, SE44=NCE80_SE44, SE66=NCE80_SE66,
          Qm=NCE80_Qm, 
)




# PZT-5H - prawdopodobnie ceramika twarda
# http://piezomat.org/materials/128
# "A new topological structure for the Langevin-type ultrasonic transducer"; Xiaolong Lu, Junhui Hu, Hanmin Peng, Yuan Wang
# Material class: Lead Zirconate Titanate
# https://pl.scribd.com/doc/34864965/PZT-Material-Properties


PZT_5H_damping_coefficient = 5.5e-4 # "A new topological structure for the Langevin-type ultrasonic transducer"
PZT_5H_cp = 491. # Heat capacity, J/(kg*K)
PZT_5H_λ = 1.5 # Thermal conductivity, W / (m * K)
PZT_5H_ρ = 7600. # 7500. # kg/m3 - "A new topological structure for the Langevin-type ultrasonic transducer"
PZT_5H_Tc = 195. # *C
PZT_5H_d31 = -265e-12 # m/V # http://www.bostonpiezooptics.com/ceramic-materials-pzt -280.0e-12 # C/N
PZT_5H_d33 = 585e-12 # m/V # http://www.bostonpiezooptics.com/ceramic-materials-pzt http://www.bostonpiezooptics.com/ceramic-materials-pz 593.e-12 # 690.0e-12 # C/N - "A new topological structure for the Langevin-type ultrasonic transducer"
PZT_5H_d15 = 730e-12 # m/V # http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_g31 = -8.5e-3 # # http://www.bostonpiezooptics.com/ceramic-materials-pzt #-9.50e-3 # Vm/N
PZT_5H_g15 = 29e-3 # # http://www.bostonpiezooptics.com/ceramic-materials-pzt #-9.50e-3 # Vm/N
PZT_5H_g33 = 19.7e-3 # Vm/N
PZT_5H_tan_delta = 0.020
PZT_5H_eT33 = 3300.
PZT_5H_kt = 0.51
PZT_5H_kp = 0.65 # 0.69 # http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_k31 = 0.39 # http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_k33 = 0.75
PZT_5H_k15 = 0.680
PZT_5H_YE11 = 6.4e10 # 6.e10 # modul Younga, N/m2 #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_YE33 = 5.e10 # modul Younga, N/m2
PZT_5H_Qm = 65. #70. #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_Qe = 40. #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_Np = 1950. # Hzm
PZT_5H_Nt = 1765. #2000. # Hz-m #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_Ns = 1100. # Hz-m #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_Nr = 1950. # Hz-m #http://www.bostonpiezooptics.com/ceramic-materials-pzt
PZT_5H_ν = 0.31 # liczba Poissona
PZT_5H_c11 = 12.7e10 # N/m2

# http://onlinelibrary.wiley.com/doi/10.1002/9781119991151.app5/pdf
PZT_5H_SE11 = 16.5e-12 # m2/N
PZT_5H_SE12 = -4.78e-12 # m2/N
PZT_5H_SE13 = -8.45e-12 # m2/N
PZT_5H_SE33 = 20.7e-12 # m2/N
PZT_5H_SE55 = 43.5e-12 # m2/N
PZT_5H_SE66 = 42.6e-12 # m2/N

PZT_5H_d31 = -274e-12 # m/V
PZT_5H_d33 = 593e-12 # m/V
PZT_5H_d15 = 741e-12 # m/V

PZT_5H_eT11pe0 = 3130.
PZT_5H_eT33pe0 = 3400.

PZT_5H_eT11 = PZT_5H_eT11pe0*e0
PZT_5H_eT33 = PZT_5H_eT11*e0

# porownane z Timoshenko beam theory
PZT_5H_cE11 = 60.6e6 # Pa
PZT_5H_e31 = -16.6 # C/m2
PZT_5H_eS33 = 25.55e-9 # F/m
PZT_5H_cE55 = 23.0e6 # Pa

# Kirchhoff plate
PZT_5H_cE12 = 19.2e6 # Pa
PZT_5H_cE66 = 23.5e6 # Pa
PZT_5H_e31 = -23.4 # C/m2

PZT_5H_KT33 = 3400.
PZT_5H_KT11 = 3130.
PZT_5H_KS11 = 1700.

PZT_5H_e33 = 23.3 # C/m https://pl.scribd.com/doc/34864965/PZT-Material-Properties
PZT_5H_e15 = 17.0 # C/m2

PZT_5H_h33 = 1.8e9 # V/m #PZT_5H_e33/PZT_5H_eS33
PZT_5H_h31 = -0.505e9 # V/m
PZT_5H_h15 = 1.13e9 # V/m

PZT_5H_cD11 = 130e9 # Pa
PZT_5H_cD33 = 157e9 # Pa
PZT_5H_cD12 = 82.8e9 # Pa
PZT_5H_cD13 = 72.2e9 # Pa
PZT_5H_cD44 = 42.2e9 # Pa
PZT_5H_cD66 = 0.5*(PZT_5H_cD11-PZT_5H_cD12) # N*m^-2 , Dynamics of Mechatronics Systems: Modeling, Simulation, Control ... Autorzy Jan Awrejcewicz,Donat Lewandowski,Paweł


PZT_5H = Piezoceramic(name = "",damping_coefficient= PZT_5H_damping_coefficient,
          cp=PZT_5H_cp, λ=PZT_5H_λ, ν=PZT_5H_ν, #α3=PIC181_α3, α1=PIC181_α1,\
          ρ= PZT_5H_ρ, Tc= PZT_5H_Tc, eT33= PZT_5H_eT33, tan_delta= PZT_5H_tan_delta,
          kp=PZT_5H_kp, kt=PZT_5H_kt, k31=PZT_5H_k31, k33=PZT_5H_k33, #k15=PZT_5H_k15,
          d31=PZT_5H_d31, d33=PZT_5H_d33, #d15=PZT_5H_d15,
          e31=PZT_5H_e31, e33=PZT_5H_e33, #e15=PZT_5H_e15,
          g31=PZT_5H_g31, g33=PZT_5H_g33,
          Np=PZT_5H_Np, Nt=PZT_5H_Nt, #N1=PZT_5H_N1, N3=PZT_5H_N3, 
          SE11=PZT_5H_SE11, SE33=PZT_5H_SE33, #SE12=PZT_5H_SE12, SE13=PZT_5H_SE13, SE44=PZT_5H_SE44, SE66=PZT_5H_SE66,
          #SD11=PZT_5H_SD11, SD12=PZT_5H_SD12, SD13=PZT_5H_SD13, SD33=PZT_5H_SD33, SD44=PZT_5H_SD44, SD66=PZT_5H_SD66,
          cE11=PZT_5H_cE11, cE12=PZT_5H_cE12, cE66=PZT_5H_cE66, cE55=PZT_5H_cE55, #cE13=PZT_5H_cE13, cE33=PZT_5H_cE33, cE44=PZT_5H_cE44,
          cD11=PZT_5H_cD11, cD12=PZT_5H_cD12, cD13=PZT_5H_cD13, cD33=PZT_5H_cD33, cD44=PZT_5H_cD44, cD66=PZT_5H_cD66,
          Qm=PZT_5H_Qm, #TK_S33=PZT_5H_TK_S33,
          #eS11=PZT_5H_eS11,
          eS33=PZT_5H_eS33,
          h33=PZT_5H_h33,
          YE11= PZT_5H_YE11, YE33= PZT_5H_YE33
)



# quartz

quartz = Piezoceramic(name = "",d33 = 2.3e-12, # m/V
                        h33 = 4.9e9 # V/m
)

# pzt:
pzt = Piezoceramic(
        name = "",
        ρ = 7475.,   #kg/m3
        c = 3803., # m/s
)


# PZT4
# A   simple   electric   circuit   for   teaching   one - dimensional   characterization   of   piezoelectric plates
# F. J. Arnold, L. L. Bravo - Roger, M. S. Gonçalves, M. Grilo
PZT4_vz = 4600. # m/s
PZT4_Z0 = 34.5e6 # kg/m2s
PZT4_eS33 = 12.e-9 # F/m
PZT4_cD33 = 15.9e10 # N/m2
PZT4_h33 = 26.8e8 # N/c
PZT4_kt = 0.41
PZT4_CT0 = 3.15e-9 # F
PZT4_C0 = 2.62e-9 # F
PZT4_ρ = 7600. # kg/m3 jako PZT401

PZT4 = Piezoceramic(name = "",vz = PZT4_vz,
        Z0 = PZT4_Z0,
        eS33 = PZT4_eS33,
        cD33 = PZT4_cD33,
        h33 = PZT4_h33,
        kt = PZT4_kt,
        CT0 = PZT4_CT0,
        C0 = PZT4_C0,
        ρ = PZT4_ρ
)



###########################################################################################
#                functions
###########################################################################################


#---------------------------------------------------------------------------------- 





        
                                                                                                                                             
#----------------------------------------------------------------------------------


###########################################################################################
#                main
###########################################################################################

if __name__ == "__main__":
        pass
        #duralumin_PA6.Z_shaft_Two_port_network_plot(f=np.arange(15e3,100e3,1), l=60.0e-3, D=51e-3, useMori=False, use_speed_of_sound_approximation=True)

        """f, Z1, Z2, Z3 = duralumin_PA6.Z1_Z2_Z3_shaft_Two_port_network(f=np.arange(15e3,100e3,1), D=51e-3, d=0, l=60e-3, useMori=False, use_speed_of_sound_approximation=True)

        pyplot.semilogy(f*1e-3, np.absolute(Z1), label="Z1, Z2")
        pyplot.semilogy(f*1e-3, np.absolute(Z3), label="Z3")

        pyplot.legend(loc=1)

        pyplot.xlabel("f, kHz")
        pyplot.ylabel("|Z|")
        
        pyplot.show()"""

        #print( PIC181.external_force(20e3,A,5600,5600,5e-3,10.) )

        #print( PIC181.calculation_of_ceramic_parameters( fs=419.695e3, fp=454e3, D=18.98e-3*2, d=7.48e-3*2, l=5e-3, C0=1.74e-9, m=35.9e-3) )

        #f_c2, Z_c2, O_c2 = np.loadtxt(r"D:\Tomasz\dokumentacje\ceramiki\PIC181\pic181 - nr4 - D50 d19.6 t5\ceramika 2\skan_Hz.dat", unpack=True)
        """f_c2, Z_c2, O_c2 = np.loadtxt(r"D:\Tomasz\dokumentacje\przetworniki\mocy\przetwornik do papy papierowej\ceramiki PIC181 D50d19.5L5\ceramika9\skan_Hz.dat", unpack=True)

        f0, Z0 = PIC181.electrical_impedance_axial_vibrations_Kogut_round_ceramics(f=np.arange(400e3,500e3,10))

        print( PIC181.calculation_of_ceramic_parameters( fs=415543, fp=427374, D=49.89e-3, d=19.54e-3, l=5e-3, C0=3.64e-9, m=65.425e-3) )

        f1, Z1 = PIC181.electrical_impedance_axial_vibrations_Kogut_round_ceramics(f=np.arange(400e3,500e3,10))

        Zc = PIC181.model_Mason(f=f0, D=50e-3, d=19e-3, l=5e-3, p=1)[0]

        #pyplot.semilogy(f0*1e-3, np.absolute(Z0) )
        pyplot.semilogy(f1*1e-3, np.absolute(Z1) )
        pyplot.semilogy(f0*1e-3, np.absolute(Zc) )

        pyplot.semilogy(f_c2*1e-3, Z_c2 )
        
        #f_r, Z_r = PIC181.electrical_impedance_radial_vibrations_Kogut_round_ceramics()
        #pyplot.semilogy(f_r*1e-3, np.absolute(Z_r) )

        pyplot.xlabel("f, kHz")
        pyplot.ylabel("|Z|, $\Omega$")
        
        pyplot.show()"""

        

        

        """f_r, Z_r = PZT8.electrical_impedance_radial_vibrations_Kogut_round_ceramics(f=np.arange(1e3,450e3,1),D=50e-3, d=20e-3, l=6e-3)
        pyplot.semilogy(f_r*1e-3, np.absolute(Z_r) )

        pyplot.xlabel("f, kHz")
        pyplot.ylabel("|Z|, $\Omega$")
        
        pyplot.show()"""

        #duralumin_PA6.Z_cone_plot(f=np.arange(15e3,60e3,1), D1=51e-3, D2=40e-3, l=40e-3, useMori=True)
        #duralumin_PA6.Z_shaft_plot(f=np.arange(15e3,45e3,1), D=50.9e-3, l=121.1e-3, useMori=True)
        #duralumin_PA6.Z_shaft_Two_port_network_plot(f=np.arange(15e3,45e3,1), D=50.9e-3, l=121.1e-3, useMori=True)
        #duralumin_PA6.Z_shaft_plot(f=np.arange(15e3,50e3,1), D=70.0e-3, l=34.87e-3, useMori=True)
        #duralumin_PA6.Z_shaft_plot(f=np.arange(15e3,50e3,1), D=51.0e-3, l=50e-3, useMori=True)
        #duralumin_PA6.Z_shaft_Two_port_network_plot(f=np.arange(15e3,50e3,1), D=51e-3, l=50e-3, useMori=True)

        #PIC181.model_Mason_plot(f=np.arange(15e3,600e3,100), p=1)
        #PIC181.electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot(f=np.arange(15e3,1e6,100))

        #assert round(duralumin_PA6.resonant_frequency_shaft_fi50(40e-3),2) == round(duralumin_PA6.resonant_frequency_shaft(40e-3,51e-3),2)

        #print( duralumin_PA6.resonant_frequency_shaft(50e-3, 121.1e-3) )

        #f, Z = duralumin_PA6.Z_shaft_plot(f=np.arange(15e3,80e3,1), D=50e-3, l=29.6e-3, useMori=True)
        #print( duralumin_PA6.resonant_frequency_shaft(D=90.9e-3, l=40.9e-3) )

        #steel_1_2379.Z_shaft_Two_port_network_plot(f=np.arange(15e3,50e3,1), D=51.0e-3, l=51.9e-3, useMori=True)
        #duralumin_PA6.Z_shaft_Two_port_network_plot(f=np.arange(15e3,50e3,1), D=51.0e-3, l=60e-3, useMori=True, Z1_load=liquids.water.Z)
        #duralumin_PA6.Z_shaft_plus_cone_Two_port_network_plot()
        #duralumin_PA6.Z_cone_plot(useMori=True)
        #duralumin_PA6.Z_cone_Two_port_networke_plot(useMori=True)

        #f_c, Z_c, O_c = np.loadtxt(r"D:\Tomasz\dokumentacje\ceramiki\PIC181\ceramika z kablem przylutowanym do elektrody\zanurzony w wodzie\skan_Hz.dat", unpack=True)
        #PIC181.calculation_of_ceramic_parameters( fs=418426, fp=465294, D=49.89e-3, d=19.54e-3, l=5e-3, C0=3.64e-9, m=65.425e-3)
        #Zc = PIC181.model_Mason(f=f_c, D=50e-3, d=19.6e-3, l=5e-3, p=1, Z1_load=liquids.water, Z2_load=liquids.water)[0]

        #a,b = np.polyfit([f_c[0], f_c[int(len(f_c)/3)],f_c[-1]], [ 693.500739  ,   46.64904784,   27.48263197], 1)
        #y = a*f_c + b

        #pyplot.semilogy(f_c, np.absolute(Z_c*np.exp(1j*O_c)))
        #pyplot.semilogy(f_c, np.absolute(Zc))
        #pyplot.semilogy(f_c, np.absolute(y))

        #pyplot.show()

        #PIC255.model_Mason_plot(f=np.arange(700e3,1.2e6,1),D=9.95e-3, d=0, l=2.0e-3,p=1)
        #print( PIC255.calculation_of_ceramic_parameters( fs=1055000, fp=1118750, D=9.95e-3, d=0e-3, l=2e-3, C0=.530e-9, m=1.198e-3) )
        #PIC255.model_Mason_plot(f=np.arange(700e3,1.2e6,1),D=9.95e-3, d=0, l=2.0e-3,p=1)

        #PZT8.model_Mason_plot(f=np.arange(10e3,450e3,10), D=50e-3,d=20e-3,l=6e-3)
        #SonoxP8.model_Mason_plot(f=np.arange(10e3,500e3,10), D=38e-3,d=15e-3,l=5e-3)

        #f=np.arange(40e3,50e3,10)
        #f, Ze_czwornik, Ze11, Ze12, Ze13 = duralumin_PA6.Z_shaft_plus_cone_Two_port_networke(f=f, D=51.e-3, d=0, l=20e-3, D1=51e-3, D2=69.7e-3, L=47.1e-3, useMori=True, Z1_load=0, Z2_load=0)
        #pyplot.semilogy(f, np.absolute(Ze_czwornik))
        #pyplot.show()

        #duralumin_PA6.Z_shaft_plus_cone_Two_port_network_plot(f=np.arange(15e3,50e3,1), D=51e-3, d=0e-3, l=20e-3, D1=51e-3, D2=70e-3, L=47.1e-3, useMori=True, use_speed_of_sound_approximation=False, Z1_load=0, Z2_load=0)
        

        
        #PIC181.model_Mason_plot(f=np.arange(390e3,480e3,10), D=50e-3, d=9.8e-3, l=5e-3)
        #PIC181.electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot(f=np.arange(390e3,480e3,100),D=50e-3, l=5e-3)
        #PIC181.model_Mason_plot(f=np.arange(200e3,4e6,10), D=25e-3, d=9e-3, l=1.9e-3)
        #PIC181.electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot(f=np.arange(350e3,1600e3,100),D=25e-3,l=1.9e-3)
        #SonoxP8.electrical_impedance_axial_vibrations_Kogut_round_ceramics_plot(f=np.arange(15e3,500e3,100),D=38e-3, l=5e-3)

        # rura

        """material = Material(
                        ρ = 7900.,      # density, kg/m3
                        E = 2.e11,   # adopted, Young module, Pa
                        ν = 0.28,       # Poisson
                )

        D = 30e-3
        h = 3e-3
        L = 300e-3

        print( steel_316.cylindrical_shell__radial_vibration_with_axial_symetry(D=D, h=h, L=L ) )
        print( steel_316.cylindrical_shell__vibration( D=D, h=h, L=L ) )
        print( steel_316.cylindrical_shell__vibration( D=D, h=h, L=L, m=1, n=1 ) )
        print( steel_316.cylindrical_shell__vibration( D=D, h=h, L=L, m=1, n=0 ) )
        print( steel_316.cylindrical_shell__vibration( D=D, h=h, L=L, m=0, n=1 ) )"""
        
        
        #D = [70e-3,]
        #L = 1e-3*np.array([45,])

        #for DD,LL in zip(D,L):
                #print(DD*1e3, LL*1e3)

                #f = duralumin_PA6.resonant_frequency_shaft(D=DD, l=LL, skok=1)
                #duralumin_PA6.c = f*2.*LL
                #print(duralumin_PA6.c)
                        
                #f, Z, Z1, Z2, Z3 = duralumin_PA6.Z_shaft_Two_port_network_plot(f=np.arange(15e3,60e3,10), D=DD, l=LL, useMori=False)
                #print(duralumin_PA6.c)
                #print("%.3f" % duralumin_PA6.resonant_frequency_shaft(D=DD, l=LL, skok=1))

        #D = 10.15e-3
        #L = 100.9e-3
        """D = 3.e-3
        L = 21.75e-3
        print(macor.c)
        f = macor.resonant_frequency_shaft(D=D, l=L, skok=1)
        macor.c = f*2.*L
        print(macor.c)"""


        #duralumin_PA6.Z_shaft_plot(f=np.arange(90e3,100e3,1), D=40e-3, d=0, l=26.5*2e-3, useMori=False)
        #duralumin_PA6.Z_shaft_plot(f=np.arange(90e3,120e3,1), D=40e-3, d=0, l=49e-3, useMori=False)
        #duralumin_PA6.Z_shaft_plot(f=np.arange(10e3,30e3,1), D=55e-3, d=0, l=112.0e-3, useMori=True)

        #duralumin_PA6.sonotroda_schodkowa( )

        #duralumin_PA6.Z_cone_plot(f=np.arange(15e3,30e3,1), D1=51e-3, D2=40e-3, l=123.7e-3, useMori=False)
        #duralumin_PA6.Z_cone_Two_port_networke_plot(f=np.arange(15e3,30e3,1), D1=51e-3, D2=40e-3, l=123.7e-3, useMori=False)
        
