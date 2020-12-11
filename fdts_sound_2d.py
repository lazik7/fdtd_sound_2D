# -*- coding: utf-8 -*-


from FreeCAD import Gui
from PySide import QtGui, QtCore
import csv

from FreeCAD import Gui
from FreeCAD import Base
import Fem
import FreeCAD, FreeCADGui, Part, os, math

from ctypes import *
import numpy as np
from matplotlib.pyplot import figure, draw, pause
from matplotlib import pyplot as plt
import datetime

__dir__ = os.path.dirname(__file__)
iconPath = os.path.join( __dir__, 'Icons' )


fdtd_dx = 0.1
fdtd_cmax = 8000.
fdtd_time_max = 10

fdtd_x = 20
fdtd_y = 20

fdtd_material_number = 0

fdtd_materials = {"water": [1000,2.25,2.25,2.25,0],"steel": [7800,268.5,104.4,268.2,82]}

fdtd_material0 = "water"


Carte = np.zeros(int((fdtd_x/fdtd_dx)*(fdtd_y/fdtd_dx)), dtype=np.intc)
CarteGen = [20100,]
        
TableLegerete = (1.0/(fdtd_materials[fdtd_material0][0]/1000.))*np.ones(100, dtype=np.double)
TableC11 = fdtd_materials[fdtd_material0][1]*np.ones(100, dtype=np.double)
TableC12 = fdtd_materials[fdtd_material0][2]*np.ones(100, dtype=np.double)
TableC22 = fdtd_materials[fdtd_material0][3]*np.ones(100, dtype=np.double)
TableC33 = fdtd_materials[fdtd_material0][4]*np.ones(100, dtype=np.double)

fdtd_plot_max = 0.5




class Time_max():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'time.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+x", 
                'MenuText': "dx",
                'ToolTip' : "dx"}

    def Activated(self):

        global fdtd_time_max
        
        fdtd_time_max_old = fdtd_time_max
        fdtd_time_max = 0

        while fdtd_time_max==0:
            try:
                fdtd_time_max = float(str(QtGui.QInputDialog.getText(None, "Get text", "fdtd_time_max, (old=%f) = " % fdtd_time_max_old)[0]))
            except:
                pass

        
    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True




class plotVmax():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'plotVmax.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+x", 
                'MenuText': "dx",
                'ToolTip' : "dx"}

    def Activated(self):

        global fdtd_plot_max
        
        fdtd_plot_max_old = fdtd_plot_max
        fdtd_plot_max = 0

        while fdtd_plot_max==0:
            try:
                fdtd_plot_max = float(str(QtGui.QInputDialog.getText(None, "Get text", "fdtd_plot_max, (old=%f) = " % fdtd_plot_max_old)[0]))
            except:
                pass

        
    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True






class MakeGen():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'MakeGen.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+a", 
                'MenuText': "Area",
                'ToolTip' : "Area"}

    def Activated(self):

        global fdtd_x
        global fdtd_y
        global fdtd_dx

        global CarteGen

        CarteGen = []
        
        dxp2 = fdtd_dx/2.

        nx = int(fdtd_x/fdtd_dx)
        ny = int(fdtd_y/fdtd_dx)

        point_list = np.zeros([nx, ny])

        volume = Gui.Selection.getSelectionEx()[0]

        for i in range(nx):
            for j in range(ny):
                point=FreeCAD.Vector(i*fdtd_dx+dxp2,j*fdtd_dx+dxp2,0)

                if volume.Object.Shape.isInside(point, fdtd_dx/10, False):
                    CarteGen.append( int(ny*i+j) )

            print( "%.1f" % (i/nx*100) )
 
                
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True






class MakeMesh():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'MakeMesh.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+a", 
                'MenuText': "Area",
                'ToolTip' : "Area"}

    def Activated(self):

        global fdtd_x
        global fdtd_y
        global fdtd_dx

        global TableLegerete
        global TableC11
        global TableC12
        global TableC22
        global TableC33

        global Carte
        global fdtd_material_number

        fdtd_material_number += 1

        TableLegerete[fdtd_material_number] = (1.0/(fdtd_materials["steel"][0]/1000.))
        TableC11[fdtd_material_number] = fdtd_materials["steel"][1]
        TableC12[fdtd_material_number] = fdtd_materials["steel"][2]
        TableC22[fdtd_material_number] = fdtd_materials["steel"][3]
        TableC33[fdtd_material_number] = fdtd_materials["steel"][4]

        dxp2 = fdtd_dx/2.

        nx = int(fdtd_x/fdtd_dx)
        ny = int(fdtd_y/fdtd_dx)
        

        volume = Gui.Selection.getSelectionEx()[0]

        print("len(c) = %s" % len(Carte) )

        #np.arange(0,fdtd_y,fdtd_dx)
        #np.arange(0,fdtd_x,fdtd_dx)

        for i in range(nx):
            for j in range(ny):
                point=FreeCAD.Vector(i*fdtd_dx+dxp2,j*fdtd_dx+dxp2,0)

                if volume.Object.Shape.isInside(point, fdtd_dx/10, False):
                    try:
                        Carte[ int(ny*i+j) ] = fdtd_material_number
                    except:
                        print(i,j, nx*i+j, len(Carte), nx, ny )

            print( "%.1f" % (i/nx*100) )

            

        
 
                
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True





class Area():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'Area.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+a", 
                'MenuText': "Area",
                'ToolTip' : "Area"}

    def Activated(self):

        global fdtd_x
        global fdtd_y
        global Carte
        global fdtd_dx

        
        x_old = fdtd_x
        fdtd_x = 0

        while fdtd_x==0:
            try:
                fdtd_x = float(str(QtGui.QInputDialog.getText(None, "Get text", "x, mm, (old=%f) = " % x_old)[0]))
            except:
                pass


        
        y_old = fdtd_y
        fdtd_y = 0

        while fdtd_y==0:
            try:
                fdtd_y = float(str(QtGui.QInputDialog.getText(None, "Get text", "y, mm, (old=%f) = " % y_old)[0]))
            except:
                pass

        m = Fem.FemMesh()
        m.addNode(0,0,0)
        m.addNode(fdtd_x,0,0)
        m.addNode(fdtd_x,fdtd_y,0)
        m.addNode(0,fdtd_y,0)

        m.addVolume([1,2,3,4])

        Fem.show(m)

        nx = int(fdtd_x/fdtd_dx)
        ny = int(fdtd_y/fdtd_dx)

        Carte = np.zeros(int(nx*ny), dtype=np.intc)
        fdtd_material_number = 0
        
 
                
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True







try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_DockWidget(object):
    def setupUi(self, DockWidget, setRowCount=3, setColumnCount=6):
        DockWidget.setObjectName(_fromUtf8("DockWidget"))
        DockWidget.resize(267, 136)
        DockWidget.setFloating(True)
        self.dockWidgetContents = QtGui.QWidget()
        self.dockWidgetContents.setObjectName(_fromUtf8("dockWidgetContents"))
        self.gridLayout = QtGui.QGridLayout(self.dockWidgetContents)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))


        self.table = QtGui.QTableWidget()
        self.table.setRowCount(setRowCount)
        self.table.setColumnCount(setColumnCount)

        button_add = QtGui.QPushButton('add')
        button_add.clicked.connect(self.add)

        button_remove = QtGui.QPushButton('remove')
        button_remove.clicked.connect(self.remove)

        button_save = QtGui.QPushButton('save')
        button_save.clicked.connect(self.save)
        
        self.vLayout = QtGui.QVBoxLayout()
        self.vLayout.addWidget(self.table)
        self.vLayout.addWidget(button_add)
        self.vLayout.addWidget(button_remove)
        self.vLayout.addWidget(button_save)
        
        
        self.gridLayout.addLayout(self.vLayout, 1, 0, 1, 1)

        DockWidget.setWidget(self.dockWidgetContents)

        self.retranslateUi(DockWidget)
        QtCore.QMetaObject.connectSlotsByName(DockWidget)

    def retranslateUi(self, DockWidget):
        DockWidget.setWindowTitle(_translate("DockWidget", "Mass calculator", None))


    def put_in_table(self):
        global fdtd_materials
        global fdtd_material0
        

        self.table.setItem(0, 0, QtGui.QTableWidgetItem("Material name"))
        self.table.setItem(0, 1, QtGui.QTableWidgetItem("q[kg/m3]"))
        self.table.setItem(0, 2, QtGui.QTableWidgetItem("C11"))
        self.table.setItem(0, 3, QtGui.QTableWidgetItem("C12"))
        self.table.setItem(0, 4, QtGui.QTableWidgetItem("C22"))
        self.table.setItem(0, 5, QtGui.QTableWidgetItem("C33"))


        for i, material in enumerate(fdtd_materials.keys()):
            self.table.setItem(i+1, 0, QtGui.QTableWidgetItem( material ))
            for j in range(len(fdtd_materials[material])):
                self.table.setItem(i+1, j+1, QtGui.QTableWidgetItem( "%s" % fdtd_materials[material][j] ))


    def save(self):
        global fdtd_materials

        for i in range(1,self.table.rowCount()):
            material = self.table.item(i, 0).text()
            for j in range(1, len(fdtd_materials[material])):
                fdtd_materials[material][j-1] = float(self.table.item(i, j).text())


    def add(self):
        self.table.setRowCount(self.table.rowCount()+1)


    def remove(self):
        self.table.setRowCount(self.table.rowCount()-1)
            

        
       



dlg = QtGui.QDockWidget()
dlg.ui = Ui_DockWidget()
dlg.ui.setupUi(dlg)
Gui.getMainWindow().addDockWidget(QtCore.Qt.RightDockWidgetArea, dlg)
dlg.setFloating(True)
dlg.hide()



class Materials():
  """Display a calculator for needed screw holes"""

  def GetResources(self):
    FreeCAD.Console.PrintLog("Getting resources\n")
    icon = os.path.join( iconPath , 'Materials.png')
    return {'Pixmap'  : icon , # the name of a svg file available in the resources
            'MenuText': "Mass calculator" ,
            'Accel' : "Shift+M",
            'ToolTip' : "Calc mass"}
 
  def Activated(self):
    if dlg.isHidden():
      dlg.show()
    else:
      dlg.hide()

    dlg.ui.put_in_table()    
      
   
  def IsActive(self):
    return True







class Cmax():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'cmax.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+c", 
                'MenuText': "dx",
                'ToolTip' : "dx"}

    def Activated(self):

        global fdtd_cmax
        c_old = fdtd_cmax
        fdtd_cmax = 0

        while fdtd_cmax==0:
            try:
                fdtd_cmax = float(str(QtGui.QInputDialog.getText(None, "Get text", "dx, mm, (old=%f) = " % c_old)[0]))
            except:
                pass         
                
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True




class Dx():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'dx.png')
        

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+x", 
                'MenuText': "dx",
                'ToolTip' : "dx"}

    def Activated(self):

        global fdtd_dx
        global Carte
        
        dx_old = fdtd_dx
        fdtd_dx = 0

        while fdtd_dx==0:
            try:
                fdtd_dx = float(str(QtGui.QInputDialog.getText(None, "Get text", "dx, mm, (old=%f) = " % dx_old)[0]))
            except:
                pass

        Carte = np.zeros(int((fdtd_x/fdtd_dx)*(fdtd_y/fdtd_dx)), dtype=np.intc)
                
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True






class Fdts_sound_2d_calc():
    """My new command"""

    from os.path import join, dirname    
    __dir__ = dirname(__file__)
    mojmod_icons_path =  join( __dir__, 'Icons')
    Icon = join( mojmod_icons_path , 'Fdts_sound_2d_calc.png')
    Lib = join( __dir__ , 'lib_fdtd_spund2D.so')
    

    def GetResources(self):
        return {'Pixmap'  : self.Icon, 
                'Accel' : "Shift+F", 
                'MenuText': "FDTD",
                'ToolTip' : "FDTD"}

    def Activated(self):

        global fdtd_dx
        global fdtd_x
        global fdtd_y
        global fdtd_cmax
        global fdtd_materials
        global fdtd_material0

        global TableLegerete
        global TableC11
        global TableC12
        global TableC22
        global TableC33

        global Carte
        global CarteGen

        global fdtd_plot_max

        global fdtd_time_max

        
                        
        fg = figure()
        ax = fg.gca()

        cmax = fdtd_cmax/1000.
        dx = fdtd_dx
        dt = 0.99*((dx/cmax)/(2**0.5))

        f = 1 # MHz
        t = fdtd_time_max

        IndiceTempsMax = int(t/dt+1)

        self.libsound = np.ctypeslib.load_library( self.Lib, self.__dir__)

        array_1d_double = np.ctypeslib.ndpointer( dtype=np.double, ndim=1, flags='CONTIGUOUS')
        array_1d_int = np.ctypeslib.ndpointer( dtype=np.intc, ndim=1, flags='CONTIGUOUS')

        self.libsound.calc.restype = None
        self.libsound.calc.argtypes = [c_double, c_int, c_int, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, array_1d_double, c_int, c_double, array_1d_int]


        DeltaTsurH = dt/dx
        X = int(fdtd_x/fdtd_dx)
        Y = int(fdtd_y/fdtd_dx)

        Vx = np.zeros((X+1)*Y, dtype=np.double)
        Vy = np.zeros(X*(Y+1), dtype=np.double)
        V = np.zeros(X*Y, dtype=np.double)
        Txx = np.zeros(X*Y, dtype=np.double)
        Txy = np.zeros(X*Y, dtype=np.double)
        Tyy = np.zeros((X+1)*(Y+1), dtype=np.double)

        

        
        #-------------------------------------------------------- plot



        Varray = np.array(np.split(np.array(V), X))
        plot_max = fdtd_plot_max
        #plot_max = 0.0004
        h = ax.imshow(Varray.transpose(),vmin=0, vmax=plot_max,  extent=[0, X*dx, 0, Y*dx])  # set initial display dimensions
        #ax.axvline(x=15, color="red")

        #h.set_cmap("Blues_r")
        cbar = fg.colorbar(h)
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")

        draw(), pause(1e-6)

        


        for IndiceTemps in range(IndiceTempsMax):
            self.libsound.calc(DeltaTsurH, X-1, Y, TableLegerete, Vx, Vy, V, Txx, Txy, Tyy, TableC11, TableC12, TableC22, TableC33, IndiceTemps, dt, Carte)

            #Txx[int(len(Carte)/2 + self.X/2)] += np.exp(-((IndiceTemps*dt-3)/0.5)**2)*np.sin(2.*np.pi*1*IndiceTemps*dt)
            #Tyy[int(len(Carte)/2 + self.X/2)] += np.exp(-((IndiceTemps*dt-3)/0.5)**2)*np.sin(2.*np.pi*1*IndiceTemps*dt)
          
            for i in CarteGen:
                Txx[i] += np.exp(-((IndiceTemps*dt-3)/0.5)**2)*np.sin(2.*np.pi*1*IndiceTemps*dt)
                Tyy[i] += np.exp(-((IndiceTemps*dt-3)/0.5)**2)*np.sin(2.*np.pi*1*IndiceTemps*dt)

            

            if round(IndiceTemps*dt*10,1) % 1 == 0:
                Varray = np.array(np.split(np.array(V), X))
                ax.set_title("time %.1f $\mu$s" % (IndiceTemps*dt))
                h.set_data( Varray.transpose()[::-1] )
                draw(), pause(1e-6)


        

        
        

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True


