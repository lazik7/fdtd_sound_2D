# -*- coding: utf-8 -*-
###################################################################################
#
#  
#  
###################################################################################




class Fdtd_sound_2D(Workbench):

    MenuText = "fdtd_sound_2D"
    ToolTip = "fdtd sound 2D"

    from os.path import join, dirname
    import moj_mod_locator
    mojmod_icons_path =  join( dirname(moj_mod_locator.__file__), 'Icons')
    Icon = join( mojmod_icons_path , 'FNLogo.svg')

    def Initialize(self):
        "This function is executed when FreeCAD starts"

        import fdts_sound_2d
        
        FreeCADGui.addCommand('fdts_sound_2d_culc', fdts_sound_2d.Fdts_sound_2d_calc())
        FreeCADGui.addCommand('dx', fdts_sound_2d.Dx())
        FreeCADGui.addCommand('cmax', fdts_sound_2d.Cmax())
        FreeCADGui.addCommand('materials', fdts_sound_2d.Materials())
        FreeCADGui.addCommand('area', fdts_sound_2d.Area())
        FreeCADGui.addCommand('makemesh', fdts_sound_2d.MakeMesh())
        FreeCADGui.addCommand('makegen', fdts_sound_2d.MakeGen())
        FreeCADGui.addCommand('plotVmax', fdts_sound_2d.plotVmax())
        FreeCADGui.addCommand('time_max', fdts_sound_2d.Time_max())
        
        self.list = ["dx", "cmax", "area", "materials", 'makemesh', "makegen", "plotVmax", "time_max","fdts_sound_2d_culc"]
        self.appendToolbar("fdts_sound_2d_tolobar",self.list)
        
    def Activated(self):
        "This function is executed when the workbench is activated"
        return

    def Deactivated(self):
        "This function is executed when the workbench is deactivated"
        return

    def ContextMenu(self, recipient):
        "This is executed whenever the user right-clicks on screen"
        # "recipient" will be either "view" or "tree"
        self.appendContextMenu("My commands",self.list) # add commands to the context menu

    def GetClassName(self): 
        # this function is mandatory if this is a full python workbench
        return "Gui::PythonWorkbench"

       
Gui.addWorkbench(Fdtd_sound_2D())









