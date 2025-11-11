import os
import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from scipy import constants as const
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
# Loading console code to verify initialisation:
print("=======================")
print("Initialising program...")
print("=======================")
print("Loading variables...")
# To determine the base directory of the .ico file:
basedir = os.path.dirname(__file__)
# Define initial basic variables
constantg = const.G
g = constantg
constantc = const.c
c = constantc
Solar_Mass = 1.98847e30

def Schwarzschild_Radius(M,g,c):
    # Gives the Schwarzschild Radius for a given mass (M)
    Schwarzschild_Radius = (2*g*M/c**2) # Just Rs = 2GM/c^2
    return Schwarzschild_Radius

def Gravitational_Redshift(M,g,c,d):
    # Gives the Gravitational Redshift for a given mass (M) at a chosen distance
    Gravitational_Redshift = 1/(np.sqrt(1-(2*g*M)/(d*c**2))) - 1 # Equation of z = 1/sqrt(1- 2GM/rc^2) - 1
    return Gravitational_Redshift

def Gravitational_Time_Dilation(M,g,c,d):
    # Gives the Gravitational Redshift for a given mass (M) at a chosen distance
    Gravitational_Time_Dilation = (np.sqrt(1-(2*g*M)/(d*c**2))) # Equation of z = 1/sqrt(1- 2GM/rc^2) - 1
    return Gravitational_Time_Dilation

def Escape_Velocity(M,g,d):
    # Gives the Gravitational Redshift for a given mass (M) at a chosen distance
    Escape_Velocity = (np.sqrt((2*g*M)/(d))) # Equation of z = 1/sqrt(1- 2GM/rc^2) - 1
    return Escape_Velocity

def convert_to_metric(chosen_mass,g,c):
    # Need to convert the scale to metric -> do this by just by using the Schwarzchild Radius equation
    # Rs = 2GM/c^2 to convert known quantities into a metric distance.
    metric_distance = (2*g*Solar_Mass*chosen_mass)/c**2
    return metric_distance

def newtonian_potential(mass, r, g):
    r = np.maximum(r, 1e-8) # To prevent dividing by zero
    new_pot = - (g * mass) / r # Just V = -GM/r
    return new_pot

print("Loading Tkinter GUI...")

class window(tk.Tk):
    def __init__(self):
        super().__init__() # Keeps initial conditions set here in other methods like the UI.
        self.title("Newtonian Potential Simulator, by John Ray")
        width= self.winfo_screenwidth() 
        height= self.winfo_screenheight()
        # Sets the Tkinter window size to be fullscreen
        self.geometry("%dx%d" % (width, height))

        # Changeable parameters:

        self.chosen_mass = tk.DoubleVar(value=100) # Default chosen mass is 1000 solar masses.
        self.Hubble_constant = tk.DoubleVar(value=70) # Default H0 value is 70
        self.G = tk.DoubleVar(value=g) # Default G value is ~ 6.67e-11 Nm^2/Kg^2.
        self.C = tk.DoubleVar(value=c) # Default C Value is ~ 3e8 m/s.
        self.D = tk.DoubleVar(value=1e6) # Default value of 1,000,000 metres
        self.X = tk.IntVar(value=250)
        self.Y = tk.IntVar(value=250) # Grid spacing is set to 250 x 250 length units
        self.Z = tk.DoubleVar(value=1e-8)
        self.R = tk.DoubleVar(value=100)
        self.scale = tk.StringVar(value="Kilometres (km)")
        self.cmap = tk.StringVar(value="Black-orange")
        self.frame = ttk.Frame(self, padding = 8)
        self.frame.pack(side=tk.LEFT, fill=tk.Y,expand=True)
        self.UI() # Loads the UI function.
        print("Loading UI...")
        self.plot() # Plots the graph (repeats when replot button is pressed).
        print("Plotting graph...")
        self.UIValues()

    def UI(self):

        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        ProgramTitle = tk.Label(self.frame,text="Newtonian Potential Simulator")
        ProgramTitle.config(font=("Helvetica", 12, "bold", "underline"))
        ProgramTitle.pack()
        ProgramSubTitle = tk.Label(self.frame,text="by John Ray")
        ProgramSubTitle.config(font=("Helvetica", 10))
        ProgramSubTitle.pack()

        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        ttk.Label(self.frame, text="Value of object mass (in units of Solar Masses):").pack()
        ttk.Entry(self.frame, textvariable=self.chosen_mass).pack()
        ttk.Label(self.frame, text="Value of the Gravitational Constant (in units of N^2/Kg^2):").pack()
        ttk.Entry(self.frame, textvariable=self.G).pack()
        ttk.Label(self.frame, text="Value of the speed of light (in units of metres/second):").pack()
        ttk.Entry(self.frame, textvariable=self.C).pack()
        ttk.Label(self.frame, text="Custom distance for values (in units of metres):").pack()
        ttk.Entry(self.frame, textvariable=self.D).pack()
        
        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        ttk.Label(self.frame, text="Units of length on axis:").pack()
        ttk.OptionMenu(self.frame, self.scale, "Kilometres (km)","Kilometres (km)", "Metres (m)", "Light Years (ly)", "Astronomical Units (AU)", "Parsecs (pcs)").pack()
        ttk.Label(self.frame, text="Colour map of graph:").pack()
        ttk.OptionMenu(self.frame, self.cmap, "Orange-black","Orange-black","White-black","Purple-red","Green-grey","Blue-pink","Blue-green","Red-yellow").pack()
        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        ttk.Label(self.frame, text="Grid X-spacing:").pack()
        ttk.Entry(self.frame, textvariable=self.X).pack()
        ttk.Label(self.frame, text="Grid Y-spacing:").pack()
        ttk.Entry(self.frame, textvariable=self.Y).pack()
        ttk.Label(self.frame, text="Max X/Y scale:").pack()
        ttk.Entry(self.frame, textvariable=self.R).pack()
        ttk.Label(self.frame, text="Max Z scale:").pack()
        ttk.Entry(self.frame, textvariable=self.Z).pack()

        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        ttk.Button(self.frame, text="Update", command=lambda:(self.plot(), self.UIValues())).pack()

        ttk.Separator(self.frame).pack(fill=tk.X, pady=10)

        plot_graph = ttk.Frame(self)
        plot_graph.pack(side=tk.RIGHT)
        self.fig = plt.Figure(figsize=(10,8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_box_aspect((1,1,0.6))

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_graph)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, plot_graph)
        toolbar.update()
        self.canvas._tkcanvas.pack(fill=tk.BOTH, expand=True)

    def plot(self):
        self.ax.cla()
        CustomMass = float(self.chosen_mass.get())
        CustomG = float(self.G.get())
        global CustomC
        CustomC = float(self.C.get())
        global CustomD
        CustomX = int(self.X.get())
        CustomY = int(self.Y.get())
        CustomR = float(self.R.get())
        CustomZ = float(self.Z.get())

        Length = convert_to_metric(CustomMass, CustomG, CustomC)
        Lengthz = convert_to_metric(CustomMass, CustomG, CustomC)
        if self.scale.get().startswith("Metres"):
            divide = float(1)
        elif self.scale.get().startswith("Kilometres"):
            Length = Length /1e3
            divide = float(1e3)          
        elif self.scale.get().startswith("Light"):
            Length = Length / (CustomC*24*3600*365)
            divide = float(1/(CustomC*24*3600*365))
        elif self.scale.get().startswith("Astronomical"):
            Length = Length * (149597870700)
            divide = float(1/149597870700)
        elif self.scale.get().startswith("Parsecs"):
            Length = Length / (3.085677581e16)
            divide = float(3.085677581e16)
        else:
            Length = Length # Defaults to metres if error.
        CustomD = float(self.D.get())
        # Same for colour maps
        if self.cmap.get().startswith("Purple-red"):
            cmap2 = str("Spectral")
        elif self.cmap.get().startswith("Orange-black"):
            cmap2 = str("copper")        
        elif self.cmap.get().startswith("White-black"):
            cmap2 = str("gray")  
        elif self.cmap.get().startswith("Blue-pink"):
            cmap2 = str("cool_r")  
        elif self.cmap.get().startswith("Blue-green"):
            cmap2 = str("winter_r")
        elif self.cmap.get().startswith("Green-grey"):
            cmap2 = str("Accent_r")
        elif self.cmap.get().startswith("Red-yellow"):
            cmap2 = str("autumn_r") 
        else:
            cmap2 = str("copper")
        KiloMass = Solar_Mass * CustomMass
        MinR = 0.5 * Length
        print(f"MinR = {MinR}")
        MaxR = CustomR * Length
        print(f"MaxR = {MaxR}")
        r = np.linspace(MinR, MaxR, CustomX)         
        phi = np.linspace(0, 2*np.pi, CustomY)
        R, PHI = np.meshgrid(r, phi, indexing='xy')
        MinRZ = 0.5 * Lengthz
        MaxRZ = CustomR * Lengthz
        rZ = np.linspace(MinRZ, MaxRZ, CustomX)
        RZ, PHIZ = np.meshgrid(rZ, phi, indexing='xy')
        X = R * np.cos(PHI); Y = R * np.sin(PHI)
        Potential = newtonian_potential(KiloMass, R, CustomG)
        PotentialZ = newtonian_potential(KiloMass, RZ, CustomG)
        global Radii
        Radii = Schwarzschild_Radius(KiloMass, CustomG, CustomC)
        global Redshift
        Redshift = Gravitational_Redshift(KiloMass, CustomG, CustomC, CustomD)
        global Dilation
        Dilation = Gravitational_Time_Dilation(KiloMass, CustomG, CustomC, CustomD)
        global EscapeVel
        EscapeVel = Escape_Velocity(KiloMass, CustomG, CustomD)
        print("Schwarzschild Radius = " + f"{Radii}")
        print("Gravitational Redshift = " + f"{Redshift}")
        print("Time Dilation = " + f"{Dilation}")
        print("Escape Velocity = " + f"{EscapeVel}")
        print(CustomD)
        Z = PotentialZ * CustomZ
        # mask infs
        Z = np.where(np.isfinite(Z), Z, np.nan)
        self.ax.plot_surface(X, Y, Z, cmap=cmap2, rstride=4, cstride=4, antialiased=True, alpha=0.9, label="Curvature")
        self.ax.scatter(0, 0, 0, color='black', s=100, label='Black Hole')
        title = ("""Newtonian potential graph for a Schwarzschild Black Hole of
        mass = {:.3g} M☉, c = {:.1f} m/s, and G = {:} N*m^2/Kg^2""").format(CustomMass, CustomC, CustomG)
        self.ax.set_title(title)
        if self.scale.get().startswith("Metres"):
            self.ax.set_xlabel("Metres (m)")
            self.ax.set_ylabel("Metres (m)")

        elif self.scale.get().startswith("Kilometres"):
            self.ax.set_xlabel("Kilometres (km)")
            self.ax.set_ylabel("Kilometres (km)")
           
        elif self.scale.get().startswith("Light"):
            self.ax.set_xlabel("Light Years (ly)")
            self.ax.set_ylabel("Light Years (ly)")

        elif self.scale.get().startswith("Astronomical"):
            self.ax.set_xlabel("Astronomical Units (AU)")
            self.ax.set_ylabel("Astronomical Units (AU)")

        elif self.scale.get().startswith("Parsecs"):
            self.ax.set_xlabel("Parsecs (pcs)")
            self.ax.set_ylabel("Parsecs (pcs)")

        else:
            self.ax.set_xlabel("Metres (m)") # Again, defaults to metres if error.
            self.ax.set_ylabel("Metres (m)")

        self.ax.set_zlabel("Potential (J/Kg)")
        lim = MaxR
        self.ax.set_xlim(-lim, lim); self.ax.set_ylim(-lim, lim)
        self.canvas.draw()

    def UIValues(self):
        if hasattr(self,"UI1"):
            self.UI1.config(text="         Schwarzschild Radius = " + f"{Radii}" + " m")
        else:
            self.UI1 = ttk.Label(self.frame, text="         Schwarzschild Radius = " + f"{Radii}" + " m")
            self.UI1.pack(anchor="w")
            ttk.Separator(self.frame).pack(fill=tk.X, pady=5)
        if CustomD > Radii:
            if hasattr(self,"UI2"):
                self.UI2.config(text="         Gravitational Redshift = " + f"{Redshift}")
            else:
                self.UI2 = ttk.Label(self.frame, text="         Gravitational Redshift = " + f"{Redshift}")
                self.UI2.pack(anchor="w")
                ttk.Separator(self.frame).pack(fill=tk.X, pady=5)
            if hasattr(self,"UI3"):
                self.UI3.config(text="             Time dilation ratio = " + f"{Dilation}")
            else:
                self.UI3 = ttk.Label(self.frame, text="             Time dilation ratio = " + f"{Dilation}")
                self.UI3.pack(anchor="w")
                ttk.Separator(self.frame).pack(fill=tk.X, pady=5)
        if CustomC > EscapeVel:
            if hasattr(self,"UI4"):
                self.UI4.config(text="           Escape Velocity = " + f"{EscapeVel}" + " m/s")
            else:
                self.UI4 = ttk.Label(self.frame, text="           Escape Velocity = " + f"{EscapeVel}" + " m/s")
                self.UI4.pack(anchor="w")
                ttk.Separator(self.frame).pack(fill=tk.X, pady=5)

        else:
            self.UI2.config(text="  Cannot determine Redshift (Custom distance > Rs)")
            self.UI3.config(text="Cannot determine Time dilation (custom distance > Rs)")
            self.UI4.config(text="             Escape Velocity exceeds the speed of light")             


print("=======================")
print("    Program loaded!")
print("=======================")

if __name__ == "__main__":
    app = window()
    app.iconbitmap(os.path.join(basedir, "icon.ico"))
    app.mainloop()
print("=======================")
print("    Program ended!")
print("   Thanks for using!")
print("=======================")