import math
from numpy import arange
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import meshgrid
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import cm,imshow,contour,clabel,colorbar,axis,title,show
from bs4.builder._htmlparser import major

file_first = open ("reaction_force_1.200", "r")
file_second = open ("reaction_force_1.220", "r")
file_third = open ("reaction_force_1.240", "r")
file_fourth = open ("reaction_force_1.245", "r")
file_fifth = open ("reaction_force_1.250", "r")
file_sixth = open ("reaction_force_1.300", "r")
file_seventh = open ("reaction_force_1.510", "r")
file_eightth = open ("reaction_force_1.650", "r")
file_nighth = open ("reaction_force_2.0", "r")
file_tenth = open ("reaction_force", "r")
file_eleventh = open ("cycle_load_lambda2.0", "r")
file_twelevth = open ("cycle_unload_lambda2.0", "r")
file_thirteenth = open ("cycle_reload_lambda2.0", "r")

Displacement=[]
Displacement_unload=[]
Displacement_reload=[]
Force_first=[]
Force_second=[]
Force_third=[]
Force_fourth=[]
Force_fifth=[]
Force_sixth=[]
Force_seventh=[]
Force_eightth=[]
Force_nighth=[]
Force_tenth=[]
Force_eleventh=[]
Force_twelevth=[]
Force_thirteenth=[]

for line in file_first:
    splitLine=line.split(" ")
    #Displacement.append(float(splitLine[0]))
    Force_first.append(float(splitLine[1]))
    
for line in file_second:
    splitLine=line.split(" ")
    Force_second.append(float(splitLine[1]))

for line in file_third:
    splitLine=line.split(" ")
    Force_third.append(float(splitLine[1]))

for line in file_fourth:
    splitLine=line.split(" ")
    Force_fourth.append(float(splitLine[1]))
    
for line in file_fifth:
    splitLine=line.split(" ")
    Force_fifth.append(float(splitLine[1]))
    
for line in file_sixth:
    splitLine=line.split(" ")
    Force_sixth.append(float(splitLine[1]))
    
for line in file_seventh:
    splitLine=line.split(" ")
    Force_seventh.append(float(splitLine[1]))
    
for line in file_eightth:
    splitLine=line.split(" ")
    Force_eightth.append(float(splitLine[1]))

for line in file_nighth:
    splitLine=line.split(" ")
    Force_nighth.append(float(splitLine[1]))

for line in file_tenth:
    splitLine=line.split(" ")
    Force_tenth.append(float(splitLine[1]))

for line in file_eleventh:
    splitLine=line.split("	")
    Displacement.append(float(splitLine[0]))
    Force_eleventh.append(float(splitLine[1]))

for line in file_eleventh:
    splitLine=line.split("	")
    Displacement.append(float(splitLine[0]))
    Force_eleventh.append(float(splitLine[1]))

for line in file_twelevth:
    splitLine=line.split("	")
    Displacement_unload.append(float(splitLine[0]))
    Force_twelevth.append(float(splitLine[1]))

for line in file_thirteenth:
    splitLine=line.split("	")
    Displacement_reload.append(float(splitLine[0]))
    Force_thirteenth.append(float(splitLine[1]))

plt.subplot(1,2,1)

plt.plot(Displacement, Force_eleventh, "k", linewidth=3, linestyle='-')
plt.plot(Displacement_unload, Force_twelevth, "b", linewidth=3, linestyle='-')
plt.plot(Displacement_reload, Force_thirteenth, "r", linewidth=2, linestyle='-.')

plt.xticks(color='k', size=15)
plt.yticks(color='k', size=15)
plt.legend(('Load', 'Unload', "Reload") ,shadow=False, loc=(0.05, 0.85), handlelength=1.5, fontsize=16)
plt.xlabel("Displacement ($\AA$)",{'color': 'k', 'fontsize': 22})
plt.ylabel("Force$(eV/\AA)$",{'color': 'k', 'fontsize': 22})

plt.subplot(1,2,2)
plt.plot(Displacement, Force_tenth, "k", linewidth=3, linestyle=":")
plt.plot(Displacement, Force_tenth, "b", linewidth=2, linestyle=":")
plt.plot(Displacement, Force_tenth, "r", linewidth=1, linestyle="-.")

plt.xticks(color='k', size=15)
plt.yticks(color='k', size=15)

plt.xlabel("Displacement ($\AA$)",{'color': 'k', 'fontsize': 22})
plt.ylabel("Force$(eV/\AA)$",{'color': 'k', 'fontsize': 22})


#plt.legend(('$(1) K_F=3.544$', '$(2) K_F=3.708$', "$(3) K_F=3.769$", "$(4) K_F=3.568$", "$(5) K_F=1.992$", "$(6) K_F=1.259$", "$(7) K_F=0.424$", "(8) Energy Minimization"  )\
#,shadow=False, loc=(0.05, 0.65), handlelength=1.5, fontsize=16)

plt.legend(('Load', 'Unload', "Reload") ,shadow=False, loc=(0.05, 0.85), handlelength=1.5, fontsize=16)

plt.show()
    

file.close()

#plt.plot(Displacement, Force_first, "k", linewidth=1, linestyle='-')
#plt.plot(Displacement, Force_second, "b", linewidth=1, linestyle='-')
#plt.plot(Displacement, Force_seventh, "m", linewidth=4, linestyle='-')
#plt.plot(Displacement, Force_eightth, "d", linewidth=5, linestyle='-')
#plt.plot(Displacement, Force_nighth, "b", linewidth=6, linestyle='-')
#plt.plot(Displacement, Force_tenth, "y", linewidth=7, linestyle=':')


#plt.plot(Displacement, Force_fourth, "y", linewidth=2, linestyle='-')
#plt.plot(Displacement, Force_sixth, "r", linewidth=3, linestyle='-')
