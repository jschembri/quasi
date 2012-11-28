# Quasi 1D code
# created by Jeremy Schembr on November 8, 2012

from sys import argv
from matplotlib.patches import Patch
from pylab import *
import subprocess


if len(argv) < 2:
   calculated_time= "There is none"
else:
   script, calculated_time = argv

data= subprocess.Popen('./pressure', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate()
data = data[0].split(",")
x_values = data[data.index('X Value Start')+1:data.index("X Value End")]
for i in range(0,len(x_values)):
   x_values[i] = float(x_values[i])

y_area_values = data[data.index('Area Start')+1:data.index("Area End")]
for i in range(0,len(y_area_values)):
   y_area_values[i] = float(y_area_values[i])



y_pressure_values = data[data.index('Pressure Start')+1:data.index("Pressure End")]
for i in range(0,len(y_pressure_values)):
   y_pressure_values[i] = float(y_pressure_values[i])



plt.subplot(211)
plt.plot(x_values, y_pressure_values, mfc='red', ms=12, label='Velocity',linewidth=4)
plt.ylabel('Pressure')

plt.subplot(212)
plt.plot(x_values, y_area_values, mfc='red', ms=12, label='Area',linewidth=4)
plt.ylabel('Area')




plt.xlabel('X-Axes')



#plt.axis([0,4,0,1.5])

show()

