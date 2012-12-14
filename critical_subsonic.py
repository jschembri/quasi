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

dict = {"Subsonic": './subsonic'}
for key,value in dict.items():

	data= subprocess.Popen('%s %s' % (value, calculated_time), shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate()
	data = data[0].split(",")
	x_values = data[data.index('X Value Start')+1:data.index("X Value End")]
	for i in range(0,len(x_values)):
		x_values[i] = float(x_values[i])
	y_values = data[data.index('Y Value Start')+1:data.index("Y Value End")]
	for i in range(0,len(y_values)):
		y_values[i] = float(y_values[i])

	y_area_values = data[data.index('Area Start')+1:data.index("Area End")]
	for i in range(0,len(y_area_values)):
		y_area_values[i] = float(y_area_values[i])

	y_velocity_values = data[data.index('Velocity Start')+1:data.index("Velocity End")]
	for i in range(0,len(y_velocity_values)):
		y_velocity_values[i] = float(y_velocity_values[i])

	y_pressure_values = data[data.index('Pressure Start')+1:data.index("Pressure End")]
	for i in range(0,len(y_pressure_values)):
		y_pressure_values[i] = float(y_pressure_values[i])

	y_mach_values = data[data.index('Mach Start')+1:data.index("Mach End")]
	for i in range(0,len(y_mach_values)):
		y_mach_values[i] = float(y_mach_values[i])

	y_temperature_values = data[data.index('Temperature Start')+1:data.index("Temperature End")]
	y_temperature_values = data[data.index('Temperature Start')+1:data.index("Temperature End")]

	plt.subplot(311)
	plt.plot(x_values, y_mach_values, ms=12, label='Numeric',linewidth=4)
	plt.ylabel('Mach')

	plt.subplot(312)
	plt.plot(x_values, y_area_values, ms=12, label='area',linewidth=4)
	plt.ylabel('Area')

	plt.subplot(313)
	plt.plot(x_values, y_values, ms=12, label='Numeric',linewidth=4)
	plt.ylabel('Density')
	plt.xlabel('X-Axes')


#data= subprocess.Popen('./analytic_subsonic', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate()
#data = data[0].split(",")
#x_analytic_values = data[data.index('X Value Start')+1:data.index("X Value End")]
#for i in range(0,len(x_analytic_values)):
#	x_analytic_values[i] = float(x_analytic_values[i])
#y_analytic_values = data[data.index('Mach Start')+1:data.index("Mach End")]
#for i in range(0,len(y_analytic_values)):
#	y_analytic_values[i] = float(y_analytic_values[i])
#y_row_values = data[data.index('Row Start')+1:data.index("Row End")]
#for i in range(0,len(y_row_values)):
#	y_row_values[i] = float(y_row_values[i])
#y_pressure_values = data[data.index('Pressure Start')+1:data.index("Pressure End")]
#for i in range(0,len(y_pressure_values)):
#	y_pressure_values[i] = float(y_pressure_values[i])


#y_temperature_values = data[data.index('Temperature Start')+1:data.index("Temperature End")]

#plt.subplot(311)
#plt.plot(x_analytic_values, y_analytic_values, ms=12, label="Analytic",linewidth=4)
#plt.ylabel('Mach Number')
#plt.legend()

#plt.subplot(313)
#plt.plot(x_analytic_values, y_row_values, ms=12, label='Analytic',linewidth=4)
#plt.ylabel('Density')
#plt.legend()
#plt.xlabel('X-Axes')

#plt.subplot(513)
#plt.plot(x_analytic_values, y_temperature_values, ms=12, label='analytic',linewidth=4)
#plt.ylabel('Temperature')

#plt.subplot(514)
#plt.plot(x_analytic_values, y_pressure_values, ms=12, label='analytic',linewidth=4)
#plt.ylabel('Pressure')




#plt.axis([0,4,0,1.5])

show()

