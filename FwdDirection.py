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

dict = {"FiniteVolume": './finiteVolume', "FiniteDifference": './finiteDifference'}
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




	plt.subplot(411)
	plt.plot(x_values, y_area_values, ms=12, label='%s' % key,linewidth=4)
	plt.ylabel('Area')


	plt.subplot(412)
	plt.plot(x_values, y_values, ms=12, label='%s' % key,linewidth=4)
	plt.ylabel('Density')

	plt.subplot(413)
	plt.plot(x_values, y_velocity_values, ms=12, label='%s' % key,linewidth=4)
	plt.ylabel('Velocity')

	plt.subplot(414)
	plt.plot(x_values, y_pressure_values, ms=12, label='%s' % key,linewidth=4)
	plt.ylabel('Pressure')


plt.legend()
plt.xlabel('X-Axes')



#plt.axis([0,4,0,1.5])

show()

