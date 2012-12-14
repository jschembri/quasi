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


	x_iter_values = data[data.index('R1 Start')+1:data.index("R1 End")]
	x_R2_values = data[data.index('R2 Start')+1:data.index("R2 End")]
	x_R3_values = data[data.index('R3 Start')+1:data.index("R3 End")]
	y_res_values = data[data.index('iteration_list Start')+1:data.index("iteration_list End")]

	y_error = data[data.index('residual_list Start')+1:data.index("residual_list End")]

#	plt.subplot(211)
#	plt.plot(y_res_values, y_error, ms=12, label='R1',linewidth=4)
#	plt.ylabel('Error')
#	plt.yscale('log')


#	plt.subplot(212)
	plt.plot(y_res_values, x_iter_values, ms=12, label='R1',linewidth=4)
	plt.plot(y_res_values, x_R2_values, ms=12, label='R2',linewidth=4)
	plt.plot(y_res_values, x_R3_values, ms=12, label='R3',linewidth=4)
	plt.ylabel('Residual')
	plt.yscale('log')
	plt.legend()
	plt.xlabel('X-Axes')



show()

