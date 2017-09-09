#! /usr/bin/env python
import os
import subprocess
import shutil
import itertools
import numpy as np

update1 = [ ]
update2 = [ ]

error1 = [ ]
error2 = [ ]
error3 = [ ]
radius = []


# for k in range ( 1 ) :
#     Nx = 2**k + 1
#     options = ['-N ', str(Nx)]
#     error.append(os.system('./boundaryvalues '+' '.join(options)))
#     samplesize.append(Nx)
# print zip(error, samplesize)

#timestepping took 53.233215 seconds clocktime 
#timestepping took 7 seconds walltime 
#updating took 28.403000 seconds clocktime 


radiussample=range(1,8+1)

for (radius,solver_opt) in itertools.product(radiussample, ["pre", "opt", "slow"]):
	cmd = './naivefd --radius ' + str(radius) + " --solver " + str(solver_opt)
	print cmd
	out = subprocess.check_output(cmd, shell=True) 
	#print out
	splitstring = out.split("\n")
	for linenumber, lines in enumerate(splitstring):
		start = (lines.find("timestepping took"))
		if start > -1:
			if solver_opt=="slow":
				print(float(lines.split()[2]))
				error1.append(float(lines.split()[2]))
			if solver_opt=="pre":
				print(float(lines.split()[2]))
				error2.append(float(lines.split()[2]))
			if solver_opt=="opt":
				print(float(lines.split()[2]))
				error3.append(float(lines.split()[2]))
			

			print(lines)
		start2 = (lines.find("updating took"))
		if start2 > -1:
			if solver_opt=="pre":
				print(float(lines.split()[2]))
				update1.append(float(lines.split()[2]))
			if solver_opt=="opt":
				print(float(lines.split()[2]))
				update2.append(float(lines.split()[2]))
			print(lines)
	#shutil.copy("soln_nx=100_ny=100_nz=100_nt=50.png", "test/heatmap" +cmd + ".png")

from pylab import legend, plot, loglog, show, title, xlabel, ylabel, savefig
#plot(samplesize, error, '--k')
loglog(radiussample, error1, '--k', radiussample, error2, 'r--', radiussample, error3, 'r')
legend(['slow', 'precompute', 'opt'])
title('Scaling of various implementations')
xlabel('Radius N')
ylabel('Compute time')
show()


#plot(samplesize, error, '--k')
loglog(radiussample, update1, '--k', radiussample, update2, 'r--')
legend(['precompute', 'opt'])
title('Time spent writing to memory with new method')
xlabel('Radius' )
ylabel('Compute time')
show()


from matplotlib import pyplot

width=.35 
fig, ax = pyplot.subplots()
base1 = [error_i- update_i for error_i, update_i in zip(error2, update1)]
rects1 = ax.bar(radiussample, base1, width, color='maroon')
rects1plus = ax.bar(radiussample, update1, width, color='lightcoral',
             bottom=base1)

shifted = [radius + width for radius in radiussample]
base2 = [error_i- update_i for error_i, update_i in zip(error3, update2)]
rects2 = ax.bar(shifted, base2, width, color='navy')
rects2plus = ax.bar(shifted, update2, width, color='skyblue',
             bottom=base2)
ax.set_ylabel('Time for each step')
ax.set_title('Makeup of runtimes', y=1.05)
#ax.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
#ax.yticks(np.arange(0, 81, 10))
#ax.legend((rects1[0], rects2[0], rects1plus[0], rects2plus[0]), ('Precompute-update', 'Opt-update', 'Precompute-total', 'Opt-total'))
ax.legend((rects1[0], rects2[0], rects1plus[0], rects2plus[0]), ('Compute time:Pre-compute', 'Compute time:New Method', 'Update time: Pre-compute', 'Update time:New Method'), ncol=2, fancybox=True, shadow=True, loc='upper center', bbox_to_anchor=(0.5, 1.05))
pyplot.show()

'''

#plot(work, reserror, '--k')
loglog(work, reserror, '--k', samplesize, compareres, 'r--')
legend(['cheby precision vs work', '1/h^2'])
title('u(x) = (2x^2-2)*(2y^2-2) on patch [-1, 0]x[-1, 0]')
xlabel('Problem Size $N$')
ylabel('Error (N)')
show()

# plot(samplesize, highreserror, '--k')
# legend(['fft accuracy'])
# title('fft')
# xlabel('Problem Size $N$')
# ylabel('High Resolution Error (N)')
# show()
#loglog(sizes, times)
#title('SNES ex5')
#xlabel('Problem Size $N$')
#ylabel('Time (s)')
#show()

'''
