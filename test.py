#! /usr/bin/env python
import os
import subprocess
import shutil
import itertools

error = [ ]
samplesize = []
highreserror = []
reserror = []
work = []
compare = []
compareres =[]



# for k in range ( 1 ) :
#     Nx = 2**k + 1
#     options = ['-N ', str(Nx)]
#     error.append(os.system('./boundaryvalues '+' '.join(options)))
#     samplesize.append(Nx)
# print zip(error, samplesize)

#timestepping took 53.233215 seconds clocktime 
#timestepping took 7 seconds walltime 
#updating took 28.403000 seconds clocktime 


#for k in range (1,8):
for (radius,solver_opt) in itertools.product(range(1,8+1), ["pre", "opt", "slow"]):
	cmd = './naivefd --radius ' + str(radius) + " --solver " + str(solver_opt)
	print cmd
	out = subprocess.check_output(cmd, shell=True) 
	splitstring = out.split("\n")
	for linenumber, lines in enumerate(splitstring):
		start = (lines.find("timestepping took"))
		if start > -1:
				print(float(lines.split()[2]))
				#error.append(float(lines[start + len("error vector "):]))
				#compare.append(1.0/(Nx)**2)
				#samplesize.append(Nx)
				print(lines, linenumber)
	#shutil.copy("soln_nx=100_ny=100_nz=100_nt=50.png", "test/heatmap" +cmd + ".png")

'''
from pylab import legend, plot, loglog, show, title, xlabel, ylabel, savefig
#plot(samplesize, error, '--k')
loglog(samplesize, error, '--k', samplesize, compare, 'r--')
legend(['cheby accuracy', '1/h^2'])
title('u(x) = (2x^2-2)*(2y^2-2) on patch [-1, 0]x[-1, 0]')
xlabel('Problem Size $N$')
ylabel('Error (N)')
show()

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
