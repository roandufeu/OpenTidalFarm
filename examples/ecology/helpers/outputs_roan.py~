#!/usr/bin/python
import sys
import re
import matplotlib.pyplot as plt
import argparse
from matplotlib import rc

parser = argparse.ArgumentParser(description='Create a convergence plot from \
                                 OpenTidalFarm stdout file.')
parser.add_argument('file', type=str, help='the filename containing the \
                    OpenTidalFarm stdout')

args = parser.parse_args()

sed_scaling_factor = 1 #5
eco_scaling_factor = 1 #0.4

f = open(args.file, "r")
# f = open(args.file, "output.txt")
# Read the output file of the optimisation run and parse it
# This script assumes that scipy.slsqp was used.

found_slsqp = False
# The iteration numbers
it = []
last_it = 0
# The functional values
func_pow = []
func_eco = []
func_sed = []
# The number of functional evaluations
func_evals = []
finished = False
# Rescale the  functional values before plotting
# Note: this does not yet include the automatic scaling performed by the turbine
# code. We will recover that factor from the output file
rescale = 10**-6
rescale_sed = rescale / sed_scaling_factor
rescale_eco = rescale / eco_scaling_factor

it.append(int(1))

for line in f:
    if "automatic_scaling:" in line:
        m = re.match(r".*: ([0-9|\.]+)", line)
        print "Found rescaling factor: ", m.group(1)
        rescale_pow = rescale/float(m.group(1))
    if "  NIT    FC           OBJFUN            GNORM" in line:
        found_slsqp = True

    if re.match(r".* ([0-9]+) \s+ ([0-9]+) \s+ ([0-9|\.|E\+\-]+) \s+ ([0-9|\.|E\+\-]+)", line):
        m = re.match(r".* ([0-9]+) \s+ ([0-9]+) \s+ ([0-9|\.|E\+\-]+) \s+ ([0-9|\.|E\+\-]+)", line)
        it.append(int(m.group(1))+1)
        func_evals.append(int(m.group(2)))

    if "power=" in line:
        m = re.match(r".* ([0-9|\.]+)", line)
        if m.group(1) != '0.0':
            func_pow.append((it[-1], float(m.group(1))*rescale_pow))

    if "ecology =" in line:
        m = re.match(r".* ([0-9|\.]+)", line)
        if m.group(1) != '0.0':
            func_eco.append((it[-1], float(m.group(1))*rescale_eco))

    if "sediment =" in line:
        m = re.match(r".* ([0-9|\.]+)", line)
        if m.group(1) != '0.0':
            func_sed.append((it[-1], float(m.group(1))*rescale_sed))

f.close()

if not found_slsqp:
    print "No SLSQP output found. Please supply the stdout record of an \
    OpenTidalFarm simulation which used the SLSQP optimisation algorithm."
    sys.exit(1)


repeats_pow = []
repeats_eco = []
repeats_sed = []
power = []
ecology = []
sediment = []

# from IPython import embed; embed()

if len(func_pow) > 0:
    for i in range(1, it[-2]):
        for j in range(0, len(func_pow)):
            if func_pow[j][0] == i:
                repeats_pow.append(func_pow[j][1])
        if len(repeats_pow) > 0:
            power.append(repeats_pow[-1])
        repeats_pow = []
    print "Power output of initial layout: ", power[0]
    print "Power output of final layout: ", power[-1]
    print "Relative power increase: ", power[-1]/power[0]
    print "Number of iterations: ", len(power)
    # Produce plot
    scaling = 0.7
    rc('text', usetex=True)
    plt.figure(1, figsize=(scaling*7., scaling*4.))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.plot(power, color='black')
    plt.ylabel(r"Power production [MW]")
    plt.xlabel(r"Optimisation iteration")
    plt.savefig("power.pdf")
    plt.close()

if len(func_eco) > 0:
    for i in range(1, it[-2]):
        for j in range(0, len(func_eco)-1):
            if func_eco[j][0] == i:
                repeats_eco.append(func_eco[j][1])
        if len(repeats_eco) > 0:
            ecology.append(repeats_eco[-2])
        repeats_eco = []
# from IPython import embed; embed()
    print "Initial ecology: ", ecology[0]
    print "Final ecology: ", ecology[-1]
    print "Relative change in ecology: ", ecology[-1]/ecology[0]
    print "Number of iterations: ", len(ecology)
    # Produce plot
    scaling = 0.7
    rc('text', usetex=True)
    plt.figure(1, figsize=(scaling*7., scaling*4.))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.plot(ecology, color='black')
    plt.ylabel(r"Ecology")
    plt.xlabel(r"Optimisation iteration")
    plt.savefig("ecology.pdf")
    plt.close()

if len(func_sed) > 0:
    for i in range(1, it[-2]):
        for j in range(0, len(func_sed)):
            if func_sed[j][0] == i:
                repeats_sed.append(func_sed[j][1])
        if len(repeats_sed) > 0:
            sediment.append(repeats_sed[-2])
        repeats_sed = []
    print "Initial sediment: ", sediment[0]
    print "Final sediment: ", sediment[-1]
    print "Relative change in sediment: ", sediment[-1]/sediment[0]
    print "Number of iterations: ", len(sediment)
    # Produce plot
    scaling = 0.7
    rc('text', usetex=True)
    plt.figure(1, figsize=(scaling*7., scaling*4.))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.plot(sediment, color='black')
    plt.ylabel(r"Sediment")
    plt.xlabel(r"Optimisation iteration")
    plt.savefig("sediment.pdf")
    plt.close()
