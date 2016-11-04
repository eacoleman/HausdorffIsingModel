#!/usr/bin/env python

import numpy
from numpy import *
from pylab import *  # plotting library
from scipy import weave
from scipy.weave import converters

"""
Edits made by Evan A. Coleman
Summary:

    - Implemented magnetic field strength element to class
    - Changed calculation of Boltzmann elements on the fly
    - Got rid of unused methods

The actual functionality modifications I make amount to ~5 lines

"""

class Ising2D (object):

    """Class that describes equilibrium statistical mechanics of the two-dimensional Ising model"""

    def __init__(self, L=32, temperature=10.0, h=0.0):

        self.L = L    # length of the grid
        self.N = L**2 # number of particles
        self.h = h    # magnetic interaction strength

        self.temperature = temperature

        self.w = zeros(9) # store Boltzmann weights
        self.w[8] = exp(-8.0/self.temperature)
        self.w[4] = exp(-4.0/self.temperature)

        self.state = ones((self.L, self.L), int) # initially all spins up
        self.energy = - 2 * self.N
        self.magnetization = self.N

        self.reset()


    def reset(self):

        self.monteCarloSteps = 0
        self.acceptedMoves = 0
        self.energyArray = array([], int)
        self.magnetizationArray = array([], int)

    def rasterMonteCarloStep(self):

        N = self.N
        L = self.L
        w = self.w
        h = self.h

        # EDIT:
        #  - calculate the magnetic interaction term
        magInteraction = exp(h/self.temperature)

        # setup helpers
        state = self.state
        acceptedMoves = array([self.acceptedMoves], int)
        energy = array([self.energy], int)
        magnetization = array([self.magnetization], int)

        randomArray = numpy.random.random(N)

        # run code to flip spins where energy / randomness favors
        # EDIT:
        #  - if statement now takes into account difference of energy in Boltzmann terms
        code = """
        int L1 = L-1;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {

                int dE = 2*state(i, j) * (state((i+1)%L, j) + state((i+L1)%L, j) + state(i, (j+1)%L) + state(i, (j+L1)%L));

                if ( (dE+h*state(i, j)) <= 0
                     || w(abs(dE))*pow(magInteraction,-state(i, j)) > randomArray(i*L+j)) {
                    acceptedMoves(0) += 1;
                    int newSpin = -state(i, j);
                    state(i, j) = newSpin;
                    energy(0) += dE;
                    magnetization(0) += 2*newSpin;
                }
            }
        }
        """

        # EDIT:
        #  - add necessary items to compiler environment
        weave.inline(code, ['state', 'L', 'w', 'h', 'magInteraction', 'randomArray', 'acceptedMoves', 'energy', 'magnetization'],
                type_converters=converters.blitz, compiler='gcc')


        self.state = state
        self.acceptedMoves = acceptedMoves[0]
        self.energy = energy[0]
        self.magnetization = magnetization[0]

        self.energyArray = append(self.energyArray, self.energy)
        self.magnetizationArray = append(self.magnetizationArray, self.magnetization)
        self.monteCarloSteps += 1


    def steps(self, number = 100):
        for i in range(number):
            self.rasterMonteCarloStep()

    # Observables
    def meanEnergy(self):
        return self.energyArray.mean() / self.N

    def specificHeat(self):
        return (self.energyArray.std() / self.temperature)**2 / self.N

    def meanMagnetization(self):
        return self.magnetizationArray.mean() / self.N

    def susceptibility(self):
        return (self.magnetizationArray.std())**2 / (self.temperature * self.N)

    def observables(self):
        print "\nTemperature = ", self.temperature
        print "Mean Energy = ", self.meanEnergy()
        print "Mean Magnetization = ", self.meanMagnetization()
        print "Specific Heat = ", self.specificHeat()
        print "Susceptibility = ", self.susceptibility()
        print "Magnetic Field Strength = ", self.h
        print "Monte Carlo Steps = ", self.monteCarloSteps, " Accepted Moves = ", self.acceptedMoves


    # Visual snapshot of state
    def plot(self):

        x = indices((self.L, self.L))[0,:,:]
        y = indices((self.L, self.L))[1,:,:]
        color = []
        for i in range(self.L):
            for j in range(self.L):
                if self.state[i, j] == 1:
                    color.append([1.0, 0.0, 0.0])
                else:
                    color.append([0.0, 0.0, 1.0])
        scatter(x.ravel(), y.ravel(), s=30.0, c=color, edgecolors='none')



model = Ising2D(temperature=0.1, L=32,h=1)
figure(1)
model.plot()

model.steps(number=10000)
model.reset()

# do a second pass
model.steps(number=10000)
model.observables()
figure(2)
model.plot()

# do a third pass
model.steps(number=10000)
model.observables()
figure(3)
model.plot()


show()
