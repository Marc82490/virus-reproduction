import random
import pylab
import numpy as np

random.seed(0)                                         # Debugging statement
''' 
Begin helper code
'''

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. You can use NoChildException as is, you do not need to
    modify/add any code.
    """

'''
End helper code
'''

#
# PROBLEM 1
#
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb

    def getMaxBirthProb(self):
        """
        Returns the max birth probability.
        """
        return self.maxBirthProb

    def getClearProb(self):
        """
        Returns the clear probability.
        """
        return self.clearProb

    def doesClear(self):
        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """
        if self.getClearProb() > random.random():
            return True
        else:
            return False
    
    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """
        if self.maxBirthProb * (1 - popDensity) > random.random():
            return SimpleVirus(self.maxBirthProb, self.clearProb)
        else:
            raise NoChildException


class Patient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """    

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the maximum virus population for this patient (an integer)
        """
        self.viruses = viruses
        self.maxPop = maxPop

    def getViruses(self):
        """
        Returns the viruses in this Patient.
        """
        return self.viruses

    def getMaxPop(self):
        """
        Returns the max population.
        """
        return self.maxPop

    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
        returns: The total virus population (an integer)
        """
        return len(self.viruses)

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.   
        
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        # Check each virus to see if it is removed in this time step.
        for virus in self.viruses[:]:
            if virus.doesClear():
                self.viruses.remove(virus)
        
        # Calculate the current population density.        
        popDensity = self.getTotalPop() / self.getMaxPop()
        
        # If the max population density has not been reached, check each remaining
        #  virus to see if it reproduces in this time step.
        if popDensity <= 1:
            for virus in self.viruses[:]:
                try:
                    self.viruses.append(virus.reproduce(popDensity))
                except NoChildException:
                    pass
        
        return self.getTotalPop()


#
# PROBLEM 2
#
def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph for problem 3 (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    """
    # Initialize array of zeroes for each time step.
    average_population_size = np.zeros( (1,300) )
    # Run each trial.
    for trial in range(numTrials):
        # Initialize the viruses.
        viruses = []
        for virus in range(numViruses):
            viruses.append(SimpleVirus(maxBirthProb, clearProb))

        # Initialize the patient.    
        patient = Patient(viruses, maxPop)
        
        # Initialize a list to record the virus population sizes over time in this trial.
        population_size = []
        
        # Update the patient and record the virus population size at each time step.
        for timestep in range(300):
            patient.update()
            population_size.append(patient.getTotalPop())
        
        # Convert the population list into a numpy array.
        trial_population = np.array(population_size)
        # Add this trial's population sizes to the sum at each time step.
        average_population_size += trial_population
    
    # Divide the population size at each time step by the number of trials
    #  in order to get the average population size at each time step.
    average_population_size /= float(numTrials)
    # Unnests the the array from the outer list, so that it can be plotted with pylab.
    raveled = np.ravel(average_population_size)
    
    pylab.plot(raveled.tolist(), label = 'Simple Virus')
    pylab.title('Average Virus Population for ' + str(numTrials) + ' Trials')
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population Size')
    pylab.legend()
    pylab.show()

#~ print(simulationWithoutDrug(100, 1000, 0.1, 0.05, 100))
#
# PROBLEM 3
#
class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """   

    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)       

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """
        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        self.resistances = resistances
        self.mutProb = mutProb

    def getResistances(self):
        """
        Returns the resistances for this virus.
        """
        return self.resistances

    def getMutProb(self):
        """
        Returns the mutation probability for this virus.
        """
        return self.mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        return self.getResistances().get(drug, False)

    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.maxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """
        # Check the virus's resistance to each active drug.
        drugs_copy = activeDrugs[:]
        for drug in drugs_copy:
            # If the virus is resistant to that drug, remove the drug from the list.
            if self.isResistantTo(drug):
                drugs_copy.remove(drug)

        # If there are no more active drugs, the virus can reproduce.
        if not drugs_copy:
            if self.maxBirthProb * (1 - popDensity) > random.random():
                # Initialize the offspring's resistance dictionary.
                offspring_resistances = {}
                # For each resistance the parent has, determine whether the
                #  offspring will inherit the same status or mutate. Adds the
                #  resulting status to the offspring's resistance dictionary.
                for drug in self.getResistances():
                    if self.isResistantTo(drug):
                        if (1 - self.mutProb) > random.random():
                            offspring_resistances[drug] = True
                        else:
                            offspring_resistances[drug] = False
                    else:
                        if self.mutProb > random.random():
                            offspring_resistances[drug] = True
                        else:
                            offspring_resistances[drug] = False
                return ResistantVirus(self.maxBirthProb, self.clearProb, offspring_resistances, self.mutProb)
            else:
                raise NoChildException
        
        raise NoChildException
            

class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """
        Patient.__init__(self, viruses, maxPop)
        self.prescriptions = []

    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """
        if newDrug not in self.prescriptions:
            self.prescriptions.append(newDrug)

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """
        return self.prescriptions

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        resistant_viruses = 0
        for virus in self.getViruses():
            resistant = True
            for drug in drugResist:
                if not virus.isResistantTo(drug):
                    resistant = False
            if resistant:
                resistant_viruses += 1               
        
        return resistant_viruses

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """
        # Check each virus to see if it is removed in this time step.
        for virus in self.viruses[:]:
            if virus.doesClear():
                self.viruses.remove(virus)
        
        # Calculate the current population density.        
        popDensity = self.getTotalPop() / self.getMaxPop()
        
        # If the max population density has not been reached, check each remaining
        #  virus to see if it reproduces in this time step.
        if popDensity <= 1:
            for virus in self.viruses[:]:
                try:
                    self.viruses.append(virus.reproduce(popDensity, self.getPrescriptions()))
                except NoChildException:
                    pass
        
        return self.getTotalPop()


#
# PROBLEM 4
#
def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials):
    """
    Runs simulations and plots graphs for problem 5.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
    """
    # Initialize two arrays of zeroes for each time step.
    avg_total_pop_size = np.zeros( (1,300) )
    avg_resist_pop_size = np.zeros( (1,300) )
    
    # Run each trial.
    for trial in range(numTrials):
        # Initialize the viruses.
        viruses = []
        for virus in range(numViruses):
            viruses.append(ResistantVirus(maxBirthProb, clearProb, resistances, mutProb))

        # Initialize the patient.    
        patient = TreatedPatient(viruses, maxPop)
        
        # Initialize a list to record the total virus population sizes over time.
        total_population_size = []
        
        # Initialize a list to record the resistant virus population sizes over time.
        resistant_population_size = []
        
        # Update the patient and record the virus population sizes at each time step.
        for timestep in range(300):
            # On the 150th time step, add a prescription for 'guttagonol'.
            if timestep == 150:
                patient.addPrescription('guttagonol')
            patient.update()
            total_population_size.append(patient.getTotalPop())
            resistant_population_size.append(patient.getResistPop(['guttagonol']))
        
        # Convert the trial's population lists into numpy arrays.
        trial_total_pop = np.array(total_population_size)
        trial_resist_pop = np.array(resistant_population_size)
        # Add this trial's population sizes to the sum at each time step.
        avg_total_pop_size += trial_total_pop
        avg_resist_pop_size += trial_resist_pop
    
    # Divide the population sizes at each time step by the number of trials
    #  in order to get the average population sizes at each time step.
    avg_total_pop_size /= float(numTrials)
    avg_resist_pop_size /= float(numTrials)
    # Unnests the the arrays from the outer list, so that they can be plotted with pylab.
    raveled_total_pop = np.ravel(avg_total_pop_size)
    raveled_resist_pop = np.ravel(avg_resist_pop_size)
    
    pylab.plot(raveled_total_pop.tolist(), label = 'Total Population')
    pylab.plot(raveled_resist_pop.tolist(), label = 'Resistant Population')
    pylab.title('Average Virus Population for ' + str(numTrials) + ' Trials')
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population Size')
    pylab.legend()
    pylab.show()
    
#~ simulationWithDrug(100, 1000, 0.1, 0.05, {'guttagonol': False}, 0.005, 100)
