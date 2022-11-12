"""
This code calculates and plots the infromation transffered by a single
quantum state measured in the correct basis in our protocol.

This calculation is done under the two assumptions:
1. Low Amplification - we assume that the squeezing in the non-linear crystals
    was low (small alpha) such that the quantum state will contain one
    (signal, idler) pair at most.
    Here we will show also the results for large alpha but these results
    should be viewed carefully.
2. Randomaly Chosen Phase - we assume that Alice created the state with a
    random phase between 0, pi/2, pi and 3pi/2. Specifically, we assume that
    in the correct basis (Without loss of generality, [0, pi]) the
    probabilities are:
    >> P(phi = 0) = P(phi = pi) = 0.5 
    These probabilities are the standard probabilities for most QKD protocols.
"""

import numpy as np
import matplotlib.pyplot as plt

import os


EPS = 1e-5

# Alice's transmission probabilities
P_ZERO = 0.5
P_PI = 0.5

def phase_measurement_probs(loss, alpha, zero_phase):
    """
    Calculates the probabilities to measure different cases for a given phase.

            Parameters:
                    loss (float): A float between 0 and 1 represents the loss
                        in the setup (the loss in the intensity).
                    alpha (float): A float between 0 and 1 represents the 
                        field squeezing.
                    zero_phase (boolean) : A boolean represents Alice's phase.
                        zero_phase = True -> phase = 0
                        zero_phase = False -> phase = pi
            Returns:
                    probs (np.array): A numpy array with the probabilities to
                        measure the 4 different cases: 
                        (signal, idler) = (1, 1), (1, 0), (0, 1) or (0, 0)
    """
    # Calculating the transmission coefficient
    trans = 1 - loss
    # Calculating the unnormalized probabilities
    p_1_1 = (trans + 1) ** 2 if zero_phase else (trans - 1) ** 2
    p_1_0 = loss * trans
    p_0_1 = loss * trans
    p_0_0 = (loss - (1/alpha)) ** 2
    # Normalizing the probabilities
    sum_p = p_1_1 + p_1_0 + p_0_1 + p_0_0
    p_1_1 /= sum_p
    p_1_0 /= sum_p
    p_0_1 /= sum_p
    p_0_0 /= sum_p
    # Creating the final array
    probs = np.array([p_1_1, p_1_0, p_0_1, p_0_0])
    return probs

def single_state_information(loss, alpha):
    """
    Calculates the infromation transffered by a single quantum state.

            Parameters:
                    loss (float): A float between 0 and 1 represents the loss
                        in the setup (the loss in the intensity).
                    alpha (float): A float between 0 and 1 represents the 
                        field squeezing.
            Returns:
                    information (float): The calculated information (in bits).
    """
    # The probabilities to measure each case in each transmission 
    z_measure_probs = P_ZERO*phase_measurement_probs(loss, alpha, zero_phase=True)
    pi_measure_probs = P_PI*phase_measurement_probs(loss, alpha, zero_phase=False)

    # The total measurement probabilities
    measure_probs = z_measure_probs + pi_measure_probs

    # The mutual information between the input and output of the channel
    z_info = z_measure_probs * np.log2(z_measure_probs / (P_ZERO * measure_probs))
    pi_info = pi_measure_probs * np.log2(pi_measure_probs / (P_PI * measure_probs))

    information = np.sum(z_info + pi_info)
    
    return information

def calculate_simulation(data_path, loss_samples=25, alpha_samples=25):
    """
    Calculates the single states information for multiple losses and gains.
    The results are saved in the data dir for future uses.

            Parameters:
                    data_path (str): The path to which the simulation data
                        will be saved.
                    loss_samples (int): The number of losses to calculate.
                    alpha_samples (int): The number of gains to calculate.
    """
    # Initializing the calculation space.
    loss_array = np.linspace(EPS, 1-EPS, loss_samples)
    alpha_array = np.linspace(EPS, 1, alpha_samples)

    # Initializing the results array.
    info_array = np.zeros((alpha_samples, loss_samples))

    # Calculating the information over the calculation space.
    for loss_ind, loss in enumerate(loss_array):
        for alpha_ind, alpha in enumerate(alpha_array):
            info_array[alpha_ind, loss_ind] = single_state_information(loss, alpha)

    # Saves the calculations.
    np.save(data_path, info_array)

def simulate(data_path, save_path, reprocess=False):
    """
    Calculates the information and saves the results in the given path.

            Parameters:
                    data_path (str): The path in which the simulation data is
                        saved.
                    save_path (str): The path to which the final graph will be
                        saved.
                    reprocess (boolean): A boolean telling whether to
                        recalculate the information if it was already
                        calculated before. Default as False.
    """
    # Testing whether to run the simulation calculation.
    if not os.path.exists(data_path) or reprocess:
        calculate_simulation(data_path)

    # Loading the simulation data.
    info_array = np.load(data_path)

    # Plotting the simulation results.
    plt.figure()
    plt.imshow(info_array[::-1, :], extent=[0, 1, 0, 1], aspect="auto", 
               cmap='jet', interpolation='bilinear')
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=16)
    plt.xlabel('Loss', fontsize=24)
    plt.ylabel(r'$\alpha_{\omega}$', fontsize=24)
    plt.title('Single Measurement Information', fontsize=30)
    plt.xticks(fontsize=16, rotation=0)
    plt.yticks(fontsize=16, rotation=0)
    fig = plt.gcf()
    fig.set_size_inches((11, 8.5), forward=False)

    # Saving the figure to the given path.
    plt.savefig(save_path, dpi=500)