"""
This code calculates and plots the infromation transffered in the correct
basis during an integration time in our protocol.

To calculate the integration time information we simulate a large number
of quantum states (usally 100) with very low gain such that the average number
of photons during the integration time will remain constant.

This calculation is done under the three assumptions:
1. Low Amplification - we assume that the squeezing in the non-linear crystals
    was low (small alpha) such that each quantum state will contain one
    (signal, idler) pair at most. Here, in contrast to the single state
    information this assumption is automatically applied for a small number of
    photons since we are using a large number of quantum states.
2. Randomaly Chosen Phase - we assume that Alice created the state with a
    random phase between 0, pi/2, pi and 3pi/2. Specifically, we assume that
    in the correct basis (without loss of generality, [0, pi]) the
    probabilities are:
    >> P(phi = 0) = P(phi = pi) = 0.5 
    These probabilities are the standard probabilities for most QKD protocols.
"""

import numpy as np
import matplotlib.pyplot as plt

import os


EPS = 1e-10

# Alice's transmission probabilities.
P_ZERO = 0.5
P_PI = 1 - P_ZERO


def single_state_measurement_probs(loss, alpha, zero_phase):
    """
    Calculates the measurement probabilities of a single state.

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
    # Calculating the transmission coefficient.
    trans = 1 - loss

    # Calculating the unnormalized probabilities.
    p_1_1 = (trans + 1) ** 2 if zero_phase else (trans - 1) ** 2
    p_1_0 = loss * trans
    p_0_1 = loss * trans
    p_0_0 = (loss - (1/alpha)) ** 2

    # Normalizing the probabilities.
    sum_p = p_1_1 + p_1_0 + p_0_1 + p_0_0
    p_1_1 /= sum_p
    p_1_0 /= sum_p
    p_0_1 /= sum_p
    p_0_0 /= sum_p

    # Creating the final array.
    probs = np.array([p_1_1, p_1_0, p_0_1, p_0_0])
    return probs


def create_probability_table(loss, n_photons, n_states, zero_phase):
    """
    Calculates the measurement probabilities of the entire integration time.

            Parameters:
                    loss (float): A float between 0 and 1 representing the
                        loss in the setup (the loss in the intensity).
                    n_photons (float): A float representing the average number
                        of photons in the integration time.
                    n_states (int): The number of simulated quantum states 
                        during the integration time.
                    zero_phase (boolean) : A boolean represents Alice's phase.
                        zero_phase = True -> phase = 0
                        zero_phase = False -> phase = pi
            Returns:
                    prob_table (np.array): A 2d numpy array with the 
                        probabilities to measure the different number of
                        photons in the (signal, idler).
    """
    # Calculating the squeezing factor.
    alpha = pow(n_photons / (4 * n_states), 0.5)
    
    # Calculating the single state probabilities.
    single_state_probs = single_state_measurement_probs(loss, 
        alpha, 
        zero_phase)

    # Initializing the probabilities table.
    prob_table = np.zeros((n_states, n_states))
    prob_table[0, 0] = P_ZERO if zero_phase else P_PI
    prob_table += EPS

    # Update the probabilities for the different quantum states.
    for i in range(n_states):
        new_prob_table = np.zeros((n_states, n_states))
        new_prob_table[1:,1:] += single_state_probs[0] * prob_table[:-1,:-1]
        new_prob_table[1:,:] += single_state_probs[1] * prob_table[:-1,:]
        new_prob_table[:,1:] += single_state_probs[2] * prob_table[:,:-1]
        new_prob_table[:,:] += single_state_probs[3] * prob_table[:,:]
        prob_table = new_prob_table

    return prob_table


def calc_information(loss, n_photons, n_states):
    """
    Calculates the infromation transffered during an integration time.

            Parameters:
                    loss (float): A float between 0 and 1 representing the
                        loss in the setup (the loss in the intensity).
                    n_photons (float): A float representing the average number
                        of photons in the integration time.
                    n_states (int): The number of simulated quantum states 
                        during the integration time.
            Returns:
                    information (float): The calculated information (in bits).
    """
    
    # Creating the measurement probability table for each phase
    z_prob_table = create_probability_table(loss, n_photons, n_states, True)
    pi_prob_table = create_probability_table(loss, n_photons, n_states, False)

    # The total measurement probabilities
    measure_probs = z_prob_table + pi_prob_table

    # Calculate the mutual information between the input and output of the channel
    z_info = z_prob_table * np.log2(z_prob_table / (P_ZERO * measure_probs))
    pi_info = pi_prob_table * np.log2(pi_prob_table / (P_PI * measure_probs))

    information = np.sum(z_info+pi_info)

    return information


def calculate_simulation(data_path, loss_samples=100, photons_samples=100,
    n_states=100, max_photons=10):
    """
    Calculates the integration time information for multiple losses and gains.
    The results are saved in the data dir for future uses.

            Parameters:
                    data_path (str): The path to which the simulation data
                        will be saved.
                    loss_samples (int): The number of losses to calculate.
                    photons_samples (int): The number of photons to calculate.
                    n_states (int): The number of quantum states to simulate.
                    max_photons (float): The maximal average number of photons
                        in integration time.

    """
    # Initializing the calculation space.
    loss_array = np.linspace(EPS, 1-EPS, loss_samples)
    n_photons_array = np.linspace(EPS, max_photons, photons_samples)

    # Initializing the results array.
    info_array = np.zeros((photons_samples, loss_samples))

    # Calculating the information over the calculation space.
    for loss_ind, loss in enumerate(loss_array):
        for alpha_ind, n_photons in enumerate(n_photons_array):
            info_array[alpha_ind, loss_ind] =\
                calc_information(loss, n_photons, n_states)

    # Saves the calculations.
    with open(data_path, 'wb') as f:
        np.save(f, max_photons)
        np.save(f, info_array)


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
    with open(data_path, 'rb') as f:
        max_photons = np.load(f)
        info_array = np.load(f)

    # Plotting the simulation results.
    plt.figure()
    plt.imshow(info_array[::-1, :], extent=[0, 1, 0, max_photons], aspect="auto", 
               cmap='jet', interpolation='bilinear')
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=16)
    plt.xlabel('Loss', fontsize=24)
    plt.ylabel('Number of photons', fontsize=24)
    plt.title('Integration Time Information', fontsize=30)
    plt.xticks(fontsize=16, rotation=0)
    plt.yticks(fontsize=16, rotation=0)
    fig = plt.gcf()
    fig.set_size_inches((11, 8.5), forward=False)

    # Saving the figure to the given path.
    plt.savefig(save_path, dpi=500)
