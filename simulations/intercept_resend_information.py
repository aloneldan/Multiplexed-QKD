"""
This code calculates and plots the amount of infromation eve can receive by
using intercept-resend attack on a single quantum state.

This calculation is done under the two assumptions:
1. Low Amplification - we assume that the squeezing in the non-linear crystals
    was low (small alpha) such that the quantum state will contain one
    (signal, idler) pair at most. Here we will calculate the information
    for a small average number of photons since we are only interested in
    eve's phase and so this assumption will be valid.
2. Randomaly Chosen Phase - we assume that Alice created the state with a
    random phase between 0, pi/2, pi and 3pi/2 (P = 0.25).
"""


import numpy as np
import matplotlib.pyplot as plt

import os


EPS = 1e-10

RANDOM_PHASE_PROB = 0.25
MAX_EVE_PHASE = np.pi/2


def single_state_information(n_photons, phi_eve):
    """
    Calculates the information eve receives from an intercept-resend attack.

            Parameters:
                    n_photons (float): The average number of photons in
                        Alice's quantum state.
                    phi_eve (float): A float between 0 and 2*pi representing
                        eve's measurement basis (i.e. her phase).
            Returns:
                    information (float): The calculated information (in bits).
    """
    # Calculating the probabilities to measure 0 or 1 photon.
    phi_alice = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
    p_1 = n_photons * (np.cos((phi_eve + phi_alice) / 2) ** 2) 
    p_0 = 1 - p_1

    p_0 *= RANDOM_PHASE_PROB
    p_1 *= RANDOM_PHASE_PROB

    # Calculating the information.
    info_0 = p_0 * np.log2(p_0 / (RANDOM_PHASE_PROB * np.sum(p_0)))
    info_1 = p_1 * np.log2(p_1 / (RANDOM_PHASE_PROB * np.sum(p_1)))

    information = np.sum(info_0) + np.sum(info_1)

    return information


def calculate_simulation(data_path, phase_samples=100, n_photons=0.1):
    """
    Calculates the single states information for multiple losses and gains.
    The results are saved in the data dir for future uses.

            Parameters:
                    data_path (str): The path to which the simulation data
                        will be saved.
                    phase_samples (int): The number of phase samples.
                    n_photons (float): The average number of photons in
                        Alice's quantum state.
    """
    # Initializing the calculation space.
    phases = np.linspace(0, MAX_EVE_PHASE, phase_samples)

    # Initializing the results array.
    information = np.zeros(phase_samples)

    # Calculating the information over the calculation space.
    for phi_ind, phi in enumerate(phases):
        information[phi_ind] = single_state_information(n_photons, phi)

    # Saves the calculations.
    with open(data_path, 'wb') as f:
        np.save(f, n_photons)
        np.save(f, phases)
        np.save(f, information)


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
        n_photons = np.load(f)
        phases = np.load(f)
        information = np.load(f)

    # Plotting the simulation results.
    plt.figure()
    plt.plot(phases, information)
    plt.title('Information vs. Phase for N = {}'.format(n_photons), fontsize=30)
    plt.xlabel('Phase', fontsize=24)
    plt.ylabel('Information [bits]', fontsize=24)
    plt.xticks(ticks=[0, np.pi/8, np.pi/4, 3*np.pi/8, np.pi/2],
        labels=['0', r'$\pi$/8', r'$\pi$/4', r'3$\pi$/8', r'$\pi$/2'])
    plt.xticks(fontsize=16, rotation=0)
    plt.yticks(fontsize=16, rotation=0)
    fig = plt.gcf()
    fig.set_size_inches((11, 8.5), forward=False)

    # Saving the figure to the given path.
    plt.savefig(save_path, dpi=500)


