"""
This code calculates and plots the sensitivity to intercept-resend attack.

This calculation is done by calculating the probability for Eve to guess
Alice's basis correctly since if she knows it, she can retransmit the original
state with very high probability (close to 1).
Here we assume that Eve guesses Alice's basis using the following method:
1.  Eve measures Alice's transmission in the original bases, each half of the
    integration time in a different basis (or a quarter of the integration
    time in each phase). She does that since the original bases gives her
    the maximal amount of information (as can be seen in the
    intercept-resend information simulation).
2.  Eve calculates which basis had the highest probability to produce her
    measurement results and guesses this basis. If the two basis had the same
    probability - she guesses a random basis.

This method is also called Maximum Likelihood Estimation (or MLE) for the
suggested measurement method.

Since it is difficult to analytically calculate the probabilities for this
method we give a numerical estimation by simulating 100 quantum states per
integration time. We will also assume without loss of generality that Alice
transmits a state with phi_alice = 0.

This calculation is done under the three assumptions:
1. Low Amplification - we assume that the squeezing in the non-linear crystals
    was low (small alpha) such that each quantum state will contain one
    (signal, idler) pair at most. Here, in contrast to the single state
    information this assumption is automatically applied for a small number of
    photons since we are using a large number of quantum states.
2. Randomaly Chosen Phase - we assume that Alice created the state with a
    random phase between 0, pi/2, pi and 3pi/2 (P = 0.25).
"""

import numpy as np
import matplotlib.pyplot as plt

import itertools as its

import os


ERR = 1e-5


def single_state_measurement_prob(n_photons):
    """
    Calculates the measurement probabilities of a single state.

            Parameters:
                    n_photons (float): The average number of photons in the 
                        quantum state.
            Returns:
                    p_1 (np.array): The probabilities to measure 1 photon in
                        each basis (0, pi/2, pi and 3*pi/2).
    """
    phi_eve = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
    p_1 = n_photons * (np.cos(phi_eve/2) ** 2)
    return p_1


def calculate_probability_table(n_photons, n_states):
    """
    Calculates the measurement probabilities of the entire integration time.

            Parameters:
                    n_photons (float): A float representing the average number
                        of photons in the integration time.
                    n_states (int): The number of simulated quantum states 
                        during the integration time.
            Returns:
                    p_table (np.array): A (n_states X 4) numpy array with
                        the probabilities to measure the different number of
                        photons in each basis.
    """
    # Initialize the probabilities table.
    p_table = np.zeros((n_states//4, 4))
    p_table[0, :] = 1

    # Calculate the single state probabilities.
    p_1 = single_state_measurement_prob(n_photons / n_states)
    p_0 = 1 - p_1

    # Update the probabilities for the different quantum states.
    for basis in range(4):
        for state in range(n_states // 4):
            p_table[1:, basis] =\
                p_0[basis]*p_table[1:, basis] + p_1[basis]*p_table[:-1, basis]
            p_table[0, basis] *= p_0[basis]

    return p_table


def calculate_correct_guess_prob(n_photons, n_states):
    """
    Calculates Eve's probability to guess Alice's basis correctly.

            Parameters:
                    n_photons (float): A float representing the average number
                        of photons in the integration time.
                    n_states (int): The number of simulated quantum states 
                        during the integration time.
            Returns:
                    guess_prob (float): Eve's correct guess probability.
    """

    # Create iterator over all the possible measurements results.
    measurement_ops = its.product(range(n_states // 4), repeat=4)

    # Calculating the probabilities table.
    prob_table = calculate_probability_table(n_photons, n_states)

    # For each measurement option, calculate the probability to measure it in
    # the two bases and take the basis with the higher probability.
    guess_prob = 0
    for state in measurement_ops:

        phase_probs = []
        for guessed_phi in range(4):
            prob = 1
            for ind in range(4):
                prob *= prob_table[state[ind], (ind + guessed_phi) % 4]
            phase_probs.append(prob)

        # The real probability to produce the measured result.
        real_prob = phase_probs[0]

        # The probabilities to produce the measurement state in both bases.
        basis_0_prob = (phase_probs[0] + phase_probs[2]) * 0.5
        basis_1_prob = (phase_probs[1] + phase_probs[3]) * 0.5

        basis_0_prob = round(basis_0_prob / ERR) * ERR
        basis_1_prob = round(basis_1_prob / ERR) * ERR
        
        # Guessing the basis.
        if basis_0_prob > basis_1_prob:
            # Guessed the basis correctly.
            guess_prob += real_prob
        elif basis_0_prob == basis_1_prob:
            # Had 50% to guess the basis correctly.
            guess_prob += real_prob * 0.5

    return guess_prob


def calculate_simulation(data_path, photons_samples=50, n_states=100,
    max_photons=10):
    """
    Calculates the single states information for multiple losses and gains.
    The results are saved in the data dir for future uses.

            Parameters:
                    data_path (str): The path to which the simulation data
                        will be saved.
                    photons_samples (int): The number of photons samples.
                    n_states (int): The number of quantum states to simulate.
                    max_photons (float): The maximal average number of photons
                        in integration time.
    """
    # Initializing the calculation space.
    photons_array = np.linspace(max_photons/photons_samples, max_photons,
        photons_samples)

    # Initializing the results array.
    correct_guess_probs = np.zeros(photons_samples)

    # Calculating the information over the calculation space.
    for ind, photons in enumerate(photons_array):
        correct_guess_probs[ind] =\
            calculate_correct_guess_prob(photons, n_states)

    # Saves the calculations.
    with open(data_path, 'wb') as f:
        np.save(f, photons_array)
        np.save(f, correct_guess_probs)


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
        photons_array = np.load(f)
        correct_guess_probs = np.load(f)

    # Plotting the simulation results.
    plt.figure()
    plt.plot(photons_array, correct_guess_probs)
    plt.title('Sensitivity to Intercept-Resend Attack', fontsize=30)
    plt.xlabel('$N_{\omega}^{max}$', fontsize=24)#Average Number of Photons
    plt.ylabel('$P_{Eve}$', fontsize=24)#Correct Basis Guess Probability
    plt.xticks(fontsize=16, rotation=0)
    plt.yticks(fontsize=16, rotation=0)

    fig = plt.gcf()
    fig.set_size_inches((11, 8.5), forward=False)

    # Saving the figure to the given path.
    plt.savefig(save_path, dpi=500)
