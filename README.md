# Multiplexed_QKD

This repository is an attachment to my academic work about multiplexed qkd. In my work, I proposed a new QKD protocol based on an unseeded SU(1,1) interferometer, similar to BB84. This repository contains 4 main simulation codes that were used during my work to further analyze the properties and the capabilities of the new suggested protocol, and their results. These simulations are:
1. *single_state_information* - Calculation of the information transfered by a single quantum state during the protocol.
2. *integration_time_information* - Calculation of the information transfered during an integration time of the detectors using the protocol and their dependency on the average number of photons and the loss in the setup.
3. *intercept_resend_information* - Calculation of the information that an attacker, Eve, can extract about the original state from an intercept-resend attack.
4. *intercept_resend_sensitivity* - Calculation of the probability of an attacker to guess correctly Alice's basis by using an intercept-resend attack for a different average numbers of photons.

Notice that simulations 2 and 4 complement each other and allows for further optimization of the protocol as one can investigate the affects of varying the average number of photons on both the protocol efficiency and its sensitivity to attacks.

For convinience, this repository is divided into 3 directories:
* *simulations* - The code directory in which the 4 simulation codes.
* *data* - The simulations raw results. One can use it to work with the simulation data without rerunning the simulations (although it takes only several minutes per simulation for the current parameters).
* *results* - The resulting figures.
