"""
This module runs all the different multiplexed qkd simulations.

For further usage, the simulation data is also saved in a given data dir
and the figures are saved in a given results dir.
"""

import simulations.single_state_information as single_state_info
import simulations.integration_time_information as integration_time_info
import simulations.intercept_resend_information as intercept_resend_info
import simulations.intercept_resend_sensitivity as intercept_resend_sen

import os


# Path for the data dir.
DATA_PATH = "data"

# Path for the results dir.
RESULTS_PATH = "results"

# Running the simulations.
if __name__ == "__main__":
	single_state_info.simulate(
		data_path=os.path.join(DATA_PATH, 'single_state_info.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'single_state_info.png'))
	integration_time_info.simulate(
		data_path=os.path.join(DATA_PATH, 'integration_time_info.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'integration_time_info.png'))
	intercept_resend_info.simulate(
		data_path=os.path.join(DATA_PATH, 'intercept_resend_info.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'intercept_resend_info.png'))
	intercept_resend_sen.simulate(
		data_path=os.path.join(DATA_PATH, 'intercept_resend_sensitivity.npy'),
		save_path=os.path.join(RESULTS_PATH, 'intercept_resend_sensitivity.png'), reprocess=True)