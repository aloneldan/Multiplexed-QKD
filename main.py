"""
This module runs all the different multiplexed qkd simulations.
"""

import simulations.single_state_information as single_state_info
import simulations.integration_time_information as integration_time_info

import os


# Path for the data dir
DATA_PATH = "data"

# Path for the results dir
RESULTS_PATH = 'results'


if __name__ == "__main__":
	single_state_info.simulate(
		data_path=os.path.join(DATA_PATH, 'single_state_info.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'single_state_info.png'))
	integration_time_info.simulate(
		data_path=os.path.join(DATA_PATH, 'integration_time_info.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'integration_time_info.png'))