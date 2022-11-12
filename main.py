import simulations.single_state_information as single_state_information

import os

# Path for the data dir
DATA_PATH = "data"

# Path for the results dir
RESULTS_PATH = 'results'

if __name__ == "__main__":
	single_state_information.simulate(
		data_path=os.path.join(DATA_PATH, 'single_state_information.npy'), 
		save_path=os.path.join(RESULTS_PATH, 'single_state_information.png'))