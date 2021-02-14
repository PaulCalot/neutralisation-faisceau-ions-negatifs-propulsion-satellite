from main import simulate, processing_only

system_cfg_path = ['cfg_files/square.cfg']
simulation_cfg_path = ['cfg_files/simulation.cfg']
processing_simulation_path = ['cfg_files/processing.cfg']

for system, simu, proc in zip(system_cfg_path, simulation_cfg_path, processing_simulation_path):
    simulate(system, simu, proc)
