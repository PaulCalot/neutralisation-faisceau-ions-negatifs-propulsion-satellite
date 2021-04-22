from main import simulate, processing_only

def simulate_perf_analysis(system):
    if(system == 'square'):
        system_cfg_path = 'cfg_files/tests/tube.cfg'
        simulation_cfg_path = 'cfg_files/tests/simulation_tube.cfg'
        processing_simulation_path = 'cfg_files/tests/processing_tube.cfg'
    elif(system == 'thruster'):
        system_cfg_path = 'cfg_files/tests/thruster.cfg'
        simulation_cfg_path = 'cfg_files/tests/simulation_thruster.cfg'
        processing_simulation_path = 'cfg_files/tests/processing_thruster.cfg'
    simulate(system_cfg_path, simulation_cfg_path, processing_simulation_path)

if __name__ == '__main__':
    # simulate_perf_analysis('square')
    simulate_perf_analysis('thruster')