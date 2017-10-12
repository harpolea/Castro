import time
import re
import numpy as np
import subprocess
import shutil
import sys
from pathlib import Path

def guesstimate_runtime(root_dir, executable, inputs_file, max_runtime=60, max_nodes=10, target_runtime=2):
    # max runtime in hours
    processes_per_node = 16

    # run prototype to guesstimate requirements for full scale run

    # first make prototype

    # now run for 5 steps
    n_calls = 5

    m = re.match('Castro(\d)d.', executable)
    if m:
        dim = int(m.group(1))
    else:
        dim = 3

    make_prototype(root_dir, root_dir + '/' + inputs_file, dim, max_step=n_calls)

    # check to see if output file from Castro exists - if so, delete
    if Path(root_dir + '/' + plot_name + '00000').is_dir():
        shutil.rmtree(root_dir + '/' + plot_name + '00000')

    bash_command = './' + executable + ' ' + inputs_file + '.proto'
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd=root_dir)

    start_time = time.time()

    output, error = process.communicate()

    exec_time = (time.time() - start_time) / n_calls

    if error is not None:
        print(f'Error!! Code: {error}\nPrinting output:\n')
        print(output.decode())
        sys.exit()

    # find size and number of steps of actual simulation
    try:
        infile = open(inputs_filename)
    except:
        sys.exit("error opening the input file")

    lines = infile.readlines()
    n_cells = None
    proto_n_cells = None
    max_step = 0
    stop_time = 0.
    for i, l in enumerate(lines):
        m = re.match('([\w\.]+)\s*\=\s*([\w\s\-\."]+)', l)
        if m:
            if m.group(1) == 'amr.n_cell':
                n_cells = m.group(2).rstrip().split()
                n_cells_copy = n_cells.copy()
                for j, v in enumerate(n_cells_copy):
                    n_cells[j] = int(v)
                n_cells = np.array(n_cells)

                max_log = np.log2(max(n_cells))
                if max_log > 5:
                    proto_n_cells = [int(v) for v in n_cells / 2**(max_log - 5)]

            elif m.group(1) == 'max_step':
                max_step = int(m.group(2).rstrip())
            elif m.group(1) == 'stop_time':
                stop_time = float(m.group(2).rstrip())

    # for this to work, need to specify problem runtime by max_step rather than stop_time as no way of finding out how timestep size used in proto. Could do this by rerouting executable's terminal output and regexing that.
    # serial runtime
    time_est = np.prod(n_cells) / np.prod(proto_n_cells) * max_step / n_calls

    scaling_factor = 0.75 # how well does performance scale as we add processes?

    # try one node
    if time_est < 0.2 * target_runtime:
        # one on one process
        n_processes = 1
        n_nodes = 1
    elif time_est / scaling_factor / processes_per_node < 2 * target_runtime:
        n_processes = int(np.ceil(time_est / scaling_factor / target_runtime))
        n_nodes = 1

        time_est = time_est / scaling_factor / n_processes

    else:

        n_processes = int(np.ceil(time_est / scaling_factor / target_runtime))

        if n_processes > processes_per_node * max_nodes:
            n_processes = int(np.ceil(time_est / scaling_factor / max_runtime))

        n_nodes = n_processes // processes_per_node
        n_processes = processes_per_node

        if n_nodes > max_nodes:
            ok = input(f'Number of nodes, {n_nodes}, is greater than maximum number of nodes given by user, {max_nodes}. Shall I continue [y/n]?')

            if ok == 'n':
                sys.exit()
        if time_est / (n_nodes * n_processes) / scaling_factor > max_runtime:
            ok = input(f'Estimated runtime, {time_est}, is greater than maximum runtime given by user, {max_runtime}. Shall I continue [y/n]?')

            if ok == 'n':
                sys.exit()
