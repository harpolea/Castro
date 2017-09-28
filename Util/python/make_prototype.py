import re
import numpy as np
import subprocess
import shutil
import sys
from pathlib import Path
from yt_plot import Simulation
from validate_inputs import check_parameter_file

def make_prototype(inputs_filename, dim):
    try:
        infile = open(inputs_filename)
    except:
        sys.exit("error opening the input file")

    check_parameter_file(inputs_filename, dim=dim)

    protofile = open(inputs_filename + '.proto', mode='w')

    lines = infile.readlines()
    set_checkpoint_output = False
    set_plot_output = False
    plot_name = None

    for i, l in enumerate(lines):
        m = re.match('([\w\.]+)\s*\=\s*([\w\s\-\."]+)', l)
        if m:
            if m.group(1) == 'amr.n_cell':
                values = m.group(2).rstrip().split()
                values_copy = values.copy()
                for j, v in enumerate(values_copy):
                    values[j] = int(v)
                values = np.array(values)
                # now we have the values, let's do something with them
                max_log = np.log2(max(values))
                if max_log > 5:
                    values = [int(v) for v in values / 2**(max_log - 5)]

                lines[i] = 'amr.n_cell = ' + str(values)[1:-1] + '\n'
            elif m.group(1) == 'max_step':
                lines[i] = 'max_step = 0\n'
            elif m.group(1) == 'stop_time':
                lines[i] = 'stop_time = 0.0'
            elif m.group(1) == 'amr.checkpoint_files_output' and m.group(2).rstrip() != '0':
                lines[i] = 'amr.checkpoint_files_output = 0\n'
                set_checkpoint_output = True
            elif m.group(1) == 'amr.plot_files_output' and m.group(2).rstrip() != '1':
                lines[i] = 'amr.plot_files_output = 1\n'
                set_plot_output = True
            elif m.group(1) == 'amr.plot_file':
                lines[i] = 'amr.plot_file = proto_' + m.group(2).rstrip() + '\n'
                set_plot_name = True
                plot_name = 'proto_' + m.group(2).rstrip()

    if not set_checkpoint_output:
        lines.append('amr.checkpoint_files_output = 0\n')
    if not set_plot_output:
        lines.append('amr.plot_files_output = 1\n')
    if plot_name is None:
        lines.append('amr.plot_file = proto_plt\n')
        plot_name = proto_plt[:-2]

    # now print to output file
    protofile.writelines(lines)
    protofile.close()
    check_parameter_file(inputs_filename + '.proto')

    return plot_name

def run_prototype(root_dir, executable, inputs_file):
    m = re.match('Castro(\d)d.', executable)
    if m:
        dim = int(m.group(1))
    else:
        dim = 3

    plot_name = make_prototype(root_dir + '/' + inputs_file, dim)
    
    # check to see if output file from Castro exists - if so, delete
    if Path(root_dir + '/' + plot_name + '00000').is_dir():
        shutil.rmtree(root_dir + '/' + plot_name + '00000')

    bash_command = './' + executable + ' ' + inputs_file + '.proto'
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd=root_dir)

    output, error = process.communicate()

    if error is not None:
        print(f'Error!! Code: {error}\nPrinting output:\n')
        print(output.decode())
        sys.exit()

    # plotting
    sim = Simulation(root_dir, plot_name)
    plot = sim.plot_at_time('prim_density', 0)
    #plot_3d = sim.plot3d_at_time('prim_density', 0, None, 'Temp')
