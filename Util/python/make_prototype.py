def make_prototype(executable, inputs_filename):
    """
    Run a scaled down version of the problem in order to check that the initial data has been set up properly
    """
    try:
        infile = open(inputs_filename)
    except:
        sys.exit("error opening the input file")

    protofile = open(inputs_filename + '.proto', mode='w')

    lines = infile.readlines()
    set_checkpoint_output = False

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
            elif m.group(1) == 'amr.plot_files_output' and m.group(2).rstrip() != '1':
                lines[i] = 'amr.plot_files_output = 1\n'
            elif m.group(1) == 'amr.plot_file':
                lines[i] = 'amr.plot_file = proto_' + m.group(2).rstrip() + '\n'


    # now print to output file

    protofile.writelines(lines)
