import pathlib
from validate_inputs import check_parameter_file

def test_check_parameter_file():
    # Tests all existing parameter files are valid

    root_dir = pathlib.Path('../../Exec')

    # trawl through folders to find all files with 'inputs' in their name
    for d in (d for d in root_dir.iterdir() if d.is_dir()):
        for c in (c for c in d.iterdir() if c.is_dir()):
            for f in (f for f in c.iterdir() if f.is_file() and 'inputs' in str(f)):
                print(f)
                if 'test_react' in str(f) or '-Tf0.1' in str(f) or 'Planck' in str(f):
                    continue
                if '3d' in str(f):
                    check_parameter_file(str(f))
                elif '2d' in str(f):
                    check_parameter_file(str(f), dim=2)
                elif '1d' in str(f):
                    check_parameter_file(str(f), dim=1)
                else:
                    try:
                        check_parameter_file(str(f))
                    except (AssertionError, TypeError):
                        try:
                            check_parameter_file(str(f), dim=2)
                        except (AssertionError, TypeError):
                            check_parameter_file(str(f), dim=1)
