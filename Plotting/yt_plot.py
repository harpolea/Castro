from yt.mods import *
import pathlib
import sys
import re
from matplotlib import animation
from matplotlib import rc_context

class Simulation(object):

    """
    A class for plotting castro simulations.
    """

    def __init__(self, location=None, prefix=None):
        """
        Create a class instance. Location is the path to the folder containing the simulation's plot folders.
        """
        if location is None:
            self.root_dir = pathlib.Path('.')
        else:
            self.root_dir = pathlib.Path(location)

        if prefix is None:
            # a list of folders containing plot data
            self.plot_folders = sorted([str(f) for f in
                self.root_dir.iterdir()
                if (f.is_dir() and
                    '.old.' not in str(f) and
                    re.fullmatch('\S+\d{5}', str(f)))])

            # exit if there are no valid plot folders in the root_dir
            if len(self.plot_folders) == 0:
                sys.exit('No plot folders in this location')

            # check to see if there are multiple sets of plot folders
            # - if so, just take the alphabetically first one
            # (printing a warning to screen to indicate have done so)
            m = sorted(list(set(re.findall('(\S+)\d{5}', ' '.join(self.plot_folders)))))

            if len(m) > 1:
                print("Found multiple sets of plot folders in location: {}".format(', '.join(m)))

                print('Using the first set of plot folders alphabetically with prefix {}'.format(m[0]))

                for i, f in enumerate(self.plot_folders):
                    if not re.fullmatch(m[0] + '\d{5}', f):
                        self.plot_folders = self.plot_folders[:i]
                        break

            # prefix to the plot folder names
            self.prefix = str(self.plot_folders[0])[:-5]
        else:
            # a list of folders containing plot data
            self.plot_folders = sorted([str(f) for f in
                self.root_dir.iterdir()
                if (f.is_dir() and
                    '.old.' not in str(f) and
                    re.fullmatch(prefix + '\d{5}', str(f)))])

            # exit if there are no valid plot folders in the root_dir
            if len(self.plot_folders) == 0:
                sys.exit('No plot folders in this location matching prefix {}'.format(prefix))

            self.prefix = prefix

        # load list of fields
        ds = load(self.plot_folders[0])
        self.fields = [f[1] for f in ds.index.field_list]

        # a list of the fields containing the (primitive) velocity
        self.velocity_fields = []

        # time series data
        self.ts = None

    def plot_at_time(self, field, t, save=True):
        """
        Plot given field at time t.
        """
        if field not in self.fields:
            raise self.InvalidFieldError(field)

        plotfilename = self.prefix + '{:05d}'.format(t)

        if (plotfilename) not in self.plot_folders:
            raise self.InvalidTimeError(t)

        # check to see if
        if self.ts is None:
            ds = load(plotfilename)
        else:
            ds = self.ts[self.plot_folders.index(plotfilename)]

        p = SlicePlot(ds, 'z', field, origin='native')
        p.set_width(1., 'unitary')

        if save:
            p.save(plotfilename + '.png')

        return p

    def set_velocity_fields(self, velocity_fields):
        """
        Set the field names corresponding to the velocity fields.
        """
        for f in velocity_fields:
            if f not in self.fields:
                raise self.InvalidFieldError(f)

        self.velocity_fields = velocity_fields

    def add_velocity_vector(self, p):
        """
        Add velocity vector to plot p.
        """
        p.annotate_quiver(self.velocity_fields[0], self.velocity_fields[1], factor=32)

        return p

    def load_time_series(self):
        """
        Load output data for entire time series.
        """
        self.ts = DatasetSeries(self.plot_folders)

    def plot_animation(self, field, animation_name=None, plot_modifier_functions=[]):
        """
        Create an animation for entire time series for given field and saves to file.

        Parameters
        ----------
        field :
            Field to be plotted
        animation_name :
            Name of file to save animation to. Defaults to an mp4
            file with same prefix as dataset.
        plot_modifier_functions :
            A list of functions to apply to plot.
        """
        if self.ts is None:
            self.load_time_series()

        plot = SlicePlot(self.ts[0], 'z', field, origin='native')
        plot.set_width(1., 'unitary')
        for fun in plot_modifier_functions:
            fun(plot)

        fig = plot.plots[field].figure

        # animate must accept an integer frame number. We use the frame number
        # to identify which dataset in the time series we want to load
        def animate(i):
            plot._switch_ds(self.ts[i])

        anim = animation.FuncAnimation(fig, animate, frames=len(self.ts))

        if animation_name is None:
            animation_name = self.prefix + '.mp4'

        with rc_context({'mathtext.fontset': 'stix'}):
            anim.save(animation_name)

        return anim

    class InvalidFieldError(Exception):
        # Raise if field does not exist in dataset
        def __init__(self, field):
            print("Field {} is invalid".format(field))
            self.message = "Field {} is invalid".format(field)

    class InvalidTimeError(Exception):
        # Raise if try to plot output data for a time for which no valid output data exists
        def __init__(self, t):
            print("No output data exists at time {}".format(t))
            self.message = "No output data exists at time {}".format(t)
