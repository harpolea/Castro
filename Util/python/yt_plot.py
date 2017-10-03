from yt.mods import *
import pathlib
import sys
import re
from matplotlib import animation, rc_context, cm
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from yt import apply_colormap
from functools import partial
import os
import shutil
from ffmpy import FFmpeg

class Simulation(object):

    """
    A class for plotting castro simulations.
    """

    def __init__(self, location=None, prefix=None):
        """
        Create a class instance. Location is the path to the folder containing the simulation's plot folders.

        Parameters
        ----------
        location : string
            path to the folder containing the simulation's plot folders
        prefix : string
            prefix used to label plot folders
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
                    re.fullmatch(str(self.root_dir) + '/' + '\S+\d{5}', str(f)))])

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
            self.prefix = str(self.root_dir) + '/' + str(self.plot_folders[0])[:-5]
        else:
            # a list of folders containing plot data
            self.plot_folders = sorted([str(f) for f in
                self.root_dir.iterdir()
                if (f.is_dir() and
                    #'.old.' not in str(f) and
                    re.fullmatch(str(self.root_dir) + '/' + prefix + '\d{5}', str(f)))])

            # exit if there are no valid plot folders in the root_dir
            if len(self.plot_folders) == 0:
                sys.exit('No plot folders in this location matching prefix {}'.format(prefix))

            self.prefix = str(self.root_dir) + '/' + prefix

        # load list of fields
        ds = load(self.plot_folders[0])
        self.fields = [f[1] for f in ds.index.field_list]

        # a list of the fields containing the (primitive) velocity
        self.velocity_fields = []

        # time series data
        self.ts = None

    def plot_at_time(self, field, t, save=True, normal_axis='z', velocity_vector=False, grids=False):
        """
        Plot given field at time t. Returns the plot object.

        Parameters
        ----------
        field : string
            name of field to be plotted
        t : int
            number of timestep we wish to plot
        save : bool
            do we save the plot or not?
        normal_axis :
            axis normal to plane of plot
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

        centre = 0.5 * (ds.domain_left_edge + ds.domain_right_edge) * np.sign(ds.domain_width)

        p = SlicePlot(ds, normal_axis, field, origin='native', center=centre)

        if str(ds.domain_width.units) == 'code_length':
            axis_units = 'unitary'
        else:
            axis_units = ds.domain_width.units

        axis_widths = tuple([(ds.domain_width.value[i], axis_units) for i in range(3) if ds.coordinates.axis_id[normal_axis] != i])

        if normal_axis == 'y': # reverse order
            axis_widths = axis_widths[::-1]

        p.set_width(axis_widths)

        if velocity_vector:
            self.add_velocity_vector(p)
        if grids:
            self.add_grids(p)

        if save:
            p.save(plotfilename +  '_' + field + '.png')

        return p

    def plot3d_at_time(self, field, t, iso_value=None, colourmap_field=None, save=True, grids=False):
        """
        Plot given field at time t. Returns the plot object.

        Parameters
        ----------
        field : string
            name of field to be plotted
        t : int
            number of timestep we wish to plot
        save : bool
            do we save the plot or not?
        iso_value : None or float
            value of iso contour. If None, will calculate it as the weighted_average_quantity of the field.
        colourmap_field : string
            name of field for colourmap. If None, use field.
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

        if iso_value is None:
            iso_value = ds.all_data().quantities.weighted_average_quantity(field, 'cell_mass')
        if colourmap_field is None:
            colourmap_field = field

        ds.periodicity = (True, True, True)
        fig = figure()
        ax = fig.gca(projection='3d')
        domain = ds.all_data()
        surface = ds.surface(domain, field, iso_value)
        p3dc = Poly3DCollection(surface.triangles, linewidth=0.)
        try:
            colours = apply_colormap(surface[colourmap_field])
            colours = colours[0,:,:] / 255.0
            colours[:,3] = 0.3
            p3dc.set_facecolors(colours)

            m = cm.ScalarMappable()
            m.set_array(colours)
            cbar = fig.colorbar(m, ax=ax)
        except ValueError:
            pass

        ax.add_collection(p3dc)

        if save:
            fig.savefig(plotfilename + '_3d.png')

        return fig

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

    def add_grids(self, p):
        """
        Add grids to plot.
        """
        p.annotate_grids()
        return p

    def load_time_series(self, force=False, setup_function=None, n_processors=True):
        """
        Load output data for entire time series.

        Parameters
        ----------
        force : bool
            Force reload time series data. Defaults to false, so time series data only loaded if not currently there.
        setup_function : callable, accepts a ds
            Function passed to DatasetSeries constructor which is called whenever a dataset is loaded. If wish to use multiple functions here, recommend putting them inside a wrapper function. For functions with multiple arguments, recommend using partial to set all arguments other than ds (or could again use a wrapper function).
        n_processors : bool or int
            Behaviour when .piter() is called on DatasetSeries object. If True or integer, will be iterated with that integer number of processors assigned to each parameter file provided in loop.
        """
        if self.ts is None or force:
            self.ts = DatasetSeries(self.plot_folders, setup_function=setup_function, parallel=n_processors)

    def plot_animation(self, field, animation_name=None, plot_modifier_functions=[], save=True, normal_axis='z'):
        """
        Create an animation for entire time series for given field and saves to file. Returns animation object

        Parameters
        ----------
        field : string
            Field to be plotted
        animation_name : string
            Name of file to save animation to. Defaults to an mp4
            file with same prefix as dataset.
        plot_modifier_functions : list of callable, accept plot object
            A list of functions to apply to plot
        save : bool
            do we save the animation or not?
        normal_axis :
            axis normal to plane of plot
        """

        if field not in self.fields:
            raise self.InvalidFieldError(field)

        self.load_time_series()

        plot = SlicePlot(self.ts[0], normal_axis, field, origin='native')

        if str(self.ts[0].domain_width.units) == 'code_length':
            axis_units = 'unitary'
        else:
            axis_units = self.ts[0].domain_width.units

        axis_widths = tuple([(self.ts[0].domain_width.value[i], axis_units) for i in range(3) if self.ts[0].coordinates.axis_id[normal_axis] != i])

        if normal_axis == 'y': # reverse order
            axis_widths = axis_widths[::-1]

        plot.set_width(axis_widths)

        # we want the colourbar to stay the same throughout - shall default to using the limits of the first dataset in the time series, but this can be overridden by using the plot modifier functions.
        colourbar_limits = self.ts[0].all_data().quantities.extrema(field)
        plot.set_zlim(field, colourbar_limits[0], colourbar_limits[1])

        for fun in plot_modifier_functions:
            fun(plot)

        fig = plot.plots[field].figure

        # animate must accept an integer frame number. We use the frame number
        # to identify which dataset in the time series we want to load
        def animate(i):
            plot._switch_ds(self.ts[i])

        anim = animation.FuncAnimation(fig, animate, frames=len(self.ts))

        if save:
            if animation_name is None:
                animation_name = self.prefix + '.mp4'

            with rc_context({'mathtext.fontset': 'stix'}):
                anim.save(animation_name)

        return anim

    def parallel_plot_animation(self, field, animation_name=None, plot_modifier_functions=[], save_frames=False, normal_axis='z'):
        """
        Create an animation for entire time series for given field and saves to file. Returns plot object

        Parameters
        ----------
        field : string
            Field to be plotted
        animation_name : string
            Name of file to save animation to. Defaults to an mp4
            file with same prefix as dataset.
        plot_modifier_functions : list of callable, accept plot object
            A list of functions to apply to plot.
        save_frames :
            don't delete frames after making animation
        normal_axis :
            axis normal to plane of plot
        """
        if field not in self.fields:
            raise self.InvalidFieldError(field)

        self.load_time_series()

        plot = SlicePlot(self.ts[0], normal_axis, field, origin='native')
        if str(self.ts[0].domain_width.units) == 'code_length':
            axis_units = 'unitary'
        else:
            axis_units = self.ts[0].domain_width.units

        axis_widths = tuple([(self.ts[0].domain_width.value[i], axis_units) for i in range(3) if self.ts[0].coordinates.axis_id[normal_axis] != i])

        if normal_axis == 'y': # reverse order
            axis_widths = axis_widths[::-1]

        plot.set_width(axis_widths)

        # we want the colourbar to stay the same throughout - shall default to using the limits of the first dataset in the time series, but this can be overridden by using the plot modifier functions.
        colourbar_limits = self.ts[0].all_data().quantities.extrema(field)
        plot.set_zlim(field, colourbar_limits[0], colourbar_limits[1])

        for fun in plot_modifier_functions:
            fun(plot)

        fig = plot.plots[field].figure

        try:
            os.mkdir(str(self.root_dir) + '/frames')
        except FileExistsError:
            pass

        for i, ds in enumerate(self.ts.piter()):
            plot._switch_ds(ds)
            plot.save(str(self.root_dir) + '/frames' + self.plot_folders[i][len(str(self.root_dir)):] + '.png')

        if animation_name is None:
            animation_name = self.prefix + '.mp4'

        # stitch frames together
        ff = FFmpeg(inputs={ str(self.root_dir) + '/frames/*.png': '-framerate 10 -pattern_type glob'},
                    outputs={animation_name : '-y -c:v libx264 -r 10'})

        # run this inside a try catch block so still deletes frames afterwards
        try:
            ff.run()
        except FFRuntimeError:
            print('Could not create animate using ffmpeg command: ')
            print(ff.cmd)

        if not save_frames:
            shutil.rmtree( str(self.root_dir) + '/frames')

    def set_display_name(self, old_name, new_name, ds):
        """
        Sets the name of the variable used in the plot.

        Parameters
        ----------
        old_name : string
            existing field name
        new_name : string
            what the field should be labelled as when we plot it
        ds : Dataset
            dataset whose label we want to change
        """
        # this was far harder to do than it should be.
        # why is there no built in function for this????
        try:
            ds.field_info
        except AttributeError:
            ds.create_field_info()

        ds.field_info[('boxlib', old_name)].display_name = new_name

        if new_name not in self.fields:
            self.fields.append(new_name)

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
