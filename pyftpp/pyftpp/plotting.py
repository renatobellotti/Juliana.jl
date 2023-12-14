import logging
from PIL import Image, ImageDraw
from typing import List, Tuple
from ipywidgets import IntSlider, RadioButtons, interact
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, from_levels_and_colors
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
from .ct import CT
from .dose import Dose
from .dvh import DVH
from .grid import Grid
from .interpolation import interpolate_data
from .structure import Structure


# Code adapted from:
# https://stackoverflow.com/a/10967471
class CTPlotter:

    DEG_TO_RAD_FACTOR = 2*np.pi / 360
    TRANSPARENT_DOSE_THRESHOLD = 0.01 # unit: Gy

    def __init__(self,
                 ct: CT,
                 dose1: Dose = None,
                 target_dose: float = None,
                 dose2: Dose = None,
                 crosshair_visible: bool = True,
                 print_legend: bool = True,
                 dose_labels=('dose1', 'dose2'),
                 dose_levels: List[float] = [],
                 zoom = None):
        '''
        Parameters
        ----------
        ct: pyftpp.CT
            The CT to plot.
        dose1: pyftpp.Dose
            Optional. If given, also display the dose distribution.
        target_dose: float
            Optional. Target dose to use for the color wash. Must be present
            if at least one dose distribution is given.
        print_legend: bool
            Whether to print the legend.
        dose2: pyftpp.Dose
            Optional. If given, display the dose distribution dose1 - dose2.
        dose_levels: list of float or None
            Edges of the bins for the color wash in Gy. Choose the edges automatically if None.
        zoom: dict
            Example:
            {
                "horizontal": [xmin, xmax],
                "vertical": [ymin, ymax],
            }
        '''
        self._dose_labels = dose_labels
        if (dose_levels is not None) and (len(dose_levels) == 0) and (dose1 is not None):
            # Default dose levels are taken from PSIPlan.
            self._dose_levels = np.array([0, 10, 50, 80, 90, 95, 100, 105, 109]) * target_dose / 100.
        else:
            self._dose_levels = dose_levels
        self._zoom = zoom
        self._initialise_grids(ct, dose1, dose2)

        self._cbar = None
        self._target_dose = target_dose

        self._structures = []
        self._colors = {}

        self._fig, self._ax = plt.subplots(1, 1,
                                           figsize=(4, 2.75),
                                           constrained_layout=True,
                                           dpi=200)

        # Cross-hair position.
        self._default_x = -1
        self._default_y = -1
        self._default_z = -1
        self._crosshair_visible = crosshair_visible
        self._print_legend = print_legend

        # Beam directions to plot.
        self._beam_angles = []

        self.unit = '%'

    def _initialise_grids(self, ct, dose1, dose2):
        # Make sure all grids are quadratic and share the same origin.
        assert ct.spacing[0] == ct.spacing[1]
        if dose1 is not None:
            assert np.isclose(ct.origin, dose1.origin, atol=1e-6).all()
            assert dose1.spacing[0] == dose1.spacing[1]
        if dose2 is not None:
            assert np.isclose(ct.origin, dose2.origin, atol=1e-6).all()
            assert dose2.spacing[0] == dose2.spacing[1]

        # Get finest grid resolution.
        finest_grid = ct.grid
        if (dose1 is not None) and dose1.grid.is_finer_than(finest_grid):
            finest_grid = dose1.grid
        if (dose2 is not None) and dose2.grid.is_finer_than(finest_grid):
            finest_grid = dose2.grid

        # Interpolate all grids to the finest resolution.
        if not np.isclose(ct.spacing, finest_grid.spacing, atol=1e-6).all():
            logging.warning('Have to interpolate the CT before plotting it; might be buggy!')
            ct = CT(
                interpolate_data(ct.data, ct.grid, finest_grid),
                finest_grid.spacing,
                finest_grid.origin,
                ct.orientation,
            )
        if (dose1 is not None) and (not np.isclose(dose1.spacing, finest_grid.spacing, atol=1e-6).all()):
            logging.warning('Have to interpolate a dose before plotting it; might be buggy!')
            dose1 = Dose(
                interpolate_data(dose1.data, dose1.grid, finest_grid),
                finest_grid.spacing,
                finest_grid.origin,
            )
        if (dose2 is not None) and (not np.isclose(dose2.spacing, finest_grid.spacing, atol=1e-6).all()):
            logging.warning('Have to interpolate a dose before plotting it; might be buggy!')
            dose2 = Dose(
                interpolate_data(dose2.data, dose2.grid, finest_grid),
                finest_grid.spacing,
                finest_grid.origin,
            )

        # Initialise the grids for this class.
        self._ct = ((ct.data + 1000) / 4096) * 255
        self._ct_origin = ct.origin
        self._ct_spacing = ct.spacing

        self._dose1 = None
        self._dose2 = None
        if dose1 is not None:
            self._dose1 = dose1.data
        if dose2 is not None:
            self._dose2 = dose2.data

    @property
    def fig(self):
        return self._fig

    def add_structure(self, structure: Structure, color: Tuple[int]):
        points = structure.points
        points = np.round((points - self._ct_origin) / self._ct_spacing)
        grid = Grid(
            self._ct.data.shape,
            np.array((1.0, 1.0, 1.0)),
            (0.0, 0.0, 0.0),
        )
        structure = Structure(structure.name, points, grid)
        self._structures.append(structure)
        self._colors[structure.name] = color

    def add_beam(self,
                 gantry_angle: float,
                 couch_angle: float,
                 isocentre: np.ndarray,
                 color: Tuple[float]):
        '''
        Parameters
        ----------
        gantry_angle: float
            Gantry angle in degrees.
        couch_angle: float
            Couch angle in degrees.
        isocentre: np.ndarray
            Isocentre of the beam. Shape: (3,)
        color: tuple of float
            rgb values in [0, 1).
        '''
        self._beam_angles.append((gantry_angle, couch_angle, isocentre, color))

    def set_crosshair(self, x, y, z):
        self._default_x = x
        self._default_y = y
        self._default_z = z

    def plot(self, x, y, z, unit, which='no_dose'):
        self._ax.clear()

        # Draw the CT.
        img = Image.fromarray(self._ct[:, :, z].T).convert('RGBA')

        # Draw the contours.
        self._draw_structures(img, z)

        # Draw the dose distribution.
        dose = self._get_dose(which)

        if (dose is not None) and (which != 'no_dose'):
            dose_img = self._draw_dose(dose[:, :, z].T, which, unit)
            img = Image.blend(img, dose_img, alpha=0.5)

        # Draw the entire picture.
        if self._zoom is not None:
            xmin, xmax = self._zoom['horizontal']
            ymin, ymax = self._zoom['vertical']
            img = img.crop((xmin, ymin, xmax, ymax))
        self._ax.imshow(img)

        # Draw the cross hair.
        if self._crosshair_visible:
            self._ax.hlines(y=x, xmin=0, xmax=self._ct.shape[1]-1)
            self._ax.vlines(x=y, ymin=0, ymax=self._ct.shape[0]-1)

        return self._fig

    def _draw_structures(self, img, z):
        # CAUTION!
        # The contours are drawn as x: horizontal, y: vertical.
        # Images are drawn as axis0: vertical, axis1: horizontal.
        draw = ImageDraw.Draw(img)
        handles = []
        for structure in self._structures:
            name = structure.name
            points = structure.points

            points = points[np.where(points[:, 2] == z)]
            if points.shape[0] > 0:
                points = np.vstack([points, points[0, :]])
            draw.line(
                list(map(tuple, points[:, :2].tolist())),
                fill=self._colors[name]
            )
            c = [channel / 255 for channel in self._colors[name]]
            handles.append(Patch(color=c, label=name))
        handles = sorted(handles, key=lambda h: h.get_label())
        if self._print_legend:
            self._ax.legend(handles=handles, fontsize=3)

    def _get_dose(self, which):
        if which == 'no_dose':
            dose = None
        elif (which == 'dose1') or (which == 'dose'):
            # 'dose' is used when there is only one dose distribution.
            dose = self._dose1
        elif which == 'dose2':
            dose = self._dose2
        elif which == 'difference':
            dose = self._dose1 - self._dose2
        else:
            raise RuntimeError(f'Unexpected argument: {which}')

        return dose
    
    def _get_cmap(self, which):
        n_colors = 10 if which == 'difference' else 12
        if which == 'no_dose':
            color_wash = None
        elif (which == 'dose1') or (which == 'dose'):
            # 'dose' is used when there is only one dose distribution.
            color_wash = cm.get_cmap('coolwarm', n_colors)
        elif which == 'dose2':
            color_wash = cm.get_cmap('coolwarm', n_colors)
        elif which == 'difference':
            color_wash = cm.get_cmap('coolwarm', n_colors)
        else:
            raise RuntimeError(f'Unexpected argument: {which}')
        
        return color_wash

    def _draw_dose(self, dose, which, unit):
        if which not in ['dose', 'dose1', 'dose2', 'difference']:
            raise RuntimeError('Cannot draw dose distribution! (violated precondition)')

        # Dose values smaller than the threshold will not be displayed.
        too_small = np.where(
            np.abs(dose) < self.TRANSPARENT_DOSE_THRESHOLD
        )

        if unit == '%':
            dose_img_matrix = dose / self._target_dose * 100

            if which == 'difference':
                min_dose = -50
                max_dose = 50
            else:
                min_dose = 0.
                max_dose = 120.

            unit_label = '% of prescribed dose'
        elif unit == 'Gy':
            dose_img_matrix = dose

            if which == 'difference':
                min_dose = -50 * self._target_dose / 100.
                max_dose = 50 * self._target_dose / 100.
            else:
                min_dose = 0.
                max_dose = 120 * self._target_dose / 100.

            unit_label = 'Gy'
        else:
            raise ValueError(f'Invalid unit: {unit}')

        if (self._dose_levels is None) or (which == 'difference'):
            cmap = self._get_cmap(which)
            norm = Normalize(
                vmin=min_dose,
                vmax=max_dose,
            )
        else:
            if unit == 'Gy':
                levels = self._dose_levels
            elif unit == '%':
                levels = self._dose_levels / self._target_dose * 100.
            else:
                raise RuntimeError(f'Unexpected unit: {unit}')
            colors = sns.color_palette(n_colors=len(levels)-1)
            cmap, norm = from_levels_and_colors(levels, colors)

        # Scale to [0, 1] and apply the color wash.
        dose_img_matrix = norm(dose_img_matrix, clip=True)
        dose_img = cmap(dose_img_matrix, bytes=True)
        # Make cells receiving than the treshold dose transparent.
        dose_img[too_small] = 0

        # Draw the color map scale.
        # Avoid duplicating the color bar on every change of the slider.
        mappable = ScalarMappable(norm=norm, cmap=cmap)
        if self._cbar is None:
            self._cbar = self.fig.colorbar(mappable, ax=self._ax)
        else:
            self._cbar.update_normal(mappable)

        if which == 'dose':
            label = f'Dose [{unit_label}]'
        elif which == 'dose1':
            label = f'{self._dose_labels[0]} Dose [{unit_label}]'
        elif which == 'dose2':
            label = f'{self._dose_labels[1]} Dose [{unit_label}]'
        elif which == 'difference':
            label = f'({self._dose_labels[0]} - {self._dose_labels[1]}) Dose [{unit_label}]'
        self._cbar.set_label(label, fontsize=6)

        return Image.fromarray(dose_img)

    def _ipython_display_(self):
        x_val = self._default_x
        y_val = self._default_y
        z_val = self._default_z

        if self._default_x == -1:
            # Below, the x and y axis are exchanged on purpose.
            # See CAUTION comment in CTPlotter.plot().
            if len(self._structures) > 0:
                # Default:
                # Center the cross hair to the center of the first structure.
                points = self._structures[0].points
                mean = np.mean(points, axis=0)
                y_val, x_val, z_val = mean.tolist()
            else:
                x_val = self._ct.shape[1] // 2
                y_val = self._ct.shape[0] // 2
                z_val = self._ct.shape[2] // 2

        # Add interactive controls.
        x = IntSlider(name='x', value=x_val, min=0, max=self._ct.shape[0])
        y = IntSlider(name='y', value=y_val, min=0, max=self._ct.shape[1])
        z = IntSlider(name='z', value=z_val, min=0, max=self._ct.shape[2])

        if self._dose1 is not None:
            if self._dose2 is None:
                options = ['no_dose', 'dose']
                default = 'dose'
            else:
                options = ['no_dose', 'dose1', 'dose2', 'difference']
                default = 'difference'
        else:
            options = ['no_dose']
            default = 'no_dose'

        which = RadioButtons(options=options, value=default)
        unit = RadioButtons(
            options=['Gy', '%'],
            value=self.unit,
        )

        return interact(self.plot, x=x, y=y, z=z, unit=unit, which=which)


class DvhPlotter:

    def __init__(self, full_dose, title, fig=None, ax=None, percent: bool = True):
        self._full_dose = full_dose
        self._title = title
        self._percent = percent

        # Data to plot.
        self._dvhs = {}
        self._colors = {}
        self._markers = []
        self._vlines = []

        # Dose range to plot.
        self._dose_range = np.arange(0., 1.5, 0.01) * full_dose

        # Create and style the figure.
        if fig is None and ax is None:
            fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)

        self._fig = fig
        self._ax = ax


    def add_curve(self, curve_label: str, dvh: DVH, color):
        self._dvhs[curve_label] = dvh
        self._colors[curve_label] = color

    def add_marker(self,
                   d: float,
                   v: float,
                   color,
                   text: str = '',
                   direction='UPPER'):
        '''
        Parameters
        ----------
        d: float
            Dose in Gy.
        v: float
            Volume as a fraction in [0, 1].
        text: str
            Text to draw as the annotation.
        color:
            Color for the marker.
        direction: str
            From where the arrow should point.
            Currently, only UPPER is supported.
        '''
        if direction != 'UPPER':
            raise NotImplementedError('Only "UPPER" markers are supported.')

        self._markers.append({
            'dose': d,
            'volume': v,
            'color': color,
            'text': text,
            'direction': direction,
        })
    
    def add_vertical_line(self, d: float, color):
        '''
        Parameters
        ----------
        d: float
            Dose in Gy at which to display the vertical line.
        '''
        self._vlines.append({
            'dose': d,
            'color': color,
        })

    def plot(self):
        for label in self._dvhs.keys():
            self._plot_curve(label)

        for marker in self._markers:
            self._plot_marker(**marker)
        
        for vline in self._vlines:
            self._plot_vline(**vline)

        self._ax.spines['top'].set_visible(False)
        self._ax.spines['right'].set_visible(False)
        if self._percent:
            self._ax.set_xlim([0, 110])
            self._ax.set_xticks(range(0, 120, 10))
        self._ax.set_ylim([0., 101])
        self._ax.grid(True)
        self._ax.legend(fontsize=16)
        if self._percent:
            xlabel = 'Dose [% of prescribed dose]'
        else:
            xlabel = 'Dose [Gy]'
        self._ax.set_xlabel(xlabel, fontsize=20)
        self._ax.set_ylabel('Volume [%]', fontsize=20)
        self._ax.tick_params(labelsize=16)
        self._ax.set_title(self._title, fontsize=16)
        self._fig.tight_layout()
        return self._fig

    def _plot_marker(self, dose, volume, color, text, direction):
        if direction != 'UPPER':
            raise NotImplementedError('Only "UPPER" markers are supported.')

        # See:
        # https://stackoverflow.com/a/27611041
        self._ax.annotate(
            text,
            xy=(dose, volume),
            xytext=(dose+15, volume + 15),
            arrowprops={
                'width': 1,
                'color': color,
            },
        )

    def _plot_vline(self, dose, color):
        if self._percent:
            dose = dose / self._full_dose * 100.
        
        self._ax.axvline(dose, color=color)

    def _plot_curve(self, label):
        dose = self._dose_range
        volume = self._dvhs[label].V(self._dose_range) * 100.
        positive = volume > 0

        if self._percent:
            dose = dose / self._full_dose * 100.

        self._ax.plot(
            dose[positive],
            volume[positive],
            label=label,
            linewidth=3,
            color=self._colors[label]
        )
