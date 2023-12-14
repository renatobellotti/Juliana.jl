from typing import List
import matplotlib.pyplot as plt
import numpy as np


class RadarPlotter:

    def __init__(self, fig=None, ax=None, axis_label_size=16,
                 tick_label_size: int = 12,
                 n_skip_ticklabels: int = 0,
                 markersize=7):
        self._names = []
        self._tick_labels = []

        self._fig = fig
        self._ax = ax
        self._initialised = False
        
        self._axis_label_size = axis_label_size
        self._tick_label_size = tick_label_size
        self._n_skip_ticklabels = n_skip_ticklabels
        self._markersize = markersize

    def add_dimension(self, name: str, tick_labels: List[str]):
        '''
        Add an axis to this radar chart.
        '''
        if len(self._tick_labels) > 0:
            assert len(self._tick_labels[0]) == len(tick_labels)
        self._names.append(name)
        self._tick_labels.append(tick_labels)

    def add_marker(self, dimension: int, value: float, text: str = ''):
        '''
        Add a red dot to the given axis at the specified position.
        '''
        value = np.interp(
            value,
            self._tick_labels[dimension],
            self._ax.get_yticks()
        )
        
        d_theta = 2 * np.pi / len(self._names)
        theta = dimension * d_theta
        self._ax.plot(
            (theta),
            (value),
            'o',
            color='red',
            markersize=self._markersize,
        )
        self._ax.annotate(text, (theta, value), color='red')

    def _annotate_dimension(self, dimension):
        '''
        Add tick and axis labels to the figure.
        '''
        theta = dimension * 2 * np.pi / len(self._names)

        if (np.pi/2 <= theta) and (theta <= 3./2.*np.pi):
            ha = 'right'
            va = 'top'
        else:
            ha = 'left'
            va = 'bottom'
    
        ticks = self._tick_labels[dimension]
        for i, (r, label) in enumerate(zip(self._ax.get_yticks(), ticks)):
            if i < self._n_skip_ticklabels:
                continue
            #self._ax.text(theta, i / len(ticks) + 0.1, label)
            if r < 0.05:
                r = 0.075
            self._ax.annotate(
                label,
                (theta, r),
                ha=ha,
                va=va,
                fontsize=self._tick_label_size,
                )

        # Label the axis.
        self._ax.text(
            theta, 
            1.15,
            self._names[dimension],
            ha=ha,
            va=va,
            fontsize=self._axis_label_size,
            #fontweight='bold',
        )

    def _init_plot(self):
        # Code based on the radar chart documentation:
        # https://matplotlib.org/stable/gallery/specialty_plots/radar_chart.html
        if self._fig is None:
            self._fig, self._ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(projection='polar'))

        # Setup grid in polar coordinates.
        n_ticks = len(self._tick_labels[0])
        # Radial ticks.
        self._ax.set_yticks(np.linspace(0, 1, n_ticks))
        self._ax.set_yticklabels([])
        # Angular ticks.
        self._ax.set_xticks(np.linspace(0, 2*np.pi, len(self._names)+1))
        self._ax.set_xticklabels([])

        for i in range(len(self._names)):
            self._annotate_dimension(i)
        
        self._ax.set_ylim([0, 1])
        self._ax.spines['polar'].set_visible(False)
        
        self._initialised = True
        
    def plot_curve(self, values, color, label=None):
        '''
        Plot a curve of the given color into this diagram.
        '''
        if not self._initialised:
            self._init_plot()
    
        scaled_values = []
        for dimension_tick_labels, value in zip(self._tick_labels, values):
            scaled = np.interp(value, dimension_tick_labels, self._ax.get_yticks())
            scaled_values.append(scaled)

        self._ax.plot(
            self._ax.get_xticks(),
            scaled_values + [scaled_values[0]],
            color=color,
            linewidth=3,
            linestyle='-',
            marker='o',
            markersize=7,
            label=label,
        )

    @property
    def fig(self):
        return self._fig
