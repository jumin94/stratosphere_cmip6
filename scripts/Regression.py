from scipy.stats import linregress
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from abc import ABC, abstractmethod

class Regression:
    """Regressed monthly resolved zonal mean anomalies onto a given yearly index. It saves the slope and corresponding pvalues (
    two-sided Wald test). """

    def __init__(self, model_name, index, anomalies, spatial_coord):
        """
        :param model_name:
        :type model_name: str
        :param index: Some yearly index
        :type index: Xarray array,  (nyears)
        :param anomalies: zonal mean anomalies
        :type anomalies: Xarray array, (ntime, nlat/npres). Coordinate names "time" and "lat"/"pres"
        """

        self._index = index
        self.name = model_name
        self._anomalies = anomalies
        # Check that their is an index entry for every year of data
        assert self.index.size == self.anomalies.shape[0]/12
        
        self.spatial_coord = spatial_coord

        self._slope, self._pvalue = self.regression()

    @property
    def anomalies(self):
        return self._anomalies

    @anomalies.setter
    def anomalies(self, anomalies):
        self._anomalies = anomalies
        self._slope, self._pvalue = self.regression()

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, index):
        self._index = index
        self._slope, self._pvalue = self.regression()

    @property
    def slope(self):
        return self._slope

    @property
    def pvalue(self):
        return self._pvalue

    def regression(self):
        """Regresses the anomalies onto the index."""

        n_spatial = self.anomalies[self.spatial_coord].size

        # Resulting arrays have the shape of (12 month, spatial resolution)
        slope_array = np.zeros((12, n_spatial))
        p_array = np.zeros((12, n_spatial))

        # For every spatial point, group the data by month
        for i in range(n_spatial):
            groups = self.anomalies.isel({self.spatial_coord:i}).groupby("time.month")
            # For every month, regress the anomaly values against the index. The anomaly values have the size nyears,
            # as they are the anomalies for the specific month at the specific latitude.
            for month in range(12):
                slope, intercept, r, p, se = linregress(self.index, groups[month + 1].data)
                slope_array[month, i] = slope
                p_array[month, i] = p

        # Save the slope and pvalue in xarrays datasets.
        slope = xr.DataArray(slope_array, coords=[np.arange(12) + 1, self.anomalies[self.spatial_coord]], dims=["month", self.spatial_coord])
        pvalue = xr.DataArray(p_array, coords=[np.arange(12) + 1, self.anomalies[self.spatial_coord]], dims=["month", self.spatial_coord])

        return slope, pvalue


class RegressionSPV(Regression, ABC):

    def __init__(self, model_name, ua10_path, variable_path, plev=None, variable_name = 'ua', spv_index_latitude=-60, **kwargs):
        """
        :param model_name:
        :type model_name: str
        :param variable_path: path to zonal mean (monthly mean) data
        :type variable_path: str
        :param ua1_path: path to zonal mean monthly mean 1hPa zonal wind data
        :type ua1_path: str
        :param plev: pressure level to use
        :type plev: int
        """
        
        self.spv_index_latitude = spv_index_latitude
        # Read in 1hPa climatology and anomalies
        climatology, anomalies = self.anomalies_read_in(variable_path, plev)
        # Calculate SPV_Index
        spv_index = self.calc_spv_index(ua10_path)
        # Do the regression
        super().__init__(model_name=model_name, index=spv_index, anomalies=anomalies[variable_name], **kwargs)
        self.climatology = climatology[variable_name]
    
    def calc_spv_index(self,datapath):
        """Read in data and calculate the SPV index as the SON mean at 60hPa (or a different supplied latitude)."""
        model = self.xarray_read_in(datapath, plev = 1000)
        model = model.sel(lat=self.spv_index_latitude, method="nearest")
        spv_index = model.groupby("time.season")['SON'].groupby('time.year').mean('time')
        spv_index = spv_index['ua'] - spv_index['ua'].mean()
        return spv_index

    def xarray_read_in(self, datapath, plev=None):
        model = xr.load_dataset(datapath)
        if plev is not None:
            model = model.sel(plev=plev, method="nearest")
        model = model.sel(time=slice("1979-01-01", "2014-12-30"))
        return model

    @abstractmethod
    def anomalies_read_in(self):
        pass


class RegressionU1SPV(RegressionSPV):
    """Subclass of Regression for handling the input for regressing zonal mean anomalies latitudinally resolved onto the
    SPV_Index. At the moment using 1979-2014 is hardcoded."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args, spatial_coord='lat', **kwargs)

    def anomalies_read_in(self, datapath, plev):
        """Read in data and calculate monthly anomalies compared to climatology."""
        model = self.xarray_read_in(datapath, plev=plev)
        if model.coords['lat'][0]<model.coords['lat'][-1]:
            model = model.sel(lat=slice(-90, 0))
        else:
            model = model.sel(lat=slice(0, -90))
        climatolotgy = model.groupby("time.month").mean("time")
        anomalies = model.groupby("time.month") - climatolotgy
        return climatolotgy, anomalies

class RegressionPresSPV(RegressionSPV):
    """Subclass of Regression for handling the input for regressing zonal mean anomalies pressure resolved onto the
    SPV_Index. At the moment using 1979-2014 is hardcoded."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args, spatial_coord='plev', **kwargs)

    def anomalies_read_in(self, datapath, spv_index_lat):
        """Read in data and calculate monthly anomalies compared to climatology."""

        model = self.xarray_read_in(datapath, plev = None)
        model = model.sel(lat=self.spv_index_latitude, method="nearest")
        climatolotgy = model.groupby("time.month").mean("time")
        anomalies = model.groupby("time.month") - climatolotgy
        return climatolotgy, anomalies
    


class RegressionERA5(RegressionSPV):

    def calc_spv_index(self, datapath):
        """Read in data and calculate the SPV index as the SON mean at 60hPa."""
        model = self.xarray_read_in(datapath, plev = None)
        model = model.sel(lat=self.spv_index_latitude, method="nearest")
        spv_index = model.groupby("time.season")['SON'].groupby('time.year').mean('time')
        spv_index = spv_index['u'] - spv_index['u'].mean()
        return spv_index

    def xarray_read_in(self, datapath, plev):
        model = xr.load_dataset(datapath)
        
        if plev is not None:
            model = model.sel(level=plev, method="nearest")
            model = model.rename({'latitude': 'lat',})
            
        model = model.sel(time=slice("1980-01-01", "2014-12-30"))
        return model

class RegressioEra5U1SPV(RegressionERA5,RegressionU1SPV):
    pass

class RegressioEra5PresSPV(RegressionERA5,RegressionPresSPV):
    pass

def plot_slope_data(ds, ax, cbar_label = '', title = None, scicbar=False, spatial_coord='lat'):
    """Plot regression slopes with stippling for significance and climatology contours.

    :param ds: Holds the slope, pvalue and climatology data
    :type ds: Instance of 'Regression' or one of its subclasses
    :param ax:
    :type ax:
    :param scicbar: Whether or not to use scientific notation for the colorbar labels
    :type scicbar: Boolean
    """
    nlevls = 14
    slope = ds.slope
    pvalue = ds.pvalue
    climatology = ds.climatology
    if title is None:
        title = ds.name

    # Check that they all have the same dimensions, i.e. (12xnlat), 12 for the 12 months.
    assert slope.shape == climatology.shape
    assert pvalue.shape == climatology.shape

    # Find the coordinates where the pvalue is smaller than 0.05
    stipp_idx_month, stipp_idx_lat = np.where(pvalue < 0.05)
    stipp_coord_lat = pvalue.coords[spatial_coord][stipp_idx_lat]
    stipp_coord_month = pvalue.coords['month'][stipp_idx_month]
    
    y = slope.coords[spatial_coord].data
    x = slope.month.data
    X, Y = np.meshgrid(x, y)

    absmax = max(abs(slope.min()), abs(slope.max()))
    # Plot the Slope
    cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                                 colors =[(0, 0, 1), 
                                                          (1, 1., 1), 
                                                          (1, 0, 0)],
                                                 N=nlevls-1,)
    
    im = ax.contourf(X, Y, slope.T, cmap=cmap, vmin=-absmax, vmax=absmax, levels=nlevls)
    cb = plt.colorbar(im, ax = ax, label=cbar_label)

    if scicbar:
        # Make scientific colorbar ticks labels
        cb.formatter.set_powerlimits((0, 0))
        cb.ax.yaxis.set_offset_position('right')                         
        cb.update_ticks()

    # Stipv the points that are significant
    ax.scatter(stipp_coord_month, stipp_coord_lat, color='black', s=1.5)

    # Plot climatology contour
    ax.contour(X, Y, climatology.T, colors='k', linewidths=1)
    
    ax.set_xlim(5, 12)

    ax.set_xlabel('Month', fontsize=14)
    if spatial_coord == 'lat.':
        ax.set_ylabel('Latitude', fontsize=14)
    elif spatial_coord == 'plev':
        ax.set_ylabel('Pressure', fontsize=14)
        ax.semilogy()
        ax.set_ylim(100001, 99)
    else:
        ax.set_ylabel(spatial_coord, fontsize=14)
    ax.set_title(title)

