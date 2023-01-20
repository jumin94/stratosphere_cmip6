from scipy.stats import linregress
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt

class Regression:
    """Regressed monthly resolved zonal mean anomalies onto a given yearly index. It saves the slope and corresponding pvalues (
    two-sided Wald test). """

    def __init__(self, model_name, index, anomalies):
        """
        :param model_name:
        :type model_name: str
        :param index: Some yearly index
        :type index: Xarray array,  (nyears)
        :param anomalies: zonal mean anomalies
        :type anomalies: Xarray array, (ntime, nlat). Coordinate names "time" and "lat"
        """

        self._index = index
        self.name = model_name
        self._anomalies = anomalies
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

        n_lat = self.anomalies['lat'].size

        # Resulting arrays have the shape of (12 month, latitudinal resolution)
        slope_array = np.zeros((12, n_lat))
        p_array = np.zeros((12, n_lat))

        # For every latitudinal point, group the data by month
        for i in range(n_lat):
            groups = self.anomalies.isel(lat=i).groupby("time.month")
            # For every month, regress the anomaly values against the index. The anomaly values have the size nyears,
            # as they are the anomalies for the specific month at the specific latitude.
            for month in range(12):
                slope, intercept, r, p, se = linregress(self.index, groups[month + 1].data)
                slope_array[month, i] = slope
                p_array[month, i] = p

        # Save the slope and pvalue in xarrays datasets.
        slope = xr.DataArray(slope_array, coords=[np.arange(12) + 1, self.anomalies['lat']], dims=["month", "lat"])
        pvalue = xr.DataArray(p_array, coords=[np.arange(12) + 1, self.anomalies['lat']], dims=["month", "lat"])

        return slope, pvalue

class RegressionU1SPV(Regression):
    """Subclass of Regression for handling the input for regressing zonal mean zonal wind anomalies at 1hPa onto the
    SPV_Index. At the moment using 1979-2014 is hardcoded."""
    def __init__(self, model_name, ua10_path, ua1_path):
        """
        :param model_name:
        :type model_name: str
        :param ua10_path: path to zonal mean (monthly mean) 10hPa zonal wind data
        :type ua10_path: str
        :param ua1_path: path to zonal mean monthly mean 1hPa zonal wind data
        :type ua1_path: str
        """
        # Read in 1hPa climatology and anomalies
        climatology, anomalies = RegressionU1SPV.anomalies_read_in(ua1_path)
        # Calculate SPV_Index
        spv_index = RegressionU1SPV.calc_spv_index(ua10_path)
        # Do the regression
        super().__init__(model_name=model_name, index=spv_index, anomalies=anomalies['ua'])
        self.climatology = climatology['ua']

    @staticmethod
    def anomalies_read_in(datapath):
        """Read in data and calculate monthly anomalies compared to climatology."""
        model = RegressionU1SPV.xarray_read_in(datapath)
        model = model.sel(lat=slice(0, -90))
        climatolotgy = model.groupby("time.month").mean("time")
        anomalies = model.groupby("time.month") - climatolotgy
        return climatolotgy, anomalies

    @staticmethod
    def calc_spv_index(datapath):
        """Read in data and calculate the SPV index as the SON mean at 60hPa."""
        model = RegressionU1SPV.xarray_read_in(datapath)
        model = model.sel(lat=-60, method="nearest")
        spv_index = model.groupby("time.season")['SON'].groupby('time.year').mean('time')
        spv_index = spv_index['ua'] - spv_index['ua'].mean()
        return spv_index

    @staticmethod
    def xarray_read_in(datapath):
        model = xr.load_dataset(datapath)
        model = model.isel(plev=0)
        model = model.sel(time=slice("1979-01-01", "2014-12-30"))
        return model

def plot_slope_data(ds, ax):
    """Plot regression slopes with stippling for significance and climatology contours.

    :param ds: Holds the slope, pvalue and climatology data
    :type ds: Instance of 'Regression' or one of its subclasses
    :param ax:
    :type ax:
    """
    slope = ds.slope
    pvalue = ds.pvalue
    climatology = ds.climatology

    # Check that they all have the same dimensions, i.e. (12xnlat), 12 for the 12 months.
    assert slope.shape == climatology.shape
    assert pvalue.shape == climatology.shape

    # Find the coordinates where the pvalue is smaller than 0.05
    stipp_idx_month, stipp_idx_lat = np.where(pvalue < 0.05)
    stipp_coord_lat = pvalue.coords['lat'][stipp_idx_lat]
    stipp_coord_month = pvalue.coords['month'][stipp_idx_month]

    y = slope.lat.data
    x = slope.month.data
    X, Y = np.meshgrid(x, y)

    absmax = max(slope.min(), slope.max(), key=abs)
    # Plot the Slope
    im = ax.contourf(X, Y, slope.T, cmap='bwr', vmin=-absmax, vmax=absmax, levels=14)
    plt.colorbar(im, ax = ax)

    # Stipple the points that are significant
    ax.scatter(stipp_coord_month, stipp_coord_lat, color='black', s=1)

    # Plot climatology contour
    ax.contour(X, Y, climatology.T, colors='k')
    
    ax.set_xlim(5, 12)

    ax.set_xlabel('Month', fontsize=14)
    ax.set_ylabel('Latitude', fontsize=14)
    ax.set_title(ds.name)

