
import xarray as xr
from matplotlib import pyplot as plt
from google.cloud import storage
import fsspec
import matplotlib.animation as animation
import pygmt
import imageio
import datetime
import seaborn as sns
import pandas as pd
import cmaps
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')

bucket_name='rangelands'
bucket = storage_client.get_bucket(bucket_name)

def load_dataset(filename, engine="scipy", *args, **kwargs) -> xr.Dataset:
    """Load a NetCDF dataset from local file system or cloud bucket."""
    with fsspec.open(filename, mode="rb") as file:
        dataset =  xr.load_dataset(file, engine=engine, mask_and_scale=True, *args, **kwargs)
    return dataset
    
    
ndvi_file = 'gs://' + bucket.name + f'/Ranch_Runs/HLR/D12/covariates_nc/D12_ndvi.nc'

ndvi = load_dataset(ndvi_file)
ndvi['ndvi'].values[ndvi['ndvi'].values<0] = 0
ndvi['ndvi'] = ndvi['ndvi'].interpolate_na(dim='time')

ndvi = ndvi.sel(time=slice('2022-03-20', '2022-10-20'))

filenames=[]
for i in range(0, len(ndvi['time'])):

  data = ndvi.sel(time=ndvi['time'][i])
  data_dist = pygmt.grd2xyz(grid=data['ndvi'], output_type="pandas")['ndvi'].astype('float64')
  
  fig, axes = plt.subplots(2,1,figsize=(12,12), gridspec_kw={'height_ratios': [3, 1]})

  axes[0].imshow(data['ndvi'].values, vmin=-.2, vmax = 1.1, cmap=cmaps.turku_map)

  axes[0].text(0.1, 0.1, pd.to_datetime(ndvi['time'][i].values).strftime('%Y-%m-%d'), fontsize=20, horizontalalignment='right',
        verticalalignment='top',
        transform=axes[0].transAxes)
  axes[0].axis('off')
  
  #sns.histplot(x = data_dist, ax=axes[1], binwidth=0.01, element='poly', palette=cmaps.turku_map, legend=False, shade=True)
  sns.kdeplot(x = data_dist, ax=axes[1], palette=cmaps.turku_map, legend=False, fill=True)
  
  # generate a gradient
  x = np.linspace(0,1,100)
  
  im = axes[1].imshow(np.vstack([x,x]), aspect='auto', extent=[*axes[1].get_xlim(), *axes[1].get_ylim()], cmap=cmaps.turku_map, zorder=10)
  path = axes[1].collections[0].get_paths()[0]
  patch = matplotlib.patches.PathPatch(path, transform=axes[1].transData)
  im.set_clip_path(patch)
  
  # Locating current axes
  divider = make_axes_locatable(axes[1])
    
  # creating new axes on the right
  # side of current axes(ax).
  # The width of cax will be 5% of ax
  # and the padding between cax and ax
  # will be fixed at 0.05 inch.
  colorbar_axes = divider.append_axes("right",
                                      size="5%",
                                      pad=0.2)
  cbar = plt.colorbar(im, cax=colorbar_axes)
  
  cbar.ax.tick_params(labelsize=20)
  cbar.ax.set_title('NDVI', fontsize=20)
  
  axes[1].yaxis.set_visible(False)
  axes[1].spines['top'].set_visible(False)
  axes[1].spines['right'].set_visible(False)
  axes[1].spines['left'].set_visible(False)
  
  axes[1].set_ylim(0,6.5)
  axes[1].set_xlim(0,1)
  axes[1].set_xlabel('NDVI', fontsize=18)
  
  filename=f'output/HLR/animation_frames/{i}.jpg'
  filenames.append(filename)
  
  plt.savefig(filename, dpi=300)
  

images=[]
for filename in filenames:
  images.append(imageio.v2.imread(filename))

for _ in range(10):
  images.append(imageio.v2.imread(filename))

imageio.mimsave('output/HLR/HLR_2022_300.gif', images, format='GIF', loop=4,duration=300)

