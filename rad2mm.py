from lics_unwrap import *
import framecare as fc
from scipy import ndimage 
from geocube.api.core import make_geocube
import geopandas as gpd
from functions import *
import daz_lib_licsar as dl
import LiCSAR_misc as misc
import sys
import rasterio
from rasterio.merge import merge

'''
This code helps to change BOI radian values to mm values regarding per swath.
Nergizci, Lazecky 28/09/23
'''

if len(sys.argv) < 3:
    print('Please provide frame and pair information: i.e python rad2mm.py 021D_05266_252525 20230129_20230210')
    sys.exit(1)

##variables
frame=sys.argv[1]
pair=sys.argv[2]
tr = int(frame[:3])

#variables from variables
batch=os.environ['BATCH_CACHE_DIR']
tif=os.path.join(batch,frame,'GEOC',pair,pair+'.geo.bovldiff.adf.tif')
outtif=os.path.join(batch,frame,'GEOC',pair,pair+'.geo.bovldiff.adf.mm.tif')
metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
primepoch = misc.grep1line('master=',metafile).split('=')[1]
path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)

##bursts geojson
gpd_bursts = fc.frame2geopandas(frame, use_s1burst=True)
gpd_bursts.to_file(frame+'.geojson', driver='GeoJSON')

###Lazecky conncomps idea
bovlpha=load_tif2xr(tif)
bovlpha = bovlpha.where(bovlpha!=0) # to get rid of zero --now will be nan
aa = bovlpha.copy()
npa = aa.values
npa[npa==0] = np.nan
mask = ~np.isnan(npa)
conncomps, ncomp = ndimage.label(mask) 
aa.values = conncomps
conncomps = aa #just for clarity
conncomps = conncomps.where(~np.isnan(bovlpha))

#geopandas dataframe of burst overlap from functions lib.
gpd_overlaps, overlap_gdf1, overlap_gdf2, overlap_gdf3 = extract_burst_overlaps(frame)

#calculate dfDC from daz_library
PRF=486.486
az_res=14 ##it can be improved extracting from par file
dfDCs=dl.get_dfDC(path_to_slcdir, f0=5405000500, burst_interval=2.758277, returnka=False, returnperswath=True)
print(dfDCs)

##rad2mm scaling factor.
scaling_factors = {1: (az_res*PRF) / (dfDCs[0]* 2 * np.pi), 2: (az_res*PRF) / (dfDCs[1]* 2 * np.pi), 3: (az_res*PRF) / (dfDCs[2]* 2 * np.pi)}
print(scaling_factors)

outbovl = bovlpha*0
for subswath in [1, 2, 3]:
    # Create a GeoDataFrame for the current subswath
    g = gpd.GeoDataFrame(
        {'bovl': globals()[f'overlap_gdf{subswath}'].index.values},
        geometry=globals()[f'overlap_gdf{subswath}'].geometry,
        crs={"init": "epsg:4326"}
    )

    # Create a GeoCube with the same spatial dimensions as 'aa'
    bovls = make_geocube(vector_data=g, like=aa.rio.set_spatial_dims(x_dim='lon', y_dim='lat'))
    bovls = bovls.rename({'x': 'lon', 'y': 'lat'})

    # Interpolate the 'bovls' data to match the spatial dimensions of 'bovlpha'
    #bovls = bovls.rio.interp_like(bovlpha)
    bovls = bovls.interp_like(bovlpha)

    # Create a binary mask where 'bovls' is multiplied by 0 and then added by 1
    bovls = bovls * 0 + 1

    # Multiply 'bovlpha' by the binary mask 'bovls' to apply the mask
    bovlphatemp = bovlpha * bovls
    if subswath in scaling_factors:
        bovlphatemp = bovlphatemp * scaling_factors[subswath]

    # add the grid values to the final output
    outbovl = outbovl + bovlphatemp
    # Export 'bovlphatemp' to a GeoTIFF file for the current subswath
    # export_xr2tif(bovlphatemp.bovl, f'subswath{subswath}.tif')

##I need the merge the subswath{sw}.tif, but I didn't manage it.    
export_xr2tif(outbovl.bovl, outtif) #.bovl, f'subswath{subswath}.tif')  
''' ML: MN, please test/check this line, I write without possibility to test it now - maybe should be outbovl.bovl?
    MN: I checked both way, the code work properly without any error but the output tiff with scale range between -3.40282e+38 and 3.40282e+38 doesn't seem merged subswath.
'''

'''
#MN's attempt to mosaic the subswaths through rasterio. But the code doesn't work properly. It produces tiff only including the first subswath.
## I put the code, maybe it can help me to find another way.
##merge the subswaths
subswath_files = ["subswath1.tif", "subswath2.tif", "subswath3.tif"]
src_files_to_mosaic = []
print(subswath_files)

for subswath_file in subswath_files:
    src = rasterio.open(subswath_file)
    src_files_to_mosaic.append(src)

print(src_files_to_mosaic)

mosaic, out_trans = merge(src_files_to_mosaic)

# Optionally, set the output file's spatial reference and affine transformation
out_meta = src.meta.copy()
out_meta.update({"driver": "GTiff",
                 "height": mosaic.shape[1],
                 "width": mosaic.shape[2],
                 "transform": out_trans})

# Specify the output file name and path
output_file = "merged_subswaths.tif"

# Create the merged GeoTIFF file
with rasterio.open(output_file, "w", **out_meta) as dest:
    dest.write(mosaic)
'''







