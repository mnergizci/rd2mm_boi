import geopandas as gpd
import pandas as pd
import os

def extract_burst_overlaps(frame, jsonpath=os.getcwd()):
    # Read GeoJSON data
    data_temp = gpd.read_file(os.path.join(jsonpath, frame + '.geojson'))

    # Change CRS to EPSG:4326
    data_temp = data_temp.to_crs(epsg=4326)

    # Extract subswath information
    if frame.startswith('00'):
        data_temp['swath'] = data_temp.burstID.str[4]
    elif frame.startswith('0'):
        data_temp['swath'] = data_temp.burstID.str[5]
    else:
        data_temp['swath'] = data_temp.burstID.str[6]

    # Divide frame into subswaths
    data_temp = data_temp.sort_values(by=['burstID']).reset_index(drop=True)
    gpd_overlaps = None
    swathdict = dict()
    # ML: a fix to handle less than 3 swaths
    for swath in data_temp.swath.unique().values:
        swdata = data_temp[data_temp.swath == swath]
        # Divide burst overlaps into odd and even numbers
        a1 = swdata.iloc[::2]
        b1 = swdata.iloc[1::2]
        # Find burst overlaps
        sw_overlaps = gpd.overlay(a1, b1, how='intersection')
        swathdict[swath] = sw_overlaps
        if type(gpd_overlaps) == type(None):
            gpd_overlaps = sw_overlaps
        else:
            gpd_overlaps = pd.concat([gpd_overlaps, sw_overlaps], ignore_index=True)
    return gpd_overlaps, swathdict
'''
    sw1 = data_temp[data_temp.swath == '1']
    sw2 = data_temp[data_temp.swath == '2']
    sw3 = data_temp[data_temp.swath == '3']

    # Divide burst overlaps into odd and even numbers
    a1 = sw1.iloc[::2]
    b1 = sw1.iloc[1::2]
    a2 = sw2.iloc[::2]
    b2 = sw2.iloc[1::2]
    a3 = sw3.iloc[::2]
    b3 = sw3.iloc[1::2]

    # Find burst overlaps
    overlap_gdf1 = gpd.overlay(a1, b1, how='intersection')
    overlap_gdf2 = gpd.overlay(a2, b2, how='intersection')
    overlap_gdf3 = gpd.overlay(a3, b3, how='intersection')

    # Merge swath overlaps
    gpd_overlaps = gpd.GeoDataFrame(pd.concat([overlap_gdf1, overlap_gdf2, overlap_gdf3], ignore_index=True))

    return gpd_overlaps, overlap_gdf1, overlap_gdf2, overlap_gdf3
'''
