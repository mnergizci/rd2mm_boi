import geopandas as gpd
import pandas as pd

def extract_burst_overlaps(frame):
    # Read GeoJSON data
    data_temp = gpd.read_file(frame + '.geojson')

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
