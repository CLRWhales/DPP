#this script contains functions to fetch AIS data from Kystdatahuset, and map onto the fiber

import requests
import json
import pandas as pd
#import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point
#from geopy.distance import geodesic
import numpy as np
#import utm
#import geopandas as gpd
from pyproj import CRS, Transformer


#%% functions
def get_token(username,password):

    reqUrl = "https://kystdatahuset.no/ws/api/auth/login"

    headersList = {
        "User-Agent": "Your Client (https://your-client.com)",
        "accept": "*/*",
        "Content-Type": "application/json"
        }

    payload = json.dumps({
        "username": username,
        "password": password
        })

    response = requests.request("POST", reqUrl, data=payload,  headers=headersList).json()

    token = response['data']['JWT']
    return token

def get_ais(bbox,start,end,minspeed,token):


    aisURL = "https://kystdatahuset.no/ws/api/ais/positions/within-bbox-time"

    headersList = {
        "User-Agent": "Your Client (https://your-client.com)",
        "Content-Type": "application/json",
        "Authorization": token
    }

    payload = json.dumps({
        "bbox": bbox,
        "start": start,
        "end": end,
        "minSpeed": 0
    })

    response = requests.request("POST", aisURL, data=payload,  headers=headersList).json()
    flag = response['success']
    rows = response['data']
    columns = ["mmsi", "datetimeUTC", "lon", "lat", "CoG", "SoG", "MsgN", "S_kph", "Sec2prev", "Dist2prev", "TrueHeading","rate_of_turn"]

    df = pd.DataFrame(rows, columns=columns)
    
    return flag, df

def get_nameMMSI(mmsi, start,end,token):
    mmsiURL = "https://kystdatahuset.no/ws/api/ais/statinfo/for-mmsis"
    headersList = {
        "User-Agent": "Your Client (https://your-client.com)",
        "Content-Type": "application/json",
        "Authorization": token
    }
    payload = json.dumps({
        "mmsiIds": mmsi,
        "start": start,
        "end": end,
        "startTime": "0001-01-01T00:00:00"
    })
    response = requests.request("POST", mmsiURL, data=payload,  headers=headersList).json()
    flag = response['success']
    rows = response['data']
    df = pd.DataFrame(rows)
    df = df.rename(columns={
        "mmsino": "mmsi",
        "imono": "imo_num",
        "shipname": "name",
        "grosstonnage": "gt",
        "shiptypegroupnor": "type",
        "breadth": "depth"
    })

    if 'name' not in df.columns:
        df['name'] = df['mmsi']
    else:
        df['name'] = df['name'].fillna('') 
        df['name'] = df.apply(
            lambda row: f"Unknown {row['mmsi']}" if row['name'] == '' or pd.isna(row['name']) else row['name'],
            axis=1
        )

    df['name'] = df['name'].astype(str)

    return flag, df

def cable2linestring(path,zone,southern = False, sp = None,buffer_size = 10000):

    line_df = pd.read_csv(path,sep = sp)
    line_points = [Point(xy) for xy in zip(line_df["UTMX"], line_df["UTMY"])]
    line = LineString(line_points)
    buffered_box = line.buffer(buffer_size).bounds
    crs_utm = CRS.from_dict({"proj": "utm", "zone": zone, "datum": "WGS84","south": southern})
    crs_wgs = CRS.from_epsg(4326)
    transformer = Transformer.from_crs(crs_utm, crs_wgs, always_xy=True)
    minx, miny, maxx, maxy = buffered_box
    corners_utm = [(minx, miny), (maxx, maxy)]
    corners_latlon = [transformer.transform(x, y) for x, y in corners_utm]

    return corners_latlon, line, crs_utm

def make_transformers(crs_utm):
    """
    Returns (wgs84_to_utm, utm_to_wgs84) pyproj Transformer instances.
    """
    crs_wgs = CRS.from_epsg(4326)
    wgs84_to_utm = Transformer.from_crs(crs_wgs, crs_utm, always_xy=True)
    utm_to_wgs84 = Transformer.from_crs(crs_utm, crs_wgs, always_xy=True)
    return wgs84_to_utm, utm_to_wgs84

def project_ais_positions(ais_df, line_utm, wgs84_to_utm, utm_to_wgs84,
                          lat_col="lat", lon_col="lon"):
    """
    For each AIS position, compute:
      - along_track_m: distance along the cable (m)
      - cross_track_m: perpendicular distance (m, unsigned)
      - proj_lat/lon: coordinates of the projection point on the cable
    """
    #df = ais_df.copy()

    # Convert AIS positions to UTM
    x, y = wgs84_to_utm.transform(ais_df[lon_col].values, ais_df[lat_col].values)
    ais_df["UTMX"], ais_df["UTMY"] = x, y

    total_len = line_utm.length

    along_distances = []
    cross_distances = []
    proj_lats = []
    proj_lons = []

    for xi, yi in zip(x, y):
        pt = Point(xi, yi)
        d_along = line_utm.project(pt)
        closest_pt = line_utm.interpolate(d_along)
        d_cross = pt.distance(closest_pt)

        proj_lon, proj_lat = utm_to_wgs84.transform(closest_pt.x, closest_pt.y)

        along_distances.append(d_along)
        cross_distances.append(d_cross)
        proj_lats.append(proj_lat)
        proj_lons.append(proj_lon)


    df = pd.DataFrame({
        "timestamp":ais_df["datetimeUTC"],
        "name":ais_df["name"],
        "along_track_m":along_distances,
        "cross_track_m":cross_distances
    })
    # df["along_track_m"] = np.array(along_distances)
    # df["cross_track_m"] = np.array(cross_distances)
    # df["along_fraction"] = df["along_track_m"] / total_len
    # df["proj_lat"] = np.array(proj_lats)
    # df["proj_lon"] = np.array(proj_lons)

    return df


# def create_ais():
#     #we need to output: timestamp, name, proportion along fiber from start
#     path = ""
#     zone = "33"
#     Southern = False
#     Sp = "\t"
#     Bsize = 10000

#     username = ""
#     password = ""
#     start = ""
#     end = ""
#     minspeed = 0
#     rev = False


#     corners_latlon, linestring,crs_utm = cable2linestring(path,zone,Southern,Sp,Bsize)

#     bbox = ",".join(f"{x:.3f}" for x in corners_latlon)

#     wgs84_to_utm, utm_to_wgs84 = make_transformers(crs_utm)
#     token = get_token(username,password)
#     flag,ais_df = get_ais(bbox,start,end,minspeed,token)
#     if not flag:
#         print("error something went wrong, do you have the right password?")
#     mmsi = ais_df['mmsi'].unique().tolist()
#     flag,info = get_nameMMSI(mmsi,start,end,token)
#     mmsi_to_name = info.drop_duplicates('mmsi').set_index('mmsi')['name']
#     ais_df['name'] = ais_df['mmsi'].map(mmsi_to_name)
#     result_inner = project_ais_positions(ais_df,linestring,wgs84_to_utm,utm_to_wgs84,lat_col='lat',lon_col='lon')
    
