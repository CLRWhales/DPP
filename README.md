# DPP
Environment for DAS process paralelle and cleaning scripts 

Once installed in your environment there are a few actions you can take to simplyfy the processingof das data:

in the comand line:
Run ExampleInI to generate a template ini file for processing instructions
Run DASprocessParalelle to use the insturctions to process data quickly and save outputs
Run DAScleaner to visualize the cleaning type results and perform coars annotation

If you need to download AIS data for cleaning purposes it is possible to do so using the collection of functions available in 
Fetchmap AIS in a script similar to below: 

```python
import dpp

path = "" #path of file containing fiber data as UTMX,UTMY
zone = "33" #UTM zone
Southern = False #northern or southern hemisphere
Sp = "\t" #separator in the file
Bsize = 10000 #size of bounding box for AIS retrieval around fiber

username = "" #AIS username for Kystdatahuset
password = "" #AIS password
start = "" #time of start, eg "202208101300"
end = "" #time of end
minspeed = 0 #minimum travel speed



corners_latlon, linestring,crs_utm = dpp.cable2linestring(path,zone,Southern,Sp,Bsize)

bbox = ",".join(f"{x:.3f}" for x in corners_latlon)

wgs84_to_utm, utm_to_wgs84 = dpp.make_transformers(crs_utm)
token = dpp.get_token(username,password)
flag,ais_df = dpp.get_ais(bbox,start,end,minspeed,token)

mmsi = ais_df['mmsi'].unique().tolist()
flag,info = dpp.get_nameMMSI(mmsi,start,end,token)
mmsi_to_name = info.drop_duplicates('mmsi').set_index('mmsi')['name']
ais_df['name'] = dpp.ais_df['mmsi'].map(mmsi_to_name)
mappedAIS = dpp.project_ais_positions(ais_df,linestring,wgs84_to_utm,utm_to_wgs84,lat_col='lat',lon_col='lon')
mappedAIS.to_csv('file/name/here')
```

Note that the distance along fiber wil be relative to the 'start' ofthe fiber file, so depending on the order, you may need to reverse it.
It may also be worth while caching the AIS data frames locally after download so that you have them if they are deleted from the server.