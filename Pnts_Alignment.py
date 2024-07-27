#
# readLandXML : read road alignment from designed CAD
#               then generate representation points over
#               the corridors
#
from pyproj import Proj, CRS
import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio as rio
import io
from pathlib import Path
import shapely
from shapely.geometry import Point,LineString
from shapely.ops import linemerge,snap
import matplotlib.pyplot as plt
import argparse , sys
gpd.options.io_engine = "pyogrio"

#DEM = r'/mnt/d/GeoData/NASADEM_HGT.001/NASADEM30m.vrt'
DEM = r'/mnt/d/GeoData/FABDEM_TH/FABDEM_Thailand.vrt'
#DEM = r'./Data/DEM_500kV_NR4-WN_WGS.tif'

class RoutePoints:
    def __init__(self):
        pass

    def PointsCorridor( self, DIV, ROW ):
        LASTDIV = None
        #import pdb; pdb.set_trace()
        SEG = self.ALIGN.geometry.to_list()
        LS = SEG[0] 
        for i, seg in enumerate( SEG[1:] ):
            LS = snap( LS,seg ,0.0001)  # coordinates transferred from LandXML have discontinuity
            LS = linemerge( [LS,seg] )
            assert LS.geom_type == 'LineString' 
        ROW_LR = ROW/2.
        nDiv,lastDiv = divmod( LS.length, DIV )
        dist = np.arange( 0, LS.length, DIV ) 
        Points = list()
        for i in range( len(dist) ):
            p0 = LS.interpolate( dist[i], normalized=False )  
            if i==len(dist)-1:
                p1 = LS.coords[-1]
                LASTDIV = True
            else:
                p1 = LS.interpolate( dist[i+1] , normalized=False )  
            div_buf = LineString( [p0,p1] ).buffer( ROW_LR, cap_style=2 )
            R = div_buf.exterior.coords[2]
            L = div_buf.exterior.coords[3]
            km,meter = divmod( dist[i],1_000 )
            STA = f'{int(km):03d}+{meter:03.0f}'
            Points.append( [ STA, p0 ] )
            Points.append( [ 'L', Point(L) ] )
            Points.append( [ 'R', Point(R) ] )
        if LASTDIV:
            km,meter = divmod( LS.length,1_000 )
            if meter>1: 
                R = div_buf.exterior.coords[1]
                L = div_buf.exterior.coords[4]
                STA = f'{int(km):03d}+{meter:03.0f}'
                Points.append( [ STA, Point(LS.coords[-1]) ] )
                Points.append( [ 'L', Point(L) ] )
                Points.append( [ 'R', Point(R) ] )
        df = pd.DataFrame( Points, columns=['Point', 'geometry' ] )
        self.dfPOINT = gpd.GeoDataFrame( df, crs='EPSG:32647', geometry=df.geometry )
        self.GetMSL()

    def GetMSL( self ):
        df = self.dfPOINT.to_crs('EPSG:4326').copy()
        with rio.open( DEM ) as dataset:
            metadata = dataset.meta
            print("Metadata:", metadata)
            # Read the dataset's data
            data = dataset.read(1)  # Read the first band
            print("Data shape:", data.shape)
            coords = list( zip( df.geometry.x, df.geometry.y ) )
            try:
                self.dfPOINT['MSL'] = np.rint( [ msl[0] for msl in dataset.sample( coords ) ] )
            except:
                print('***ERROR*** reading DEM...' )
            #import pdb;pdb.set_trace()

    def PlotMap(self, STEM ):
        Path('./CACHE').mkdir( parents=True, exist_ok=True )
        GPKG = f'CACHE/{STEM}.gpkg'
        print( f'Plotting GPKG {GPKG}... with layers Station|ROW|Alignment ... ' )
        self.ALIGN.to_file( GPKG , driver='GPKG', layer='Alignment' )
        ROW = self.dfPOINT.Point.str.contains('L|R')
        self.dfPOINT[ROW].to_file( GPKG, driver='GPKG', layer='ROW' )
        self.dfPOINT[~ROW ].to_file( GPKG, driver='GPKG', layer='Station' )
        #import pdb;pdb.set_trace()

class RouteKML( RoutePoints ):
    def __init__(self, KML ):
        align = gpd.read_file( KML ).iloc[0].geometry 
        if align.has_z:
            align = LineString([(x, y) for x, y, z in align.coords])
        self.ALIGN = gpd.GeoDataFrame( crs='EPSG:4326', geometry=[align,] )
        self.ALIGN = self.ALIGN.to_crs('EPSG:32647') 

class RouteLandXML( RoutePoints ):
    def __init__(self, landxml_file):
        self.ALIGN = self.ReadLandXML( landxml_file )
        self.CRS = self.ALIGN.crs.to_epsg()
        assert self.CRS == 32647
        print( self.ALIGN.AlignType.value_counts() )
        #self.CheckContinuity()

    def ReadLandXML(self, landxml ):
        ns = {'ns': 'http://www.landxml.org/schema/LandXML-1.2'}
        with open( landxml, 'r') as f:
            df = pd.read_xml( io.StringIO( f.read() ),  
                    xpath='.//ns:Alignment/ns:CoordGeom/*', namespaces=ns)
        def MakeLS(row):
            GeomType = { 2: 'Line', 3:'Spiral', 4:'Curve' }
            CurvPara = np.array( [row.Start, row.End, row.PI, row.Center ] )
            nStr = sum(isinstance(element, str) for element in CurvPara)
            ps = list( map(float, row.Start.split()  ) )[::-1]
            pe = list( map(float, row.End.split()   )  )[::-1]
            if nStr==2:
                LS = LineString( [ps,pe] )
            else:
                pi = list( map(float, row.End.split() ) )[::-1]
                LS = LineString( [ps,pi,pe] )
            return GeomType[nStr], LS
        df[['AlignType', 'geometry']] = df.apply( MakeLS, axis=1, result_type='expand' )
        gdf = gpd.GeoDataFrame( df , crs='EPSG:32647', geometry=df.geometry )
        return gdf

    def CheckContinuity(self):
        eps = list()
        for i in range( len(self.ALIGN)-1 ):
            ls = self.ALIGN.iloc[i:i+2].geometry
            d = ls.iloc[0].distance( ls.iloc[1] )
            eps.append( d )
        for EPS in [1E-7,1E-8,1E-9,1E-10,1E-11]:
            CHK=np.all(np.array(eps)<EPS )
            print(  f'Continuity epsilon {EPS:g} : {CHK} ' )

###################################################################################
###################################################################################
###################################################################################
if __name__=="__main__":
    print( f'DEM : {DEM}')
    parser = argparse.ArgumentParser( prog='Pnts_Alignment',
                        description='Generate points along an alignment defined by LandXML or KML',  
                        epilog='phisan.chula@gmail (Phisan Santitamnont), June 2024 ' )

    parser.add_argument('-d', '--div',type=int, default=500,
            help='division of stations along the alignment, layer="Station"' )
    parser.add_argument('-r', '--row',type=int, default=200,
            help='samplings on the left/right at the division, layer="ROW"' )
    parser.add_argument('-k', '--kml', help='specify alignment by kml' )
    parser.add_argument('-l', '--landxml', help='specify alignment by landxml' )
    args = parser.parse_args()
    print( args )
    #landxml_file = 'Data/Denchai-Chiangkhong.xml'
    if args.landxml:
        route = RouteLandXML( args.landxml )
        route.PointsCorridor( DIV=args.div, ROW=args.row )
        route.PlotMap( Path(args.landxml).stem )
    elif args.kml:
        #KML = Path(  'Data/Test01.kml' )
        #KML = Path(  'Data/Align_Srinakarin_95km.kml' )
        route = RouteKML( args.kml )
        route.PointsCorridor( DIV=args.div, ROW=args.row ) 
        route.PlotMap( Path(args.kml).stem )
    else:
        parser.print_help()
