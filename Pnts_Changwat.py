#
#
#
import pandas as pd 
import geopandas as gpd 
import numpy as np
import pygeodesy as pgd
from shapely.geometry import Point, MultiPoint, Polygon
import matplotlib.pyplot as plt
import rasterio as rio
import argparse , sys
from pathlib import Path
gpd.options.io_engine = "pyogrio"

#DEM = r'/mnt/d/GeoData/NASADEM_HGT.001/NASADEM30m.vrt'
DEM = r'/mnt/d/GeoData/FABDEM_TH/FABDEM_Thailand.vrt'

ELLPS  = pgd.datums.Ellipsoids.WGS84

class PointPolygon:
    def __init__(self, ARGS, POLY ):
        RG = ELLPS.rocGauss( POLY.centroid.y )
        self.D2M = np.radians(RG)
        grid = self.GenGrid( POLY.bounds, ARGS.spacing )
        #import pdb ; pdb.set_trace()
        PntInside = POLY.intersection( MultiPoint( grid ) )
        pnt = gpd.GeoDataFrame( crs='EPSG:4326', geometry=[PntInside,] )
        self.dfPOINT = pnt.explode(index_parts=False, ignore_index=True )
        self.dfPOINT['Point']= self.dfPOINT.index.astype( str )
        #self.PlotGrid( gdf, Poly )
        self.GetMSL()
        Path('./CACHE').mkdir( parents=True, exist_ok=True )
        GPKG = f'CACHE/PointMSL {ARGS.prov}.gpkg'
        print(f'Wrting grip point to {GPKG}...' )
        self.dfPOINT.to_file( GPKG, driver='GPKG', layer='Point' )

    def GetMSL( self ):
        print( f'GetMSL() DEM={DEM}...' )
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
        #print( f'Before {len(self.dfPOINT) }' )
        #self.dfPOINT = self.dfPOINT[ (0<=self.dfPOINT.MSL) & (self.dfPOINT.MSL<=10) ]
        #print( f'After {len(self.dfPOINT) }' )
        #import pdb;pdb.set_trace()

    def GenGrid( self, BOUNDS, SPC ):
        [minx, miny, maxx, maxy] = BOUNDS 
        spc = SPC/self.D2M
        x = np.arange( minx-spc, maxx+spc, spc )
        y = np.arange( miny-spc, maxy+spc, spc )
        X,Y = np.meshgrid( x,y)
        return np.column_stack((X.ravel(), Y.ravel()))

    def PlotGrid(self, GRID, POLY ):
        if isinstance(GRID, gpd.GeoDataFrame):
            plt.scatter(GRID.geometry.x, GRID.geometry.y )
        else:
            plt.scatter(GRID[:, 0], GRID[:, 1])
        x, y = POLY.exterior.xy 
        plt.plot(x, y, 'r')  
        plt.show()

class PointProvince( PointPolygon ): 
    def __init__(self, ARGS):
        df = gpd.read_file( 'GADM/gadm41_THA.gpkg', layer='ADM_ADM_1' )
        if ARGS.list:
            print( list(df.NAME_1) ) , sys.exit(0)
        prov = ARGS.prov.split('|')
        if len(prov)==1:
            df = df[df.NAME_1==prov[0]].explode( index_parts=True )
        else:
            df = df[df.NAME_1.isin( prov )]
            if len(prov)!=len(df):
                print( f'Expecting ... {prov}' ); print( f'Found {list(df.NAME_1)}' )
                raise Exception( "provinces missing ..." )
            df = df.dissolve()
            df = df.simplify( 300/110_000 )
            df = gpd.GeoDataFrame( {'geometry': df} )
        super().__init__( ARGS, df.iloc[0].geometry )
        #print(f'Make grid with shape {grid.shape} and intersect with province(s) {ARGS.prov}, might long ...' )
        print(f'Make grid with shape and intersect with province(s) might be long time ...' )

class PointKML( PointPolygon ): 
    def __init__(self, ARGS):
        df = gpd.read_file(ARGS.kml, driver='KML' )
        assert len(df)==1
        assert isinstance( df.iloc[0].geometry, Polygon )
        coords = df.iloc[0].geometry.exterior.coords
        poly = [coord[:2] for coord in coords]
        super().__init__( ARGS, Polygon(poly) )

######################################################
parser = argparse.ArgumentParser( prog='Pnts_Changwat',
                    description='Generate points within province(s) specified by name(s)',
                    epilog='phisan.chula@gmail (Phisan Santitamnont), June 2024 ' )
parser.add_argument('-k', '--kml', help='specify kml... ' )
parser.add_argument('-p', '--prov', help='specify province name, if neighbour provinces use "|" ... ' )
parser.add_argument('-s', '--spacing',type=int, default=500, help='Spacing of point to generate, default 500 m' )    
parser.add_argument('-l', '--list', action='store_true', help='list all province names' )
args = parser.parse_args()
print( args )

#if args.prov==None and args.list==False:
#    parser.print_help()
#    exit(1)

if args.kml:
    pp = PointKML(args)
else:
    pp = PointProvince( args )
#pp = PointPolygon( args )

if 0:
    PROV1='Bangkok Metropolis'
    PROV2='Nonthaburi'
    PROV3='Samut Prakan'
    if 1:  # single province
        PROV = PROV2
    else:  # group of provinces
        PROV = f'{PROV1}|{PROV2}|{PROV3}'

    print( f'Create point sampling for {PROV} ...' )

print( '******************** finish *****************')
