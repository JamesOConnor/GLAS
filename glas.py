### Use input HDF5 files to output a 3D model as seen by GLAS on a given day

import h5py
import numpy as np
import numpy.ma as ma
import glob

def LLHtoECEF(lat, lon, alt):
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    lat=np.deg2rad(lat)
    lon=np.deg2rad(lon)
    rad = np.float64(6378137.0)        # Radius of the Earth (in meters)
    f = np.float64(1.0/298.257223563)  # Flattening factor WGS84 Model
    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    FF     = (1.0-f)**2
    C      = 1/np.sqrt(cosLat**2 + FF * sinLat**2)
    S      = C * FF

    x = (rad * C + alt)*cosLat * np.cos(lon)
    y = (rad * C + alt)*cosLat * np.sin(lon)
    z = (rad * S + alt)*sinLat

    return (x, y, z)

def process(path,outpath):
	fns=glob.glob('%s/*.H5'%(path))
	for n,fn in enumerate(fns):
		b=fn
		indata=h5py.File(b)
		returns=np.array(indata['Data_40HZ']['Elevation_Surfaces']['d_elev'])
		lats=np.array(indata['Data_40HZ']['Geolocation']['d_lat'])
		lons=np.array(indata['Data_40HZ']['Geolocation']['d_lon'])
		RLL=np.zeros((lats.size,3))
		RLL[:,0]=lats
		RLL[:,1]=lons
		RLL[:,2]=returns
		RLL_ma=ma.masked_greater(RLL, 100000)
		RLL_onlyvalid=RLL_ma[~RLL_ma.mask]
		RLL_final=np.zeros((RLL_onlyvalid.size/3,3))
		RLL_final[:,0]=RLL_onlyvalid[::3]
		RLL_final[:,1]=RLL_onlyvalid[1::3]
		RLL_final[:,2]=RLL_onlyvalid[2::3]
		cart=np.zeros_like(RLL_final)
		for i in enumerate(RLL_final[:,0]):
			cart[i,0],cart[i,1],cart[i,2]=LLHtoECEF(RLL_final[i,0],RLL_final[i,1],RLL_final[i,2])
		cart2=np.zeros_like(cart)
		for i in enumerate(RLL_final[:,0]):
			cart2[i,0],cart2[i,1],cart2[i,2]=LLHtoECEF(RLL_final[i,0],RLL_final[i,1],0)
		np.savetxt('%s\\%s_XYZ.txt'%(outpath,n), cart)
		np.savetxt('%s\\%s_XYZ_ref.txt'%(outpath,n), cart2)
