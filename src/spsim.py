# SPSIM PYTHON MODULE

__version__ = "0.1.0"

import numpy

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
from io import StringIO
from io import BytesIO
import gzip
import Bio.PDB.PDBParser
import configparser
import matplotlib.pyplot as plt
import h5py
try:
    import jax.numpy as jnp
    import jax
except ImportError:
    print('JAX not available.')

try:
    import cupy
except ImportError:
    print('CuPy not available.')

atomsf = numpy.array([
[-1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00]
,[4.930020e-01,3.229120e-01,1.401910e-01,4.081000e-02,1.051090e+01,2.612570e+01,3.142360e+00,5.779970e+01,3.038000e-03]
,[8.734000e-01,6.309000e-01,3.112000e-01,1.780000e-01,9.103700e+00,3.356800e+00,2.292760e+01,9.821000e-01,6.400000e-03]
,[1.128200e+00,7.508000e-01,6.175000e-01,4.653000e-01,3.954600e+00,1.052400e+00,8.539050e+01,1.682610e+02,3.770000e-02]
,[1.591900e+00,1.127800e+00,5.391000e-01,7.029000e-01,4.364270e+01,1.862300e+00,1.034830e+02,5.420000e-01,3.850000e-02]
,[2.054500e+00,1.332600e+00,1.097900e+00,7.068000e-01,2.321850e+01,1.021000e+00,6.034980e+01,1.403000e-01,-1.932000e-01]
,[2.310000e+00,1.020000e+00,1.588600e+00,8.650000e-01,2.084390e+01,1.020750e+01,5.687000e-01,5.165120e+01,2.156000e-01]
,[1.221260e+01,3.132200e+00,2.012500e+00,1.166300e+00,5.700000e-03,9.893300e+00,2.899750e+01,5.826000e-01,-1.152900e+01]
,[3.048500e+00,2.286800e+00,1.546300e+00,8.670000e-01,1.327710e+01,5.701100e+00,3.239000e-01,3.290890e+01,2.508000e-01]
,[3.539200e+00,2.641200e+00,1.517000e+00,1.024300e+00,1.028250e+01,4.294400e+00,2.615000e-01,2.614760e+01,2.776000e-01]
,[3.955300e+00,3.112500e+00,1.454600e+00,1.125100e+00,8.404200e+00,3.426200e+00,2.306000e-01,2.171840e+01,3.515000e-01]
,[4.762600e+00,3.173600e+00,1.267400e+00,1.112800e+00,3.285000e+00,8.842199e+00,3.136000e-01,1.294240e+02,6.760000e-01]
,[5.420400e+00,2.173500e+00,1.226900e+00,2.307300e+00,2.827500e+00,7.926110e+01,3.808000e-01,7.193700e+00,8.584000e-01]
,[6.420200e+00,1.900200e+00,1.593600e+00,1.964600e+00,3.038700e+00,7.426000e-01,3.154720e+01,8.508860e+01,1.115100e+00]
,[6.291500e+00,3.035300e+00,1.989100e+00,1.541000e+00,2.438600e+00,3.233370e+01,6.785000e-01,8.169370e+01,1.140700e+00]
,[6.434500e+00,4.179100e+00,1.780000e+00,1.490800e+00,1.906700e+00,2.715700e+01,5.260000e-01,6.816450e+01,1.114900e+00]
,[6.905300e+00,5.203400e+00,1.437900e+00,1.586300e+00,1.467900e+00,2.221510e+01,2.536000e-01,5.617200e+01,8.669000e-01]
,[1.146040e+01,7.196400e+00,6.255600e+00,1.645500e+00,1.040000e-02,1.166200e+00,1.851940e+01,4.777840e+01,-9.557400e+00]
,[7.484500e+00,6.772300e+00,6.539000e-01,1.644200e+00,9.072000e-01,1.484070e+01,4.389830e+01,3.339290e+01,1.444500e+00]
,[8.218599e+00,7.439800e+00,1.051900e+00,8.659000e-01,1.279490e+01,7.748000e-01,2.131870e+02,4.168410e+01,1.422800e+00]
,[8.626600e+00,7.387300e+00,1.589900e+00,1.021100e+00,1.044210e+01,6.599000e-01,8.574840e+01,1.784370e+02,1.375100e+00]
,[9.189000e+00,7.367900e+00,1.640900e+00,1.468000e+00,9.021299e+00,5.729000e-01,1.361080e+02,5.135310e+01,1.332900e+00]
,[9.759500e+00,7.355800e+00,1.699100e+00,1.902100e+00,7.850800e+00,5.000000e-01,3.563380e+01,1.161050e+02,1.280700e+00]
,[1.029710e+01,7.351100e+00,2.070300e+00,2.057100e+00,6.865700e+00,4.385000e-01,2.689380e+01,1.024780e+02,1.219900e+00]
,[1.064060e+01,7.353700e+00,3.324000e+00,1.492200e+00,6.103800e+00,3.920000e-01,2.026260e+01,9.873990e+01,1.183200e+00]
,[1.128190e+01,7.357300e+00,3.019300e+00,2.244100e+00,5.340900e+00,3.432000e-01,1.786740e+01,8.375430e+01,1.089600e+00]
,[1.176950e+01,7.357300e+00,3.522200e+00,2.304500e+00,4.761100e+00,3.072000e-01,1.535350e+01,7.688050e+01,1.036900e+00]
,[1.228410e+01,7.340900e+00,4.003400e+00,2.348800e+00,4.279100e+00,2.784000e-01,1.353590e+01,7.116920e+01,1.011800e+00]
,[1.283760e+01,7.292000e+00,4.443800e+00,2.380000e+00,3.878500e+00,2.565000e-01,1.217630e+01,6.634210e+01,1.034100e+00]
,[1.333800e+01,7.167600e+00,5.615800e+00,1.673500e+00,3.582800e+00,2.470000e-01,1.139660e+01,6.481260e+01,1.191000e+00]
,[1.407430e+01,7.031800e+00,5.162500e+00,2.410000e+00,3.265500e+00,2.333000e-01,1.031630e+01,5.870970e+01,1.304100e+00]
,[1.523540e+01,6.700600e+00,4.359100e+00,2.962300e+00,3.066900e+00,2.412000e-01,1.078050e+01,6.141350e+01,1.718900e+00]
,[1.608160e+01,6.374700e+00,3.706800e+00,3.683000e+00,2.850900e+00,2.516000e-01,1.144680e+01,5.476250e+01,2.131300e+00]
,[1.667230e+01,6.070100e+00,3.431300e+00,4.277900e+00,2.634500e+00,2.647000e-01,1.294790e+01,4.779720e+01,2.531000e+00]
,[1.700060e+01,5.819600e+00,3.973100e+00,4.354300e+00,2.409800e+00,2.726000e-01,1.523720e+01,4.381630e+01,2.840900e+00]
,[1.717890e+01,5.235800e+00,5.637700e+00,3.985100e+00,2.172300e+00,1.657960e+01,2.609000e-01,4.143280e+01,2.955700e+00]
,[1.735550e+01,6.728600e+00,5.549300e+00,3.537500e+00,1.938400e+00,1.656230e+01,2.261000e-01,3.939720e+01,2.825000e+00]
,[1.717840e+01,9.643499e+00,5.139900e+00,1.529200e+00,1.788800e+00,1.731510e+01,2.748000e-01,1.649340e+02,3.487300e+00]
,[1.756630e+01,9.818399e+00,5.422000e+00,2.669400e+00,1.556400e+00,1.409880e+01,1.664000e-01,1.323760e+02,2.506400e+00]
,[1.777600e+01,1.029460e+01,5.726290e+00,3.265880e+00,1.402900e+00,1.280060e+01,1.255990e-01,1.043540e+02,1.912130e+00]
,[1.787650e+01,1.094800e+01,5.417320e+00,3.657210e+00,1.276180e+00,1.191600e+01,1.176220e-01,8.766270e+01,2.069290e+00]
,[1.761420e+01,1.201440e+01,4.041830e+00,3.533460e+00,1.188650e+00,1.176600e+01,2.047850e-01,6.979570e+01,3.755910e+00]
,[3.702500e+00,1.723560e+01,1.288760e+01,3.742900e+00,2.772000e-01,1.095800e+00,1.100400e+01,6.165840e+01,4.387500e+00]
,[1.913010e+01,1.109480e+01,4.649010e+00,2.712630e+00,8.641320e-01,8.144870e+00,2.157070e+01,8.684720e+01,5.404280e+00]
,[1.926740e+01,1.291820e+01,4.863370e+00,1.567560e+00,8.085200e-01,8.434669e+00,2.479970e+01,9.429280e+01,5.378740e+00]
,[1.929570e+01,1.435010e+01,4.734250e+00,1.289180e+00,7.515360e-01,8.217580e+00,2.587490e+01,9.860620e+01,5.328000e+00]
,[1.933190e+01,1.550170e+01,5.295370e+00,6.058440e-01,6.986550e-01,7.989290e+00,2.520520e+01,7.689860e+01,5.265930e+00]
,[1.928080e+01,1.668850e+01,4.804500e+00,1.046300e+00,6.446000e-01,7.472600e+00,2.466050e+01,9.981560e+01,5.179000e+00]
,[1.922140e+01,1.764440e+01,4.461000e+00,1.602900e+00,5.946000e-01,6.908900e+00,2.470080e+01,8.748250e+01,5.069400e+00]
,[1.916240e+01,1.855960e+01,4.294800e+00,2.039600e+00,5.476000e-01,6.377600e+00,2.584990e+01,9.280290e+01,4.939100e+00]
,[1.918890e+01,1.910050e+01,4.458500e+00,2.466300e+00,5.830300e+00,5.031000e-01,2.689090e+01,8.395710e+01,4.782100e+00]
,[1.964180e+01,1.904550e+01,5.037100e+00,2.682700e+00,5.303400e+00,4.607000e-01,2.790740e+01,7.528250e+01,4.590900e+00]
,[1.996440e+01,1.901380e+01,6.144870e+00,2.523900e+00,4.817420e+00,4.208850e-01,2.852840e+01,7.084030e+01,4.352000e+00]
,[2.014720e+01,1.899490e+01,7.513800e+00,2.273500e+00,4.347000e+00,3.814000e-01,2.776600e+01,6.687760e+01,4.071200e+00]
,[2.029330e+01,1.902980e+01,8.976700e+00,1.990000e+00,3.928200e+00,3.440000e-01,2.646590e+01,6.426580e+01,3.711800e+00]
,[2.038920e+01,1.910620e+01,1.066200e+01,1.495300e+00,3.569000e+00,3.107000e-01,2.438790e+01,2.139040e+02,3.335200e+00]
,[2.033610e+01,1.929700e+01,1.088800e+01,2.695900e+00,3.216000e+00,2.756000e-01,2.020730e+01,1.672020e+02,2.773100e+00]
,[2.057800e+01,1.959900e+01,1.137270e+01,3.287190e+00,2.948170e+00,2.444750e-01,1.877260e+01,1.331240e+02,2.146780e+00]
,[2.116710e+01,1.976950e+01,1.185130e+01,3.330490e+00,2.812190e+00,2.268360e-01,1.760830e+01,1.271130e+02,1.862640e+00]
,[2.204400e+01,1.966970e+01,1.238560e+01,2.824280e+00,2.773930e+00,2.220870e-01,1.676690e+01,1.436440e+02,2.058300e+00]
,[2.268450e+01,1.968470e+01,1.277400e+01,2.851370e+00,2.662480e+00,2.106280e-01,1.588500e+01,1.379030e+02,1.984860e+00]
,[2.334050e+01,1.960950e+01,1.312350e+01,2.875160e+00,2.562700e+00,2.020880e-01,1.510090e+01,1.327210e+02,2.028760e+00]
,[2.400420e+01,1.942580e+01,1.343960e+01,2.896040e+00,2.472740e+00,1.964510e-01,1.439960e+01,1.280070e+02,2.209630e+00]
,[2.462740e+01,1.908860e+01,1.376030e+01,2.922700e+00,2.387900e+00,1.942000e-01,1.375460e+01,1.231740e+02,2.574500e+00]
,[2.507090e+01,1.907980e+01,1.385180e+01,3.545450e+00,2.253410e+00,1.819510e-01,1.293310e+01,1.013980e+02,2.419600e+00]
,[2.589760e+01,1.821850e+01,1.431670e+01,2.953540e+00,2.242560e+00,1.961430e-01,1.266480e+01,1.153620e+02,3.582240e+00]
,[2.650700e+01,1.763830e+01,1.455960e+01,2.965770e+00,2.180200e+00,2.021720e-01,1.218990e+01,1.118740e+02,4.297280e+00]
,[2.690490e+01,1.729400e+01,1.455830e+01,3.638370e+00,2.070510e+00,1.979400e-01,1.144070e+01,9.265660e+01,4.567960e+00]
,[2.765630e+01,1.642850e+01,1.497790e+01,2.982330e+00,2.073560e+00,2.235450e-01,1.136040e+01,1.057030e+02,5.920460e+00]
,[2.818190e+01,1.588510e+01,1.515420e+01,2.987060e+00,2.028590e+00,2.388490e-01,1.099750e+01,1.029610e+02,6.756210e+00]
,[2.866410e+01,1.543450e+01,1.530870e+01,2.989630e+00,1.988900e+00,2.571190e-01,1.066470e+01,1.004170e+02,7.566720e+00]
,[2.894760e+01,1.522080e+01,1.510000e+01,3.716010e+00,1.901820e+00,9.985189e+00,2.610330e-01,8.432980e+01,7.976280e+00]
,[2.914400e+01,1.517260e+01,1.475860e+01,4.300130e+00,1.832620e+00,9.599899e+00,2.751160e-01,7.202900e+01,8.581540e+00]
,[2.920240e+01,1.522930e+01,1.451350e+01,4.764920e+00,1.773330e+00,9.370460e+00,2.959770e-01,6.336440e+01,9.243540e+00]
,[2.908180e+01,1.543000e+01,1.443270e+01,5.119820e+00,1.720290e+00,9.225900e+00,3.217030e-01,5.705600e+01,9.887500e+00]
,[2.876210e+01,1.571890e+01,1.455640e+01,5.441740e+00,1.671910e+00,9.092270e+00,3.505000e-01,5.208610e+01,1.047200e+01]
,[2.818940e+01,1.615500e+01,1.493050e+01,5.675890e+00,1.629030e+00,8.979480e+00,3.826610e-01,4.816470e+01,1.100050e+01]
,[2.730490e+01,1.672960e+01,1.561150e+01,5.833770e+00,1.592790e+00,8.865530e+00,4.179160e-01,4.500110e+01,1.147220e+01]
,[2.700590e+01,1.776390e+01,1.571310e+01,5.783700e+00,1.512930e+00,8.811740e+00,4.245930e-01,3.861030e+01,1.168830e+01]
,[1.688190e+01,1.859130e+01,2.555820e+01,5.860000e+00,4.611000e-01,8.621600e+00,1.482600e+00,3.639560e+01,1.206580e+01]
,[2.068090e+01,1.904170e+01,2.165750e+01,5.967600e+00,5.450000e-01,8.448400e+00,1.572900e+00,3.832460e+01,1.260890e+01]
,[2.754460e+01,1.915840e+01,1.553800e+01,5.525930e+00,6.551500e-01,8.707510e+00,1.963470e+00,4.581490e+01,1.317460e+01]
,[3.106170e+01,1.306370e+01,1.844200e+01,5.969600e+00,6.902000e-01,2.357600e+00,8.618000e+00,4.725790e+01,1.341180e+01]
,[3.336890e+01,1.295100e+01,1.658770e+01,6.469200e+00,7.040000e-01,2.923800e+00,8.793700e+00,4.800930e+01,1.357820e+01]
,[3.467260e+01,1.547330e+01,1.311380e+01,7.025880e+00,7.009990e-01,3.550780e+00,9.556419e+00,4.700450e+01,1.367700e+01]
,[3.531630e+01,1.902110e+01,9.498870e+00,7.425180e+00,6.858700e-01,3.974580e+00,1.138240e+01,4.547150e+01,1.371080e+01]
,[3.556310e+01,2.128160e+01,8.003700e+00,7.443300e+00,6.631000e-01,4.069100e+00,1.404220e+01,4.424730e+01,1.369050e+01]
,[3.592990e+01,2.305470e+01,1.214390e+01,2.112530e+00,6.464530e-01,4.176190e+00,2.310520e+01,1.506450e+02,1.372470e+01]
,[3.576300e+01,2.290640e+01,1.247390e+01,3.210970e+00,6.163410e-01,3.871350e+00,1.998870e+01,1.423250e+02,1.362110e+01]
,[3.565970e+01,2.310320e+01,1.259770e+01,4.086550e+00,5.890920e-01,3.651550e+00,1.859900e+01,1.170200e+02,1.352660e+01]
,[3.556450e+01,2.342190e+01,1.274730e+01,4.807030e+00,5.633590e-01,3.462040e+00,1.783090e+01,9.917220e+01,1.343140e+01]
,[3.588470e+01,2.329480e+01,1.418910e+01,4.172870e+00,5.477510e-01,3.415190e+00,1.692350e+01,1.052510e+02,1.342870e+01]
,[3.602280e+01,2.341280e+01,1.494910e+01,4.188000e+00,5.293000e-01,3.325300e+00,1.609270e+01,1.006130e+02,1.339660e+01]
,[3.618740e+01,2.359640e+01,1.564020e+01,4.185500e+00,5.119290e-01,3.253960e+00,1.536220e+01,9.749080e+01,1.335730e+01]
,[3.652540e+01,2.380830e+01,1.677070e+01,3.479470e+00,4.993840e-01,3.263710e+00,1.494550e+01,1.059800e+02,1.338120e+01]
,[3.667060e+01,2.409920e+01,1.734150e+01,3.493310e+00,4.836290e-01,3.206470e+00,1.431360e+01,1.022730e+02,1.335920e+01]
,[3.664880e+01,2.440960e+01,1.739900e+01,4.216650e+00,4.651540e-01,3.089970e+00,1.343460e+01,8.848340e+01,1.328870e+01]
,[3.678810e+01,2.477360e+01,1.789190e+01,4.232840e+00,4.510180e-01,3.046190e+00,1.289460e+01,8.600300e+01,1.327540e+01]
,[3.691850e+01,2.519950e+01,1.833170e+01,4.243910e+00,4.375330e-01,3.007750e+00,1.240440e+01,8.378810e+01,1.326740e+01]
,[-1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00]
])


element_to_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'NA': 11, 'MG': 12, 'P': 15, 'S': 16,  'CL': 17}


class Detector:
    def __init__(self, config) -> None:
        self.center_x = float(config['detector_center_x'])
        self.center_y = float(config['detector_center_y'])
        self.center_z = float(config['detector_center_z'])
        self.distance = float(config['detector_distance'])
        self.width = float(config['detector_width'])
        self.height = float(config['detector_height'])
        self.depth = float(config['detector_depth'])
        self.pixel_width = float(config['detector_pixel_width'])
        self.pixel_height = float(config['detector_pixel_height'])
        self.pixel_depth = float(config['detector_pixel_depth'])
        self.quantum_efficiency = float(config['detector_quantum_efficiency'])
        self.electron_hole_production_energy = float(config['detector_electron_hole_production_energy'])
        self.nx = round(self.width/self.pixel_width)
        self.ny = round(self.height/self.pixel_height)
        if self.pixel_depth:
            self.nz = round(self.depth/self.pixel_depth)
        else:
            self.nz = 1
            self.binning_z = 1
        # Calculate the coordinates of the center of the pixel
        # The + 0.5 accoutn for the center
        # The -0.5 is to center the detector
        self.px = ((numpy.arange(self.nx) + 0.5)/self.nx - 0.5) * self.width - self.center_x
        self.px = self.px[numpy.newaxis, :]
        self.py = ((numpy.arange(self.ny) + 0.5)/self.ny - 0.5) * self.height - self.center_y
        self.py = self.py[:, numpy.newaxis]
        self.pz = self.distance
        # Distance from the center of the pixel to the origin
        self.pr = numpy.sqrt(self.px**2 + self.py**2 +  self.pz**2)
        self.thomson_correction = None
        self.solid_angle = None
            
class Experiment:
    def __init__(self, config) -> None:
        if('experiment_wavelength' in config):            
            self.wavelength = float(config['experiment_wavelength'])
            self.photon_energy = 1.240e-6/self.wavelength # in eV
        elif('experiment_photon_energy' in config):
            self.photon_energy = float(config['experiment_photon_energy'])
            self.wavelength = 1.240e-6/self.photon_energy # in m
        else:
            raise ValueError('You need to specify the experiment_wavelength in the config file')

        self.beam_intensity = float(config['experiment_beam_intensity'])
        self.polarization = config['experiment_polarization']

        pass
        
class Options:
    def __init__(self, config):
        # Remove all the  ';' from the options
        for k in config.keys():
            config[k] = config[k].replace(';','')

        self.detector = Detector(config)
        self.experiment = Experiment(config)
        self.sf_filename = None
        self.random_orientation = int(config['random_orientation'])
        self.use_fft_for_sf = int(config['use_fft_for_sf'])
        self.use_nfft_for_sf = int(config['use_nfft_for_sf'])
        self.use_cuda = int(config['use_cuda'])
        self.use_jax = int(config['use_jax'])
        self.euler_orientation = numpy.zeros((3))
        self.euler_orientation[0] = float(config['phi'])
        self.euler_orientation[1] = float(config['theta'])
        self.euler_orientation[2] = float(config['psi'])

        self.b_factor = float(config['b_factor'])
        if 'pdb_filename' in config:
            self.pdb_filename = config['pdb_filename']
        else:
            self.pdb_filename = None

    
class Sample:
    def __init__(self, pos, Z, B=None):
        self.pos = pos
        self.Z = Z
        self.B = B
        if len(pos) != len(Z):
            raise ValueError
        self.natoms = len(pos)

class Pattern:
    def __init__(self, F, HKL_list):
        self.F = F
        self.HKL_list = HKL_list


def get_Sample_from_pdb(filename):
    parser = Bio.PDB.PDBParser(PERMISSIVE=True)
    structure = parser.get_structure('Sample', filename)
    atoms = list(structure.get_atoms())
    pos = numpy.zeros((len(atoms),3),dtype=numpy.float32)
    Z = numpy.zeros((len(atoms)),dtype=int)
    B = numpy.zeros((len(atoms)))
    for i,atom in enumerate(atoms):
        pos[i] = atom.get_coord()
        element = atom.element
        B[i] = atom.bfactor
        if element in element_to_Z:            
            Z[i] = element_to_Z[element]
        else:
            raise KeyError('%s not in element_to_Z dictionary' % element)
    # Convert to meters
    print('Read %d atoms with %d electrons' % (len(atoms), Z.sum()))
    return Sample(pos*1e-10, Z, B)


def load_pattern_from_file(det, filename):
    raise NotImplementedError

def get_HKL_list(opts):    
    if(opts.detector.nz > 1):
        raise NotImplementedError
    else:
        HKL_list = get_HKL_list_for_detector(opts.detector,opts.experiment)
    return HKL_list;

def get_HKL_list_for_detector(det, exp):
    ewald_radius = 1.0/exp.wavelength
    rx = det.px/det.pr
    ry = det.py/det.pr
    rz = det.pz/det.pr-1
    HKL_list = numpy.zeros((rx.shape[0]*rx.shape[1],3),dtype=numpy.float32)
    HKL_list[:,0] = rx.flatten() * ewald_radius
    HKL_list[:,1] = ry.flatten() * ewald_radius
    HKL_list[:,2] = rz.flatten() * ewald_radius
    return HKL_list

def apply_orientation_to_HKL_list(HKL_list, opts):
    import scipy.spatial.transform
    if(opts.random_orientation):
        rot = scipy.spatial.transform.Rotation.random()
    else:
        rot = scipy.spatial.transform.Rotation.from_euler('zxz', opts.euler_orientation, degrees=True)
    HKL_list[:,:] = rot.apply(HKL_list)
    return rot

def compute_sf(sample, HKL_list, opts):
    try:
        import cupy
        cupy_available = True
    except:
        cupy_available = False
        pass
    is3D = (opts.detector.nz > 1)
    if(opts.use_fft_for_sf):
        if(is3D is False):
            raise ValueError('Cannot use fft for 2D pattern calculation!')
        raise NotImplemented
    elif(opts.use_nfft_for_sf):
        raise NotImplementedError
    else:
        if(jax and opts.use_jax):
            pattern = jax_compute_pattern_on_list(sample, HKL_list, opts)
        elif(opts.use_cuda and cupy_available):
            pattern = cuda_compute_pattern_on_list(sample, HKL_list, opts)
        else:
            pattern = compute_pattern_on_list(sample, HKL_list, opts)
    pattern.rot = None
    return pattern    

def cuda_compute_pattern_on_list(sample, HKL_list, opts):
    F = cupy.zeros(HKL_list.shape[0],dtype=cupy.complex128)
    scattering_vector_length = cupy.sqrt(cupy.sum(HKL_list**2, axis=1))
    for i in range(sample.natoms):
        sf = scatt_factor(scattering_vector_length, sample.Z[i], opts.b_factor)
        phase = -2 * cupy.pi * cupy.dot(HKL_list,sample.pos[i])
        F += sf*cupy.exp(1.0j*phase)

    pat = Pattern(F, HKL_list)
    pat.ints = cupy.abs(F)**2   
    return pat


def compute_pattern_on_list(sample, HKL_list, opts):
    F = numpy.zeros(HKL_list.shape[0],dtype=numpy.complex64)
    scattering_vector_length = numpy.sqrt(numpy.sum(HKL_list**2, axis=1))
    points_per_percent = numpy.ceil(sample.natoms/100)
    for i in range(sample.natoms):
        if(i % points_per_percent == 0):
            print('%3.1f percent done' % (i/sample.natoms*100))
        #sf = scatt_factor(scattering_vector_length, sample.Z[i], opts.b_factor)
        sf = scatt_factor(scattering_vector_length, sample.Z[i], sample.B[i])
        phase = -2 * numpy.pi * numpy.dot(HKL_list,sample.pos[i])
        F += sf*(numpy.cos(phase)+1.0j*numpy.sin(phase))

    pat = Pattern(F, HKL_list)
    pat.ints = numpy.abs(F)**2   
    return pat

def jax_compute_pattern_on_batch(pos, Z, B, HKL, Z_cache, Z_map, d2):
    # Caching B was not worth it
    phase = -2 * jnp.pi * jnp.inner(pos, HKL)
    sf = Z_cache[Z_map[Z]]*jnp.exp(-B[:,None]*d2)
    return jnp.sum(sf*(jnp.cos(phase)), axis=0), jnp.sum(sf*(jnp.sin(phase)), axis=0)

def jax_compute_pattern_on_list(sample, HKL_list, opts):
    print('Using jax')
    print('npixels: %d' % (HKL_list.shape[0]))
    batch_size = int(min(1e8//HKL_list.shape[0], sample.natoms))
    print('batch_size: %d' % (batch_size))
    # JAX on macos has no support for complex64 so we'll use two float32
    F_r = jnp.zeros(HKL_list.shape[0],dtype=jnp.float32)
    F_i = jnp.zeros(HKL_list.shape[0],dtype=jnp.float32)

    scattering_vector_length = jnp.array(numpy.sqrt(numpy.sum(HKL_list**2, axis=1)))
    # Convert to Angstrom^-1 and square
    d2 = (scattering_vector_length*1e-10)**2
    # the 0.25 is there because the 's' used by the aproximation is 'd/2'
    d2 = (d2*0.25)

    pos = jnp.array(sample.pos)
    B = jnp.array(sample.B)
    Z = jnp.array(sample.Z)
    HKL = jnp.array(HKL_list)
    atomsf_jax = jnp.array(atomsf)
    Z_cache, Z_map = compute_sf_cache(sample, d2, atomsf_jax)
    jax_compute_pattern_on_batch_jit = jax.jit(jax_compute_pattern_on_batch)

    for i in numpy.arange(0,sample.natoms, batch_size):        
        inc_r, inc_i = jax_compute_pattern_on_batch(pos[i:i+batch_size], Z[i:i+batch_size], B[i:i+batch_size], HKL, Z_cache, Z_map, d2)
        F_r += inc_r
        F_i += inc_i
        print('%3.1f percent done' % (i/sample.natoms*100))
        
    F = numpy.zeros(HKL_list.shape[0],dtype=numpy.complex64)
    F.real = numpy.array(F_r)
    F.imag = numpy.array(F_i)
    pat = Pattern(F, HKL_list)
    pat.ints = numpy.abs(F)**2
    return pat

def compute_sf_cache(sample, d2, atomsf_jax):
    unique_Z = numpy.unique(sample.Z)
    Z_cache = []
    Z_map = numpy.zeros(int(numpy.max(unique_Z)+1), dtype=int)
    print('Computing Z cache size ',len(unique_Z))
    for Z in unique_Z:
        Z_map[Z] = len(Z_cache)
        Z_cache.append(atomsf_jax[Z,0]*jnp.exp(-atomsf_jax[Z,0+4]*d2) + 
        atomsf_jax[Z,1]*jnp.exp(-atomsf_jax[Z,1+4]*d2) + 
        atomsf_jax[Z,2]*jnp.exp(-atomsf_jax[Z,2+4]*d2) + 
        atomsf_jax[Z,3]*jnp.exp(-atomsf_jax[Z,3+4]*d2) +
        atomsf_jax[Z,8])        
    Z_cache = jnp.array(numpy.array(Z_cache))
    Z_map = jnp.array(Z_map)
    return Z_cache, Z_map

scatt_factor_hash = {}
B_factor_hash = {}
def scatt_factor(d, Z, B):
    key = '%d-%f' % (Z,B)
    global scatt_factor_hash
    if key in scatt_factor_hash:
        return scatt_factor_hash[key]

    sf = numpy.zeros_like(d)
    # Convert to Angstrom^-1 and square
    d2 = (d*1e-10)**2
    # the 0.25 is there because the 's' used by the aproxumation is 'd/2'
    for i in range(4):
        sf += atomsf[Z][i]*numpy.exp(-(atomsf[Z][i+4]+B)*d2*0.25)               
    sf += atomsf[Z][8]*numpy.exp(-B*d2*0.25);
    scatt_factor_hash[key] = sf
    return sf;    

def calculate_thomson_correction(det, exp):
    det.thomson_correction = numpy.zeros((det.nz,det.ny,det.nx))
    r0 = 2.81794e-15; # classical electron radius = e^2/(m*c^2)
    if(exp.polarization != 'ignore'):
        raise NotImplementedError
    det.thomson_correction[:] = r0**2

def calculate_pixel_solid_angle(det):
    incidence = numpy.arctan(numpy.sqrt((det.px)**2+(det.py)**2)/det.pz)
    det.solid_angle =  det.pixel_height*det.pixel_width/det.pz**2*numpy.cos(incidence)**3

def calculate_photons_per_pixel(pattern, opts):
    det = opts.detector
    exp = opts.experiment
    if(det.thomson_correction is None):
        calculate_thomson_correction(det,exp)
    if(det.solid_angle is None):
        calculate_pixel_solid_angle(det)
    det.photons_per_pixel = det.thomson_correction*det.solid_angle*pattern.ints.reshape((det.ny,det.nx))*exp.beam_intensity

def calculate_noiseless_detector_output(det, exp):
    h_c = 1.98644521e-25; # Planck's constant * the speed of light
    photon_energy = h_c/exp.wavelength
    electrons_per_photon = photon_energy/det.electron_hole_production_energy
    ADC_constant = 1
    #if det.linear_full_well and det.maximum_value:
    #    ADC_constant = det.maximum_value/det.linear_full_well

    det.noiseless_output = det.photons_per_pixel *  det.quantum_efficiency * electrons_per_photon * ADC_constant
    
def simulate_shot(sample, opts):
    if opts.sf_filename:
        pattern = load_pattern_from_file(opts.detector, opts.sf_filename)
    else:
        HKL_list = get_HKL_list(opts)
        rot = apply_orientation_to_HKL_list(HKL_list, opts)
        pattern = compute_sf(sample, HKL_list, opts)
        pattern.rot = rot
        pattern.HKL_list = HKL_list

    calculate_photons_per_pixel(pattern, opts)
    return pattern

# Convenience functions
def get_molecule_from_opts(opts):
    mol = get_molecule(opts)
    return mol
    
def get_molecule_from_atoms(atomic_numbers = None, atomic_positions = None):
    mol = alloc_mol()
    for j,(x,y,z) in zip(atomic_numbers, atomic_positions):
        add_atom_to_mol(mol, int(j), float(x), float(y), float(z))
    return mol

def get_atoms_from_molecule(mol):
    pos_img = sp_image_alloc(3, mol.natoms, 1)
    array_to_image(mol.pos, pos_img)
    pos = numpy.array(pos_img.image.real, dtype=numpy.float64)
    pos = pos.reshape((mol.natoms,3))
    #print "from molecule: Lz=%e Ly=%e Lx=%e" % (pos[:,1].max()-pos[:,1].min(),pos[:,2].max()-pos[:,2].min())
    anum_img = sp_image_alloc(mol.natoms, 1, 1)
    iarray_to_image(mol.atomic_number, anum_img)
    anum = numpy.array(anum_img.image.real, dtype=numpy.int32)
    anum = anum.reshape(anum.size)
    sp_image_free(pos_img)
    sp_image_free(anum_img)
    return anum, pos

def fetch_pdb(pdb_id):
    url = "http://www.rcsb.org/pdb/files/%s.pdb.gz" % str(pdb_id)
    filename = "./%s.pdb" % str(pdb_id)
    response = urlopen(url)
    compressedFile = BytesIO()
    compressedFile.write(response.read())
    compressedFile.seek(0)
    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
    with open(filename, 'wb') as outfile:
        outfile.write(decompressedFile.read())
    return filename

config = configparser.ConfigParser(defaults={'detector_center_x': '0', 
    'detector_center_y': '0', 'detector_center_z': '0',
    'detector_binning': '1', 'detector_binning_x': '1',
    'detector_binning_y': '1', 'detector_binning_z': '1', 'detector_depth': '0', 
    'detector_pixel_depth': '0', 'use_fft_for_sf': '0', 'use_nfft_for_sf': '0', 
    'b_factor': '0', 'origin_to_com': '0',
    'crystal_size_a': '1', 'crystal_size_b': '1', 'crystal_size_c': '1',
    'crystal_cell_a': '1', 'crystal_cell_b': '1', 'crystal_cell_c': '1', 
    'crystal_cell_alpha': '90', 'crystal_cell_beta': '90', 'crystal_cell_gamma': '90',
    'use_cuda': '1', 'number_of_patterns': 1, 'vectorize': 1, 'wavelength_samples': '5',
    'random_seed': '-1', 'output_noiseless_photons': 1, 'output_count': 1, 
    'output_scattering_factors': '0',
    'output_real_space': '0', 'verbosity_level': '1', 'experiment_polarization': 'ignore',
    'random_orientation': '0', 'phi': '0', 'theta': '0', 'psi': '0','use_jax': '1'})


# det = Detector({'detector_center_x':0, 'detector_center_y':0, 'detector_center_z':0,
#           'detector_distance':1, 'detector_width':2, 'detector_height':2,
#           'detector_depth':0, 'detector_pixel_width':0.001, 'detector_pixel_height':0.001,
#           'detector_pixel_depth':0.0,'detector_quantum_efficiency':1,
#           'detector_electron_hole_production_energy':0 })
# calculate_pixel_solid_angle(det)
# print(det.solid_angle.sum())
# print(4*numpy.pi/6)
# print(det.solid_angle.sum()-4*numpy.pi/6)


# We need to add a dummy section as required by configparser
#with open('/Users/filipe/src/spsim/examples/spsim.conf') as stream:
with open('/Users/filipe/src/spsim/examples/spsim_fast.conf') as stream:
#with open('/Users/filipe/src/spsim/examples/3wun.conf') as stream:
#with open('/Users/filipe/src/spsim/examples/1aon.conf') as stream:
    config.read_string("[top]\n" + stream.read()) 
print(list(config['top'].keys()))


opts = Options(config['top'])
#sample = get_Sample_from_pdb('/Users/filipe/src/spsim/examples/DNA.pdb')
#sample = get_Sample_from_pdb('/Users/filipe/src/spsim/examples/3wun.pdb')
sample = get_Sample_from_pdb('/Users/filipe/src/spsim/examples/1aon.pdb')
#sample = get_Sample_from_pdb('/Users/filipe/src/spsim/examples/carbon.pdb')
calculate_thomson_correction(opts.detector, opts.experiment)
calculate_pixel_solid_angle(opts.detector)
pattern = simulate_shot(sample, opts)
calculate_photons_per_pixel(pattern,opts)
calculate_noiseless_detector_output(opts.detector,opts.experiment);
f = h5py.File('../examples/py_spsim.h5','w')
f['/entry_1/data_1/data'] = opts.detector.photons_per_pixel
f['/entry_1/data_1/F'] = pattern.F
f['/entry_1/data_1/solid_angle'] = opts.detector.solid_angle
f.close()
#plt.imshow(opts.detector.noiseless_output[0]**(1/4))
#plt.colorbar()
#plt.show()
