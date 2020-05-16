"""
Get sky coordinates of cellestial body and cache it.

The SIMBAD database, used as a name resolver by astropy,
seems to be missing many constellations. The workaround is to
use the brightest star as a stand-in for the constellation.

"""
import appdirs
import astropy
import json
import math
import pydoc
from os import makedirs, path, unlink

# Custom JSON encoder for Lat / Lon
class CoordsEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, astropy.coordinates.angles.Latitude):
            return { 'value': obj.to_string(), '__class__': 'astropy.coordinates.angles.Latitude' }
        elif isinstance(obj, astropy.coordinates.angles.Longitude):
            return { 'value': obj.to_string(), '__class__': 'astropy.coordinates.angles.Longitude' }
        return super.default(self, obj)

    @staticmethod
    def decode(coords):
        def instance(d):
            c = pydoc.locate(d['__class__'])
            return c(d['value'])
        latLon = [instance(d) for d in coords]
        return astropy.coordinates.SkyCoord(*latLon)

class NameResolver:
    def __init__(self, zodiac_file='zodiac.json'):
        self.data = {}
        appName = 'arcana'
        appAuth = 'umbralist'
        self.path = appdirs.user_cache_dir(appName, appAuth)
        makedirs(self.path, exist_ok=True)
        self.path = path.join(self.path, 'stars.json')
        if path.exists(self.path):
            self._load_data()
        self._load_zodiac(zodiac_file)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._save_data()

    def _load_data(self):
        #print ('reading', self.path)
        with open(self.path) as f:
            try:
                self.data = json.load(f)
            except:                
                unlink(self.path)
                print ('unlinked cache file', self.path)

    def _save_data(self):
        #print ('writing', self.path)
        with open(self.path, 'w') as f:
            json.dump(obj=self.data, fp=f, sort_keys=True, indent=4, cls=CoordsEncoder)

    def _load_zodiac(self, filename):
        with open(filename) as f:
            zodiac = json.load(f)
        self.signs = {}
        for sign in zodiac:
            self.signs[sign['name']] = sign['stars']

    """ if name specifies a constellation, look it up by the brightest star """ 
    def _lookup_zodiac(self, name):
        if name in self.signs:
            star = self.signs[name][0]
            #print ('query-star', star)
            body = astropy.coordinates.SkyCoord.from_name(star)
            constellation = astropy.coordinates.get_constellation(body)
            print (star, 'in', astropy.coordinates.get_constellation(body))
            if body and constellation==name:
                return body

    def lookup(self, name):
        if name in self.data:
            #print ('decoding cached value', name)
            return CoordsEncoder.decode(self.data[name])
        else:
            try:
                #print ('query-name', name)
                coords = astropy.coordinates.SkyCoord.from_name(name)
            except astropy.coordinates.name_resolve.NameResolveError:
                coords = self._lookup_zodiac(name)
            # cache the coords
            self.data[name] = [coords.ra, coords.dec]
            return coords

    @staticmethod
    def azimuth(rad):
        assert rad >= 0 and rad < 2*math.pi
        return (2*math.pi-rad + math.pi/2)        

    def get_plot_angles(self, name, obstime, location):
        skyCoord = self.lookup(name)
        # get altitude and azimuth
        altAz = (skyCoord.transform_to(astropy.coordinates.AltAz(location=location, obstime=obstime)))

        # adjust Azimuth (which is measured from North) to trigonometric zero (which is East)
        return [altAz.alt.rad, self.azimuth(altAz.az.rad)]

    def all_signs(self):
        return list(self.signs.keys())