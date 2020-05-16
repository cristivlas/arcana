import argparse
import astropy.coordinates
import astropy.time
import datetime
import json
import math

from dateutil import parser as dtparser
from name_resolver import NameResolver as Stars
from os import path
from PIL import Image, ImageColor, ImageDraw, ImageFont
from signs import Icons

resource_dir = 'res'

class Geom:
    @staticmethod
    def translate(xy, coords, type=float):
        assert len(coords) % 2 == 0
        assert len(xy) == 2
        return [type(i + j) for i, j in zip(coords, int(len(coords)/2) * xy)]

    @staticmethod
    def scale(c, coords, type=float):
        return [type(c * i) for i in coords]

    @staticmethod
    def polar_to_cartesian(a, r=1, type=float):
        return Geom.scale(r, [math.cos(a), -math.sin(a)], type)

    @staticmethod
    def polar_to_xy(xy, a, r, type=float):
        return Geom.translate(xy, Geom.polar_to_cartesian(a, r, type), type)

    @staticmethod
    def box(r=1, type=float):
        return Geom.scale(r, [-1,-1,1,1], type)

    # Extend (or contract) a rectangle
    @staticmethod
    def extend(rect, pixelAmount, type=float):
        assert len(rect) == 4
        return [type(i+pixelAmount*j) for i,j in zip(rect, [-1,-1,1,1])]

def house(rad):
    deg = (math.degrees(rad)-180) % 360
    return deg * 6 / 180 + 1

def ascendant(stars, time, location, window=12):
    signs = []
    for s in stars.all_signs():
        a = stars.get_plot_angles(s, time, location)
        if a[0] > 0:
            continue # already risen, skip
        signs.append(s)

    for h in range(0, window*15):
        for s in signs:
            t = astropy.time.Time(time + datetime.timedelta(minutes=h*15))
            #print(t, s)
            a = stars.get_plot_angles(s, t, location)
            if a[0] > 0:              
                return (s,a,t)
    
#
# Helper classes for drawing the chart
#
"""
Base class for Planets chart and Signs chart.
"""
class ChartElement:
    def __init__(self, args):
        self.args = args
        self.radius = None
    
    def render(self, image, time, location):
        draw = ImageDraw.Draw(image)        
        if self.radius:
            xy = Geom.translate(Geom.scale(1/2, image.size), Geom.box(self.radius))
            color = self.args.fg        
            draw.ellipse(xy, outline=color, width=int(1.5 * self.args.antialias_scale))
        return draw

roman = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII', 8:'VIII', 9:'IX', 10:'X', 11:'XI', 12:'XII'}
        
"""
Draw axis and the twelve houses
"""
class ReferenceLines(ChartElement):
    def __init__(self, args):
        super().__init__(args)

    def render(self, image, time, location):
        size = image.size
        color = self.args.fg
        scale = self.args.antialias_scale
        draw = super().render(image, time, location)
        font = ImageFont.truetype(path.join(resource_dir, 'Praetoria D.otf'), 40 * scale)    

        # coords for the center
        x0, y0 = [int(i)/2 for i in image.size]

        # axis
        draw.line((0, y0, size[0], y0), fill=color, width=1 * scale)
        draw.line((x0, 0, x0, size[1]), fill=color, width=1 * scale)

        # houses
        radius = max(self.args.planets_radius, self.args.signs_radius) * scale
        for deg in range(0, 360, 30):
            a = math.radians(deg)
            xy = [0, 0] + Geom.polar_to_cartesian(a, radius)
            draw.line(Geom.translate((x0,y0), xy), fill=color, width = 2 *scale)
        
        radius *= .85
        for i in range(0, 12):
            a = math.radians(195 + 30 * i)
            xy = Geom.translate((x0,y0), Geom.polar_to_cartesian(a, radius))
            text = roman[i+1]
            textSize = font.getsize(text)
            xy = Geom.translate(Geom.scale(-1/2, textSize), xy)
            draw.text(xy, text, color, font)

"""
Utility class
"""
class SolarSystem:
    radii = {    
        'jupiter': 71492000,
        'mars': 3396200,
        'moon': 1738100,
        'mercury': 2439700,
        'neptune': 24764000,
        'pluto': 1188300,
        'saturn': 60268000,
        'sun': 696340000,
        'venus': 6051800,
        'uranus': 25559000
    }
    @staticmethod
    def logscale_radii(font, r):
        solarSystem = {}
        
        sunRadius = SolarSystem.radii['sun']
        plutoRadius = SolarSystem.radii['pluto']
        offs = math.log(plutoRadius)

        # scale relatively to smallest object
        f = 2 * font.getsize('Pluto')[0] / r

        for body in SolarSystem.radii:
            r = SolarSystem.radii[body]
            x = (math.log(r) - offs)/(math.log(sunRadius) - offs) + f
            solarSystem[body] = x

        return solarSystem

    @staticmethod
    def names():
        return list(SolarSystem.radii.keys())
    
"""
Plot planets on the chart
"""
class Planets(ChartElement):
    def __init__(self, args):
        super().__init__(args)
        self.radius = args.planets_radius * args.antialias_scale
    
    def render(self, image, time, location):    
        draw = super().render(image, time, location)
        scale = self.args.antialias_scale
        
        bodies = [(n, astropy.coordinates.get_body(n, time, location)) for n in SolarSystem.names()]
        # sort in reverse distance order
        bodies.sort(key=lambda x: x[1].distance, reverse=True)

        # altitude-azimuth coord frame
        altAzFrame = astropy.coordinates.AltAz(obstime=time, location=location, pressure=0)
        
        # image center coords
        center = Geom.scale(1/2, image.size)
    
        font = ImageFont.truetype(path.join(resource_dir, 'Praetoria D.otf'), 20 * scale)    

        solarSystem = SolarSystem.logscale_radii(font, self.radius)

        for (name, body) in bodies:
            # get altitude and azimuth        
            altAz = (body.transform_to(altAzFrame))
            
            if (altAz.alt.rad < 0):
                #print (name, 'not visible', altAz.alt.deg)
                continue

            # project altitude to plane
            proj = self.radius*math.cos(altAz.alt.rad)

            xy = Geom.translate(center, Geom.box(proj))
            draw.ellipse(xy, outline=self.args.fg)

            a = Stars.azimuth(altAz.az.rad)

            # plotting coords for center of body
            xy = Geom.polar_to_cartesian(a, proj)
            ###
            # draw the celestial body
            rbody = solarSystem[name] * self.radius / 3 + 2 * scale

            # bounding box
            bbox = Geom.translate(center, Geom.translate(xy, Geom.box(rbody)), int)

            crop = image.crop(bbox)
            if name == 'sun':
                for rsun in [rbody+5*scale, rbody]:
                    bbox = Geom.translate(center, Geom.translate(xy, Geom.box(rsun)), int)
                    draw.ellipse(bbox, outline=self.args.fg, width=int(1.5*scale), fill=self.args.bg)
            else:
                draw.ellipse(bbox, outline=self.args.fg, width=3*scale, fill=self.args.bg)
            crop = Image.blend(image.crop(bbox), crop, 0.45)
            image.paste(crop, bbox)

            # draw the name (text)
            textSize = font.getsize(name)
            xyText = Geom.translate(center, Geom.translate(xy, Geom.scale(-1/2, textSize)))
            draw.text(xyText, name, self.args.fg, font)


class Signs(ChartElement):
    def __init__(self, args):
        super().__init__(args)
        self.radius = args.signs_radius * args.antialias_scale
        self.icons = Icons(path.join(resource_dir, 'Astro_signs.png'), args.antialias_scale)

    def _render_asc(self, draw, scale, font, asc):
        if not asc:
            return
        text = 'Ascendant: ' + asc[0] + ' ' + str(asc[2].to_datetime())
        draw.text([10*scale, 70*scale], text, font=font, fill=self.args.fg)
    
    def _render_time_location(self, draw, scale, font, time, location):
        loc = str((location.lon, location.lat))
        draw.text([10*scale, 10*scale], 'Location: ' + loc, font=font, fill=self.args.fg)
        draw.text([10*scale, 40*scale], 'Time: ' + str(time), font=font, fill=self.args.fg)

    def render(self, image, time, location):
        scale = self.args.antialias_scale
        color = self.args.fg
        background = self.args.bg
        font = ImageFont.truetype(path.join(resource_dir, 'deutschgothic.ttf'), 20 * scale)     

        draw = super().render(image, time, location)
        center = Geom.scale(1/2, image.size)

        self._render_time_location(draw, scale, font, time, location)
        
        delayedFuncs = []
        with Stars() as stars:
            self._render_asc(draw, scale, font, ascendant(stars, time, location))                
            for index, name in enumerate(stars.all_signs()):
                name = name.capitalize()
                plot_angles = stars.get_plot_angles(name, obstime=time, location=location)
                if plot_angles[0] < 0:
                    print (name, 'is not visible', plot_angles)
                    continue
                # project altitude to plane:
                proj = self.radius * math.cos(plot_angles[0])
                xy = Geom.translate(center, Geom.box(proj))
                draw.ellipse(xy, outline=color)
                
                a = plot_angles[1] # azimuth
                print ('%s: house=%d' % (name, house(a)))
                xy = Geom.polar_to_cartesian(a, proj)
                xySign = Geom.translate(center, xy)                
                # vector
                draw.line(center + xySign, fill=color, width=int(1*scale))                            

                def draw_symbol(index=index, name=name, xySign=xySign):
                    iconMask = self.icons.get(index, self.radius/20)                
                    icon = Image.new(iconMask.mode, iconMask.size, background)
                    xyIcon = [int(i-j/2) for i,j in zip(xySign, icon.size)]                
    
                    bbox = xyIcon + [i+j for i,j in zip(xyIcon, icon.size)]                
                    bbox = Geom.extend(bbox, -6 * scale) #shrink
                    draw.ellipse(bbox, fill=color, outline=background, width=2*scale)
                    image.paste(icon, xyIcon, iconMask)
                
                    textSize = font.getsize(name)
                    xyText = [i+j for i,j in zip(xySign, [-textSize[0]/2, icon.size[1]/2])]            
                    bbox = xyText + [i+j for i,j in zip(xyText, textSize)]
                    bbox = Geom.extend(bbox, 2 * scale, int)
                    crop = image.crop(bbox)                  
                    draw.rectangle(bbox, fill=background, outline=color, width=2*scale)
                    draw.text(xyText, name, font=font, fill=color)
                    crop = Image.blend(image.crop(bbox), crop, 0.25)
                    image.paste(crop, bbox[:2])

                #delayedFuncs.append(draw_symbol)
                draw_symbol(index, name, xySign)
            
            for f in delayedFuncs:
                f()

             
def astro_chart(args, time, location):
    radius = max(args.planets_radius, args.signs_radius) + args.margin
    # image size:
    size = [2*radius] * 2
    radius *= args.antialias_scale
    # empty image to draw on:
    image = Image.new('RGBA', [2*radius] * 2, args.bg)

    ReferenceLines(args).render(image, time, location)    
    Signs(args).render(image, time, location)
    Planets(args).render(image, time, location)
    return image.resize(size, Image.LANCZOS)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='plot astrological chart')
    parser.add_argument('--time', help='observation time (time of birth)')
    parser.add_argument('--location', help='location of observer on Earth', default='Seattle, WA')
    parser.add_argument('--antialias-scale', type=int, default=4)
    parser.add_argument('--fg', help='foreground color', default='black')
    parser.add_argument('--bg', help='background color', default='wheat')
    parser.add_argument('--margin', type=int, default=50)
    parser.add_argument('--planets-radius', default=250, type=int)
    parser.add_argument('--signs-radius', default=350, type=int)
    parser.add_argument('--out')

    args = parser.parse_args()

    location = astropy.coordinates.EarthLocation.of_address(args.location)
    #astropy.coordinates.solar_system_ephemeris.set('de432s')    
    astropy.coordinates.solar_system_ephemeris.set('jpl')    

    if args.time:
        time = dtparser.parse(args.time)
    else:
        time = datetime.datetime.utcnow()
    time = astropy.time.Time(time)

    print ('Observation time:', time)
    print ('Location:', (location.lat, location.lon))
    print ('Sun in', astropy.coordinates.get_sun(time).get_constellation())

    image = astro_chart(args, time, location)
    image.show()
    if args.out:
        image.save(args.out)

