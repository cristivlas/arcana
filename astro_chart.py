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

def house(rad):
    deg = (math.degrees(rad)-180) % 360
    return deg * 6 / 180 + 1

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
            # center coords        
            xy = [int(i)/2 for i in image.size] 
            
            # ellipse coords
            xy = [i+j for i,j in zip(2*xy, [k*self.radius for k in [-1,-1,1,1]])]

            color = self.args.fg        
            draw.ellipse(xy, outline=color, width=int(1.25 * self.args.antialias_scale))
        return draw

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

        # coords for the center
        x0, y0 = [int(i)/2 for i in image.size]

        # axis
        draw.line((0, y0, size[0], y0), fill=color, width=1 * scale)
        draw.line((x0, 0, x0, size[1]), fill=color, width=1 * scale)

        # houses
        radius = max(self.args.planets_radius, self.args.signs_radius) * scale
        for deg in range(0, 360, 30):
            a = math.radians(deg)
            xy = [0, 0] + [i * radius for i in [math.cos(a), -math.sin(a)]]
            draw.line([i+j for i,j in zip(2*[x0,y0], xy)], fill=color)

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
        xyCenter = 2*[int(i)/2 for i in image.size]
        
        # load font
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
            xy = [i+proj*j for i,j in zip(xyCenter,[-1,-1,1,1])]
            draw.ellipse(xy, outline=self.args.fg)

            a = altAz.az.rad - math.pi / 2

            # plotting coords for center of body 
            xy = [proj * i for i in [math.cos(a), -math.sin(a)]]

            # vector
            #draw.line([i+j for i,j in zip(xyCenter, [0,0] + xy)], fill=self.args.fg, width=1*scale)

            ###
            # draw the celestial body
            rbody = solarSystem[name] * self.radius / 3 + 2 * scale

            # bounding box
            bbox = [int(i+j+k*rbody) for i,j,k in zip(2*xy,xyCenter,[-1,-1,1,1])]

            crop = image.crop(bbox)
            if name == 'sun':
                for rsun in [rbody+15, rbody]:
                    bbox = [int(i+j+k*rsun) for i,j,k in zip(2*xy,xyCenter,[-1,-1,1,1])]
                    draw.ellipse(bbox, outline=self.args.fg, width=scale, fill=self.args.bg)
            else:
                draw.ellipse(bbox, outline=self.args.fg, width=3*scale, fill=self.args.bg)
            crop = Image.blend(image.crop(bbox), crop, 0.45)
            image.paste(crop, bbox)

            # name
            textSize = font.getsize(name)
            draw.text([i+j-k/2 for i,j,k in zip(xy, xyCenter, textSize)], name.capitalize(), font=font, fill=self.args.fg) 

class Signs(ChartElement):
    def __init__(self, args):
        super().__init__(args)
        self.radius = args.signs_radius * args.antialias_scale
        self.icons = Icons(path.join(resource_dir, 'Astro_signs.png'), args.antialias_scale)

    def render(self, image, time, location):
        scale = self.args.antialias_scale
        color = self.args.fg
        background = self.args.bg
        font = ImageFont.truetype(path.join(resource_dir, 'deutschgothic.ttf'), 22 * scale)     

        draw = super().render(image, time, location)
        xyCenter = 2*[int(i)/2 for i in image.size]
        with Stars() as stars:
            for index, name in enumerate(stars.all_signs()):
                name = name.capitalize()
                plot_angles = stars.get_plot_angles(name, obstime=time, location=location)
                if plot_angles[0] < 0:
                    #print ('%s is not visible, altitude=%f rad' % (name, plot_angles[0]))
                    continue
                # project altitude to plane:
                proj = self.radius * math.cos(plot_angles[0])
                
                xy = [i+proj*j for i,j in zip(xyCenter,[-1,-1,1,1])]
                draw.ellipse(xy, outline=color)
                
                a = plot_angles[1]
                print ('%s: house=%d' % (name, house(a)))
                xy = [proj * i for i in [math.cos(a), -math.sin(a)]]

                xySign = [i+j for i,j in zip(xy, xyCenter)]
                
                # vector
                draw.line(xyCenter[:2] + xySign, fill=color, width=int(1*scale))                            
                
                # symbol
                iconMask = self.icons.get(index, self.radius/20)                
                icon = Image.new(iconMask.mode, iconMask.size, background)
                xyIcon = [int(i-j/2) for i,j in zip(xySign, icon.size)]                
                
                draw.ellipse(xyIcon + [i+j for i,j in zip(xyIcon, icon.size)], fill=color, outline=background, width=2*scale)
                image.paste(icon, xyIcon, iconMask)
            
                textSize = font.getsize(name)
                xyText = [i+j for i,j in zip(xySign, [-textSize[0]/2, icon.size[1]/2])]            
                draw.rectangle(xyText + [i+j for i,j in zip(xyText, textSize)], fill=background, outline=background)
                draw.text(xyText, name, font=font, fill=color)
                

             
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
    parser.add_argument('--margin', type=int, default=25)
    parser.add_argument('--planets-radius', default=250, type=int)
    parser.add_argument('--signs-radius', default=350, type=int)

    args = parser.parse_args()

    location = astropy.coordinates.EarthLocation.of_address(args.location)
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
