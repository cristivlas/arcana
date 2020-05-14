import argparse
import astropy.coordinates as cs
import astropy.time
import datetime
import math
import os.path

from astropy.constants import astropyconst40 as const
from astropy import units as u
from dateutil import parser as dtparser
from PIL import Image, ImageDraw, ImageFont

antialias_scale = 4

# radii of Solar System main bodies
solar_system_radii = {    
    #'earth': const.R_earth.value,
    'jupiter': const.R_jup.value,
    'mars': 3396200,
    'moon': 1738100,
    'mercury': 2439700,
    'neptune': 24764000,
    'pluto': 1188300,
    'saturn': 60268000,
    'sun': const.R_sun.value,
    'venus': 6051800,
    'uranus': 25559000
}

def logscale_radii(font, r):
    solar_system = {}
    f = 2 * font.getsize('pluto')[0] / r
    sun_radius = solar_system_radii['sun']
    pluto_radius = solar_system_radii['pluto']
    offs = math.log(pluto_radius)
    for body in solar_system_radii:
        r = solar_system_radii[body]
        x = (math.log(r) - offs)/(math.log(sun_radius) - offs) + f
        solar_system[body] = x
    return solar_system

def draw_reference_lines(image, draw, r, color):
    draw.line((0, r, 2*r, r), fill=color)
    draw.line((r, 0, r, 2*r), fill=color)

def plot_solar_system_bodies(names, time, location, radius, fg='black', bg='white'):
    celestialBodies = [(n.lower(), cs.get_body(n.lower(), time, location)) for n in names]
    # sort in reverse distance order
    celestialBodies.sort(key=lambda x: x[1].distance, reverse=True)

    r = radius * antialias_scale
    solar_system = logscale_radii(font, r)

    image = Image.new('RGBA', [2*r] * 2, bg)
    draw = ImageDraw.Draw(image)
    
    draw_reference_lines(image, draw, r, fg)
    draw.ellipse([0,0,2*r,2*r], outline=fg, width=antialias_scale)

    altAzFrame = cs.AltAz(obstime=time, location=location, pressure=0)

    for (name, body) in celestialBodies:
        # get altitude and azimuth        
        altAz = (body.transform_to(altAzFrame))

        if (altAz.alt.rad < 0):
            print (name.capitalize(), 'not visible', altAz.alt.deg)
            continue
        
        # project altitude to plane
        proj = r * math.cos(altAz.alt.rad)
        draw.ellipse([r - proj, r - proj, r + proj, r + proj], outline=fg)

        a = altAz.az.rad - math.pi / 2
        x,y = (proj * i for i in [math.cos(a), math.sin(a)])
        draw.line([i + r for i in [0,0,x,y]], fill=fg)

        # draw the cellestial body
        rbody = (solar_system[name]) * r / 4;

        # 2 extra pixel margin
        rbody += 2 * antialias_scale

        bbox = [int(i + r) for i in [x-rbody, y-rbody, x+rbody, y+rbody]]
        crop = image.crop(bbox)
        if name == 'sun':
            for rsun in [rbody+15, rbody]:
                bbox = [int(i + r) for i in [x-rsun, y-rsun, x+rsun, y+rsun]]
                draw.ellipse(bbox, outline=fg, width=antialias_scale, fill=bg)
        else:
            draw.ellipse(bbox, outline=fg, width=2 * antialias_scale, fill=bg)
        crop = Image.blend(image.crop(bbox), crop, 0.45)
        image.paste(crop, bbox)
        text_size = font.getsize(name)
        draw.text([i + r for i in [x-text_size[0]/2, y-text_size[1]/2]], name.capitalize(), font=font, fill=fg)        
    # anti-alias scale down
    return image.resize(2 * [radius * 2], Image.LANCZOS)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the position of planets (and Sun and Moon)')
    parser.add_argument('input', nargs='+', help='list of planets (may include Sun and Moon)')
    parser.add_argument('--fg', help='foreground color', default='white')
    parser.add_argument('--bg', help='background color', default='navy')
    parser.add_argument('--radius', default=250, type=int)
    parser.add_argument('--time')
    parser.add_argument('--location', help='location of observer on Earth; example: --location="Seattle,WA"', default='Seattle,WA')
    args = parser.parse_args()

    # https://mistifonts.com/self-deception/
    font = ImageFont.truetype(os.path.join('res', 'SelfDeceptionRegular-ALLWA.ttf'), 30 * antialias_scale)    
    #font = ImageFont.truetype(os.path.join('res', 'Praetoria D.otf'), 20 * antialias_scale)
    
    #cs.solar_system_ephemeris.set('de432s')
    cs.solar_system_ephemeris.set('jpl')

    location = cs.EarthLocation.of_address(args.location)
    #print (location)

    if 'all' in args.input:
        args.input = list(solar_system_radii.keys())

    if args.time:
        time = dtparser.parse(args.time)
    else:
        time = datetime.datetime.utcnow()
    time = astropy.time.Time(time)
    print (time)

    image = plot_solar_system_bodies(args.input, time, location, args.radius, bg=args.bg, fg=args.fg)
    image.show()

