"""
Plot sky positions of astrological signs -- or rather, plot the homonym constellations. 

Technically the program plots any celestial body for which we can get sky coords.
"""
import argparse
import astropy.coordinates
import astropy.time
import datetime
import json
import math

from dateutil import parser as dtparser
from name_resolver import NameResolver as Stars
from os import path
from PIL import Image, ImageDraw, ImageFont

# constants
antialias_scale = 4


def draw_reference_lines(image, draw, r, fg='black'):
    draw.line((0, r, 2*r, r), fill=fg)
    draw.line((r, 0, r, 2*r), fill=fg)

    # draw houses
    for deg in range(0, 360, 30):
        a = math.radians(deg)
        draw.line([r + i for i in [0, 0, r*math.cos(a), r*math.sin(a)]], fill=fg)


def plot_celestial_bodies(stars, celestialBodies, time, location, radius, fg='black', bg='white'):
    r = radius * antialias_scale

    image = Image.new('RGBA', [2*r] * 2, bg)
    draw = ImageDraw.Draw(image)
    
    draw_reference_lines(image, draw, r)
    draw.ellipse([0,0,2*r,2*r], outline=fg, width=antialias_scale)
    
    # convert datetime to astropy.Time
    time = astropy.time.Time(time)
    print ('Observation time:', time)
    print ('Location:', (location.lat, location.lon))

    for name in celestialBodies:
        name = name.capitalize()
        plot_angles = stars.get_plot_angles(name, obstime=time, location=location)
        if plot_angles[0] < 0:
            print (name, 'is not visible:', plot_angles)
            continue

        # project altitude to plane
        proj = r * math.cos(plot_angles[0])

        draw.ellipse([r - proj, r - proj, r + proj, r + proj], outline=fg)
        a = plot_angles[1]
        x,y = (proj * i for i in [math.cos(a), -math.sin(a)])        
        draw.line([i + r for i in [0,0,x,y]], fill=fg, width=1 * antialias_scale)
    
        text_size = font.getsize(name)
        draw.text([i + r for i in [x-text_size[0]/2, y-text_size[1]/2]], name, font=font, fill=fg)
        
    # anti-alias scale down
    return image.resize(2 * [radius * 2], Image.LANCZOS)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='plot celestial bodies')
    parser.add_argument('input', nargs='*', help='list of celestial bodies')
    parser.add_argument('--radius', default=300, type=int)
    parser.add_argument('--time')
    parser.add_argument('--location', help='location of observer on Earth; example: --location="Seattle,WA"', default='Seattle,WA')
    args = parser.parse_args()

    location = astropy.coordinates.EarthLocation.of_address(args.location)

    if args.time:
        time = dtparser.parse(args.time)
    else:
        time = datetime.datetime.utcnow()

    font = ImageFont.truetype(path.join('res', 'deutschgothic.ttf'), 22 * antialias_scale)     
    input = args.input
    
    with Stars() as stars:
        if not input:
            input = stars.all_signs()
        image = plot_celestial_bodies(stars, input, time, location, args.radius)
        image.show()
