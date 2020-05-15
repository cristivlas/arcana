
from PIL import Image, ImageDraw, ImageFont
from os import path
from signs import Icons
import json
import math
import argparse

class Zodiac:
    # anti-aliasing trick: draw everything as a larger image
    # then scale it back down. Try setting scale = 1 to see the difference
    scale = 4
 
    # angular step for drawing each astrological sign
    step = math.radians(30)

    """ load fonts and the symbol icons """
    def __init__(self, args, data):
        self.args = args
        self.data = data
        self.font1 = ImageFont.truetype(path.join(args.res, args.main_font), int(Zodiac.scale * args.radius / 15))
        self.font2 = ImageFont.truetype(path.join(args.res, args.secondary_font), int(Zodiac.scale * args.radius / 25))
        self.icons = Icons(path.join(args.res, 'Astro_signs.png'), Zodiac.scale)
        
    def draw(self):
        # size of final image [width, height]
        size = [2 * args.radius] * 2

        # radius for drawing at blown-up scale
        self.radius = Zodiac.scale * args.radius

        image = Image.new('RGBA', [2*self.radius] * 2, args.bg)
        canvas = ImageDraw.Draw(image)

        # draw each sign with its decans
        angle = -math.pi
        for i, sign in enumerate(self.data):
            angle -= Zodiac.step
            self._draw_sign(image, canvas, angle, sign, i)

        # draw some circles...
        for offs in [0, .2 * self.radius, .21 * self.radius]:
            boundingBox = 2 * [offs] + 2 * [2 * self.radius - offs]
            canvas.ellipse(boundingBox, fill=None, outline=self.args.fg, width=3 * Zodiac.scale)
    
        return image.resize(size, Image.LANCZOS)
    
    # draw one sector of the circle that corresponds to an astrological sign
    def _draw_sign(self, image, canvas, angle, sign, index):
        name = sign['name']

        # use the font to compute the [width, height] of the text to draw
        textSize = self.font1.getsize(name)

        # draw the name of the sign
        a = angle + Zodiac.step / 2
        cossin = [math.cos(a), math.sin(a)]

        xy = [int(0.7 * self.radius * i + self.radius - j/2) for i,j in zip(cossin, textSize)]
        canvas.text(xy, name, font=self.font1, fill=self.args.fg)

        # draw divider lines
        xy = [0, 0] + [self.radius * i for i in [math.cos(angle), math.sin(angle)]]
        canvas.line([i + self.radius for i in xy], fill=self.args.fg, width=3 * Zodiac.scale)

        # draw the icon (symbol) for the astrological sign
        iconSize = self.radius / 30 # some arbitrary scaling that looks good        
        icon = self.icons.get(index, iconSize)
        xy = [int(self.radius - iconSize/2 * Zodiac.scale + 0.575 * self.radius * i) for i in cossin]
        image.paste(icon, xy, icon)

        self._draw_decans(image, canvas, angle, sign)

    def _draw_decans(self, image, canvas, angle, sign):
        for j in range(0, 3):
            a = angle + j * Zodiac.step/3
            cossin = [math.cos(a), math.sin(a)]
            xy = [0, 0] + [0.8 * self.radius * i for i in cossin]
            canvas.line([i + self.radius for i in xy], fill=args.fg)

            planet = sign['decans'][j]['subruler']
            textSize = self.font2.getsize(planet)

            a = angle + 5 * Zodiac.step / 6 - j * Zodiac.step / 3
            xy = [.9 * self.radius * i for i in [math.cos(a), math.sin(a)]]
            canvas.text([self.radius + i - j / 2 for i,j in zip(xy, textSize)], planet, font=self.font2, fill=self.args.fg)

            # compute radius for circle around planet name
            # yes, I know the Sun is not a planet.
            r = Zodiac.step * self.radius *.45 if planet=='Sun' else max(textSize)/2
            canvas.ellipse([i + self.radius for i in [xy[0]-r, xy[1]-r, xy[0]+r, xy[1]+r]], fill=None, outline=self.args.fg, width=Zodiac.scale)

            # just a cool effect:
            canvas.ellipse([i + self.radius for i in [xy[0]-r, xy[0]-r, xy[1]+r, xy[1]+r]], fill=None, outline=self.args.fg)

            # OK. Now let's draw the names of the constellations associated with each of the decans
            constellation = sign['decans'][j]['constellation']
            name = constellation['name']

            a -= Zodiac.step / 24
            # Image rotation function want their arguments in degrees
            degrees = -math.degrees(a) % 360
            textSize = self.font2.getsize(name)

            # create a temporary image on which to draw the text
            # (we will then rotate that image and superimpose it on the main image)
            txt = Image.new(image.mode, image.size)
            xy = (int(1.2 * self.radius), self.radius - textSize[1]/2)

            # flip the text for quadrants 2 and 3 so we don't twist 
            # our necks trying to read the text -- comment the code 
            # out to see the difference
            if degrees > 90 and degrees < 270:
                xy = (txt.size[0]-xy[0]-textSize[0], xy[1])
                degrees -= 180 + 2.5
    
            ImageDraw.Draw(txt).text(xy, name, font=self.font2, fill=self.args.fg)
            txt = txt.rotate(degrees)
            image.paste(txt, (0,0), txt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw the Zodiac Wheel')
    parser.add_argument('--bg', help='background color', default='wheat')
    parser.add_argument('--res', help='specify the resources folder', default='res')
    parser.add_argument('--file', help='file containing the json data that describes the Zodiac', default='zodiac.json')
    parser.add_argument('--fg', help='foreground color', default='black')
    parser.add_argument('--main-font', help='font for astrological symbol names', default='deutschgothic.ttf')
    parser.add_argument('--out', help='output file; if none is specified then the image will not be saved')    
    parser.add_argument('--radius', help='radius for the drawing', default=400, type=int)
    parser.add_argument('--secondary-font', help='font for celestial body names', default='Praetoria D.otf')

    args = parser.parse_args()
    with open(args.file) as f:
        zodiac = Zodiac(args, json.load(f))
        image = zodiac.draw()        
        image.show()
        if args.out:
            image.save(args.out)

