from PIL import Image

class Icons:
    def __init__(self, filename, scale):        
        self.img = Image.open(filename).convert('RGBA')
        
        # image that holds the icons is assumed to be 4 x 3
        assert self.img.size[0] / self.img.size[1] == 4/3

        self.img = self.img.resize([int(i * scale) for i in [400, 300]], Image.LANCZOS)
        self.scale = scale

    def get(self, index, size, color=None):
        # calculate the row and column where the icon is in the symbols image
        row = int(index / 4)
        col = int(index % 4)
        # ... and crop it out
        cropbox = [coord * 100 * self.scale for coord in [col, row, col+1, row+1]]        
        icon = self.img.crop(cropbox).resize(2 * [int(size * self.scale)])
        if color:
            img = Image.new(icon.mode, icon.size, color)
            icon.paste(img, (0,0), icon)
        return icon
