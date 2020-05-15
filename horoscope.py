"""
Astrological info based on a Date of Birth
"""
import argparse
import json
import datetime
from dateutil import parser as dtparser

def sign_info(sign, decan):
    d = {0 : '1st', 1 : '2nd', 2 : '3rd'}
    mf = {'M':'Masculine', 'F':'Feminine'}
    info = {
        'sign': sign['name'],
        'decan': d[decan],
        'ruler': sign['decans'][decan]['subruler'],
        'constellation': sign['decans'][decan]['constellation']['name'],
        'duality': mf[sign['duality']],
        'element': sign['element'],
        'personality': sign['trait'],
        'keyword':sign['decans'][decan]['keyword']
    }
    return info

def get_astro(zodiac, time):
    year = datetime.datetime.utcnow().year
    time = datetime.datetime(year=year, month=time.month, day=time.day, hour=time.hour, minute=time.minute)
    for sign in zodiac:
        decans = sign['decans']
        assert len(decans) == 3
        for i, d in enumerate(decans):
            dates = [dtparser.parse(t) for t in d['dates']]
            if dates[0] <= time and time <= dates[1] + datetime.timedelta(days=1):
                return (sign, i)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dob')
    parser.add_argument('--file', default='zodiac.json')
    args = parser.parse_args()

    timeOfBirth = dtparser.parse(args.dob)
    with open(args.file) as f: 
        zodiacData = json.load(f)
    
    astro = get_astro(zodiacData, timeOfBirth)
    if astro:
        print (sign_info(*astro))
    else:
        print ('Unborn:', timeOfBirth)