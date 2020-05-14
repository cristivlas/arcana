import argparse
import json
import datetime
from dateutil import parser as dtparser

def get_decan(zodiac, time):
    year = datetime.datetime.utcnow().year
    time = datetime.datetime(year=year, month=time.month, day=time.day, hour=time.hour, minute=time.minute)
    for sign in zodiac:
        decans = sign['decans']
        assert len(decans) == 3
        for i, d in enumerate(decans):
            dates = [dtparser.parse(t) for t in d['dates']]
            #print (time, dates)
            if dates[0] <= time and time <= dates[1] + datetime.timedelta(days=1):
                return (sign, i)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', default='zodiac.json')
    parser.add_argument('--dob', required=True)
    args = parser.parse_args()

    timeOfBirth = dtparser.parse(args.dob)
    with open(args.file) as f: 
        zodiacData = json.load(f)
    
    decan = {0 : '1st', 1 : '2nd', 2 : '3rd'}
    astro = get_decan(zodiacData, timeOfBirth)
    if astro:
        print ('%s, %s decan' % (astro[0]['name'], decan[astro[1]]))
    else:
        print ('Unborn:', timeOfBirth)