import os
import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.transforms as tx
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos, stellarium
from skyfield.api import load_constellation_names
from skyfield.projections import build_stereographic_projection
from skyfield.framelib import ecliptic_frame

from datetime import datetime, timedelta
from pytz import timezone

import utils, getEphem

class SkyMap():    
    def __init__(self):
        self.allsky_angle = 0

        self.stardata = None
        # Download a list of constellations from Stellarium
        self.url_constellations = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
                                '/skycultures/western_SnT/constellationship.fab')        
        # Hipparcos catalog
        with load.open('hip_main.dat') as f:
            self.stardata = hipparcos.load_dataframe(f)  

        
        with load.open(self.url_constellations) as f:
            self.constellations = stellarium.parse_constellations(f)
            
        # An ephemeris from the JPL provides Sun and Earth positions.
        self.eph = load('de421.bsp')
    
    def constelattions_name_coord(self, lst):
        """
        Get the first star coordinates of the constellation, 
        converts to cartesian and adds to a dic with its name
        """
        constellation_coord = {}
        d = dict(load_constellation_names())
        for stars in self.constellations:
            star_hip = [coord[0] for coord in stars[1]]                
            try:
                ra_total = 0
                dec_total = 0
                for star in star_hip:
                    ra_total += self.stardata.loc[star].ra_hours
                    dec_total += self.stardata.loc[star].dec_degrees
                    self.stardata.loc[star].dec_degrees

                # Calc mean Right Ascension and mean declination    
                ra_mean = ra_total/len(star_hip)
                dec_mean = dec_total/len(star_hip)
                # convert to cartesian to insert on Image
                az, elev = utils.get_az_alt(ra_mean, dec_mean, lst, -22.5344)
                x, y = utils.pol2cart(90-elev, az, self.allsky_angle)
                constellation_coord[d[stars[0]]] = (x, y)    
            except:
                pass
        return constellation_coord
            
    def reload_eph(self):
        t = os.path.getmtime('de421.bsp')
        if (datetime.now() - datetime.fromtimestamp(t)) > timedelta(minutes=30):
            # An ephemeris from the JPL provides Sun and Earth positions.
            self.eph = load('de421.bsp')    
    
    def generate_constellation_lines(self, data):
        edges = [edge for name, edges in data for edge in edges]
        edges_star1 = [star1 for star1, star2 in edges]
        edges_star2 = [star2 for star1, star2 in edges]
        xy1 = self.stardata[['x', 'y']].loc[edges_star1].values
        xy2 = self.stardata[['x', 'y']].loc[edges_star2].values
        
        # The constellation lines will each begin at the x,y of one star and end
        # at the x,y of another.  We have to "rollaxis" the resulting coordinate
        # array into the shape that matplotlib expects.
        return np.rollaxis(np.array([xy1, xy2]), 1)

    def generate(self, original_img, destination_img, angles):
        while True:
            for allsky_img, skymap_img, allsky_angle in zip(original_img, destination_img, angles):
                # Timezone
                self.allsky_angle = allsky_angle
                AMS = timezone('America/Sao_Paulo')
                ts = load.timescale()

                if 6<datetime.now().hour<18:
                    #AllSky at OPD is offline
                    is_online = False
                else:
                    is_online = True

                time_at = ts.from_datetime(AMS.localize(datetime.now()))
                
                opd_local = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).at(time_at)

                #local sidereal time
                lst = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).lst_hours_at(time_at)

                position = opd_local.from_altaz(alt_degrees=90, az_degrees=0)        
                sun = self.eph['sun']
                earth = self.eph['earth']
                moon = self.eph['moon']

                try:                
                    # Make a projection and get the the x and y coordinates that each star will have on the plot.
                    projection = build_stereographic_projection(position)
                    e = earth.at(time_at)
                    star_positions = e.observe(Star.from_dataframe(self.stardata))
                    self.stardata['x'], self.stardata['y'] = projection(star_positions)

                    #Adjust fish-eye Allsky distotion                
                    azimuth, elevation = utils.get_az_alt(self.stardata.ra_hours,self.stardata.dec_degrees,lst,-22.5344)
                    self.stardata['x'], self.stardata['y'] = utils.pol2cart(90-elevation, azimuth, self.allsky_angle)
                    
                    # Create a True/False mask marking the stars bright enough to be
                    # included in our plot.  And go ahead and compute how large their
                    # markers will be on the plot.
                    limiting_magnitude = 4.5
                    bright_stars = (self.stardata.magnitude <= limiting_magnitude)
                    magnitude = self.stardata['magnitude'][bright_stars]
                    marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0
                    
                    # Get allsky from URL (640x480)
                    img = plt.imread(allsky_img)
                    fig, ax = plt.subplots(figsize=[6, 6])
                    tr = tx.Affine2D().rotate_deg(0)
                    ax.imshow(img, extent=[-320, 320, -240, 240])

                    center = (0.5, 0.5)  # Center coordinates (x, y)
                    radius = 280.4  # Radius of the circular region

                    # Filter the line segments based on their endpoints' distance from the center
                    filtered_segments = []
                    for segment in self.generate_constellation_lines(self.constellations):
                        dist1 = np.linalg.norm(segment[0] - center)
                        dist2 = np.linalg.norm(segment[1] - center)
                        if dist1 <= radius and dist2 <= radius:
                            filtered_segments.append(segment)

                    # Plot constellation lines
                    ax.add_collection(LineCollection(filtered_segments,
                                        colors='gold', linewidths=.5, zorder=1, alpha=0.4))               

                    # Plot Stars only at daytime
                    ax.scatter(self.stardata['x'][bright_stars], self.stardata['y'][bright_stars],
                        s=marker_size, color='Brown', alpha=.4)            

                    # Adjusts
                    field_of_view_degrees = 180.0
                    angle = np.pi - field_of_view_degrees / 360.0 * np.pi
                    limit = np.sin(angle) / (1.0 - np.cos(angle))

                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                    ax.patch.set_alpha(0.3)
                    ax.set_title(f"Constelações no céu do OPD UTC: {time_at.utc_strftime('%Y %B %d %H:%M')}", color='gold')

                    # Save
                    plt.xlim(-320, 320)
                    plt.ylim(-240, 240)

                    #moon
                    s = e.observe(sun).apparent()
                    m = e.observe(moon).apparent()
                    _, slon, _ = s.frame_latlon(ecliptic_frame)
                    _, mlon, _ = m.frame_latlon(ecliptic_frame)
                    phase = (mlon.degrees - slon.degrees) % 360.0

                    if 0 <= phase < 90:
                        moon_phase = 'new_moon'
                    if 180 <= phase < 270:
                        moon_phase = 'full_moon'
                    else:
                        moon_phase = 'quarter_moon'

                    #insert coordinates to insert it on image
                    constellations_coord = self.constelattions_name_coord(lst)
                    for constellation in constellations_coord:                    
                        x, y = constellations_coord[constellation]                    
                        if (-280 < x < 280) and (-220 < y < 220) and is_online:                        
                            plt.text(x, y, constellation, fontsize=4, color='gold', alpha=0.5)

                    # get planets and bright stars names/coordinates
                    planets, stars, cardinals = getEphem.astro_coordinates(self.allsky_angle)  
                    valid_planets = [(planet, planets[planet]) for planet in planets if (-280 < planets[planet][0] < 280) and (-225 < planets[planet][1] < 230)]
                    valid_stars = [(star, stars[star]) for star in stars if (-280 < stars[star][0] < 280) and (-220 < stars[star][1] < 220)]

                    for planet, (x, y) in valid_planets:
                        if 'Moon' in planet:
                            arr_img = plt.imread(f"icons/{moon_phase}.png")
                        else:
                            arr_img = plt.imread(f"icons/{planet}.png")
                        im = OffsetImage(arr_img, zoom=.1)
                        ab = AnnotationBbox(im, (x, y), frameon=False)
                        ax.add_artist(ab)

                    for star, (x, y) in valid_stars:
                        plt.text(x, y, star, fontsize=5, color='gold', alpha=0.8)
                    
                    for card in cardinals:
                        x, y = cardinals[card]
                        x = max(-310, min(310, x))
                        y = max(-225, min(225, y))
                        plt.text(x, y, card, fontsize=8, color='black', alpha=1, ha='center', va='center', fontweight='bold')
                        plt.text(x, y, card, fontsize=8, color='white', alpha=1, ha='center', va='center')


                    print(f"Creating map {time_at.utc_strftime('%H:%M')}")
                    fig.savefig(skymap_img, bbox_inches='tight', transparent=True, dpi=250)
                    plt.clf()
                    plt.close(fig) 
                                
                    self.reload_eph()
                    time.sleep(30)
                except Exception as e:
                    print("ERRO SKYMAP:", str(e))

#path of original allsky image
allsky_img = [r"C:\Users\servidorOPD\Documents\Programas\SistemaCOOPD\public\img\allsky340c.jpg", r"C:\Users\servidorOPD\Documents\Programas\SistemaCOOPD\public\img\allsky_picole.jpg"]
#path of destination skymap image
skymap_img = [r"C:\Users\servidorOPD\Documents\Programas\SistemaCOOPD\public\img\allsky340c.png", r"C:\Users\servidorOPD\Documents\Programas\SistemaCOOPD\public\img\allsky_picole.png"]
angles = [-220.5, 91.5]

skyMap = SkyMap()
skyMap.generate(allsky_img, skymap_img, angles)

