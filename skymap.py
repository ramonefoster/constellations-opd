import os
import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.transforms as tx
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection
from datetime import datetime, timedelta
from pytz import timezone

import utils, getEphem

class SkyMap():    
    def __init__(self):

        self.allsky_angle = 134.6

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

    def generate(self, allsky_img, skymap_img):
        while True:
            # self.reload_eph()
            # Timezone
            AMS = timezone('America/Sao_Paulo')
            ts = load.timescale()
            
            time_at = ts.from_datetime(AMS.localize(datetime.now()))

            opd_local = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).at(time_at)

            #local sidereal time
            lst = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).lst_hours_at(time_at)

            position = opd_local.from_altaz(alt_degrees=90, az_degrees=0)        
            sun = self.eph['sun']
            earth = self.eph['earth']

            try:                
                # Make a projection and get the the x and y coordinates that each star will have on the plot.
                projection = build_stereographic_projection(position)
                star_positions = earth.at(time_at).observe(Star.from_dataframe(self.stardata))
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

                # Plot constellation lines
                ax.add_collection(LineCollection(self.generate_constellation_lines(self.constellations),
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

                # get planets and bright stars names/coordinates
                planets, stars = getEphem.astro_coordinates(self.allsky_angle)  
                for planet in planets:                    
                    x, y = planets[planet]
                    arr_img = plt.imread(f"icons/{planet}.png")
                    if (-280 < x < 280) and (-220 < y < 220):
                        im = OffsetImage(arr_img, zoom=.1)
                        ab = AnnotationBbox(im, (x, y), frameon=False)
                        ax.add_artist(ab)
                        # plt.text(x, y, planet, fontsize=5, color='gold', alpha=0.8)
                
                for star in stars:                    
                    x, y = stars[star]
                    if (-280 < x < 280) and (-220 < y < 220):
                        plt.text(x, y, star, fontsize=5, color='gold', alpha=0.8)

                #plt.show()
                print(f"Creating map {time_at.utc_strftime('%H:%M')}")
                fig.savefig(skymap_img, bbox_inches='tight', transparent=True, dpi=250)
                plt.clf()
                plt.close(fig) 
                             
                self.reload_eph()
                time.sleep(30)
            except Exception as e:
                print("ERRO SKYMAP:", str(e))

#path of original allsky image
allsky_img = "images/allsky.jpg"
#path of destination skymap image
skymap_img = "images/skymap.png"

skyMap = SkyMap()
skyMap.generate(allsky_img, skymap_img)

