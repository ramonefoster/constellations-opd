import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from skyfield.api import Star, load, wgs84, N, S, W, E
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from datetime import datetime, timedelta
from pytz import timezone
import matplotlib.transforms as tx
import threading
import os
import time
import math

class SkyMap():    
    def __init__(self):
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

    def pol2cart(self, rho, phi):
        """Convert polar coordinates to cartesian"""
        """rho: zenith  |   phi: azimuth"""
        x = rho * np.cos(np.radians(phi+52.6+82.0))
        y = rho * np.sin(np.radians(phi+52.6+82.0))
        # x=f*rho*cos(A-?)+x0
        # y=f*rho*sin(A-?)+y0
        #f = 3.365 piexels/degree
        # z = zenith estrela
        # A = azimuth estrela
        # ? = rotacao da camera
        # x0 e y0 = coordenadas do zenith da camera
        x=3.365*x+8.5
        y=3.365*y-5
        return(x, y)

   
    def get_az_alt(self, coordRA, coordDEC, lst, latitude):
        """Convert equatorial coordinates to horizontal"""
        DEG = 180 / math.pi
        RAD = math.pi / 180.0        
        H = (lst - coordRA) * 15

        #altitude calc
        sinAltitude = (np.sin(coordDEC * RAD)) * (np.sin(latitude * RAD)) + (np.cos(coordDEC * RAD) * np.cos(latitude * RAD) * np.cos(H * RAD))
        altitude = np.arcsin(sinAltitude) * DEG #altura em graus
        #Zenith
        zenith = 90-altitude

        #azimuth calc
        y = -1 * np.sin(H * RAD)
        x = (np.tan(coordDEC * RAD) * np.cos(latitude * RAD)) - (np.cos(H * RAD) * np.sin(latitude * RAD))

        #Azimuth calc
        azimuth = np.arctan2(y, x) * DEG
        #converting neg values to pos
        if (azimuth.any() < 0):
            azimuth = azimuth + 360
        
        x, y = self.pol2cart(zenith, azimuth)

        return(x,y)
    
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

    def generate(self, allsky_path):
        while True:
            # Timezone
            AMS = timezone('America/Sao_Paulo')
            ts = load.timescale()
            if 6<datetime.now().hour<18:
                time_at = ts.from_datetime(AMS.localize(datetime(2022, 7, 4, 23, 36, 17)))
            else:
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
                self.stardata['x'], self.stardata['y'] = self.get_az_alt(self.stardata.ra_hours,self.stardata.dec_degrees,lst,-22.5344)

                # Create a True/False mask marking the stars bright enough to be
                # included in our plot.  And go ahead and compute how large their
                # markers will be on the plot.
                limiting_magnitude = 4.5
                bright_stars = (self.stardata.magnitude <= limiting_magnitude)
                magnitude = self.stardata['magnitude'][bright_stars]
                marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0

                # Get allsky from URL (640x480)
                img = plt.imread(r"C:\Users\teste\Desktop\coopd\public\img\allsky_picole.jpg")
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

                #plt.show()
                print(f"{datetime.now().hour}:{datetime.now().minute} Creating map {time_at.utc_strftime('%H:%M')}")
                fig.savefig(allsky_path, bbox_inches='tight', transparent=True, dpi=250)
                plt.clf()
                plt.close(fig)
                self.reload_eph()
                time.sleep(30)
            except Exception as e:
                print("ERRO SKYMAP:", str(e))

x = SkyMap()
x.generate(r"C:\Users\teste\Desktop\coopd\public\img\allsky_picole.png")

