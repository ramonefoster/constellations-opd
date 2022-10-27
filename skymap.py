import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from skyfield.api import Star, load, wgs84, N, S, W, E
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from datetime import datetime
from pytz import timezone
import matplotlib.transforms as tx

def pol2cart(rho, phi):
    """rho: zenith  |   phi: azimuth"""
    x = rho * np.cos(np.radians(phi+4))
    y = rho * np.sin(np.radians(phi+4))
    # x=f*rho*sin(A−Δ)+x0
    # y=f*rho*cos(A−Δ)+y0
    #f = 3.365 piexels/degree
    # z = zenith estrela
    # A = azimuth estrela
    # Δ = rotacao da camera
    # x0 e y0 = coordenadas do zenith da camera
    x=3.365*x+12
    y=3.365*y
    return(x, y)

def get_az_alt(coordRA, coordDEC, lst, latitude):
    DEG = 180 / np.pi
    RAD = np.pi / 180.0
    H = (lst - coordRA) * 15

    #altitude calc
    sinAltitude = (np.sin(coordDEC * RAD)) * (np.sin(latitude * RAD)) + (np.cos(coordDEC * RAD) * np.cos(latitude * RAD) * np.cos(H * RAD))
    altitude = np.arcsin(sinAltitude) * DEG #altura em graus

    #azimuth calc
    y = -1 * np.sin(H * RAD)
    x = (np.tan(coordDEC * RAD) * np.cos(latitude * RAD)) - (np.cos(H * RAD) * np.sin(latitude * RAD))

    #This AZposCalc is the initial AZ for dome positioning
    AZposCalc = np.arctan2(y, x) * DEG
    #converting neg values to pos
    if (AZposCalc.any() < 0):
        AZposCalc = AZposCalc + 360 

    x, y = pol2cart(90-altitude, 50+AZposCalc)  
    print((90-altitude, AZposCalc)) 

    return(x,y)

# Timezone
AMS = timezone('America/Sao_Paulo')
ts = load.timescale()
time_at = ts.from_datetime(AMS.localize(datetime.now()))

# 180 = South 0 = North
degrees = 50.0

opd_local = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).at(time_at)

#local sidereal time
lst = wgs84.latlon(-22.5344, -45.5825, elevation_m=1800).lst_hours_at(time_at)

position = opd_local.from_altaz(alt_degrees=90, az_degrees=0)

# An ephemeris from the JPL provides Sun and Earth positions.
eph = load('de421.bsp')
sun = eph['sun']
earth = eph['earth']

# Zeniths coordinates
center_star = Star(ra_hours=lst, dec_degrees=-22.5344)

# Hipparcos catalog
with load.open(hipparcos.URL) as f:
    stardata = hipparcos.load_dataframe(f)

# Download a list of constellations from Stellarium
url_constellations = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
                    '/skycultures/western_SnT/constellationship.fab')

with load.open(url_constellations) as f:
    constellations = stellarium.parse_constellations(f)

def generate_constellation_lines(data):
    edges = [edge for name, edges in data for edge in edges]
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]
    xy1 = stardata[['x', 'y']].loc[edges_star1].values
    xy2 = stardata[['x', 'y']].loc[edges_star2].values
    
    # The constellation lines will each begin at the x,y of one star and end
    # at the x,y of another.  We have to "rollaxis" the resulting coordinate
    # array into the shape that matplotlib expects.
    return np.rollaxis(np.array([xy1, xy2]), 1)

# Make a projection and get the the x and y coordinates that each star will have on the plot.
projection = build_stereographic_projection(position)
star_positions = earth.at(time_at).observe(Star.from_dataframe(stardata))
stardata['x'], stardata['y'] = projection(star_positions)

#Adjust fish-eye Allsky distotion
stardata['x'], stardata['y'] = get_az_alt(stardata.ra_hours,stardata.dec_degrees,lst,-22.5344)

# Create a True/False mask marking the stars bright enough to be
# included in our plot.  And go ahead and compute how large their
# markers will be on the plot.
limiting_magnitude = 5.0
bright_stars = (stardata.magnitude <= limiting_magnitude)
# limiting_east = float(lst)-6
# limiting_west = float(lst)+6
# visible_stars = (limiting_east<=(stardata.ra_hours)) & (limiting_west>=(stardata.ra_hours))

magnitude = stardata['magnitude'][bright_stars]
marker_size = (0.7 + limiting_magnitude - magnitude) ** 2.0

# Get allsky from URL (640x480)
img = plt.imread("http://200.131.64.237:8090/img/allsky_picole.png")
fig, ax = plt.subplots(figsize=[640, 480])
tr = tx.Affine2D().rotate_deg(0)
ax.imshow(img, extent=[-320, 320, -240, 240])

# Plot constellation lines
ax.add_collection(LineCollection(generate_constellation_lines(constellations),
                colors='gold', linewidths=1, zorder=1, alpha=0.4))

# Plot Stars only at daytime
ax.scatter(stardata['x'][bright_stars], stardata['y'][bright_stars],
    s=marker_size, color='r', alpha=.2)

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
plt.show()

fig.savefig('OPD_constellations.png', bbox_inches='tight', transparent=True)


