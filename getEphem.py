import datetime
import ephem
import utils

def get_stars_info(star, sidereal_time):       
    star = ephem.star(star)
    ra = utils.string_to_hours(str(star._ra))
    dec = utils.string_to_degrees(str(star._dec))
    lst = utils.string_to_hours(str(sidereal_time))
    az, alt = utils.get_az_alt(ra, dec, lst, -22.5344)
    return az, alt

def put_planet_name(planet, sidereal_time, allsky_angle):    
    planet.compute()  
    ra = utils.string_to_hours(str(planet.ra))
    dec = utils.string_to_degrees(str(planet.dec))
    lst = utils.string_to_hours(str(sidereal_time))
    az, alt = utils.get_az_alt(ra, dec, lst, -22.5344)
    x, y = utils.pol2cart(90-alt, az, allsky_angle)
    return x, y
    
def put_star_name(star, allsky_angle, sidereal_time):
    az, alt = get_stars_info(star, sidereal_time)
    x, y = utils.pol2cart(90-alt, az, allsky_angle)
    return x, y

def get_observer():
    OPD=ephem.Observer()
    OPD.lat='-22.5344'
    OPD.lon='-45.5825'
    OPD.date = datetime.datetime.utcnow()
    OPD.elevation = 1864 # meters
    OPD.horizon = 0 
    return OPD
    
def astro_coordinates(allsky_angle): 
    stars_list = ['Antares', 'Arcturus', 'Acrux', 'Vega', 'Sirius', 'Betelgeuse', 'Canopus', 'Capella']

    opd = get_observer()
    sidereal_time = opd.sidereal_time()

    #Planets
    mars = ephem.Mars() 
    venus = ephem.Venus()
    jupiter = ephem.Jupiter()
    saturn = ephem.Saturn()

    planets_list = [mars, venus, jupiter, saturn] 

    stars = {}
    planets = {}
    
    for planet in planets_list:
        x, y = put_planet_name(planet, sidereal_time, allsky_angle)
        planets[str(planet.name)] = (x,y)
    
    for star in stars_list:
        x, y = put_star_name(star, allsky_angle, sidereal_time)
        stars[star] = (x,y)

    return planets, stars



