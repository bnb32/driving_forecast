#driving forecast

class Weather:
    pass

def get_coordinates():
    lat=0.0
    lon=0.0
    return lat,lon

def get_base_speed(lat,lon):
    base_speed=30.0
    return base_speed

def get_weather(lat,lon):  
    weather=Weather()
    weather.snow=0
    weather.rain=0
    weather.fog=0
    return weather

def get_speed(base_speed,weather):
    speed=base_speed
    if 0<weather.snow<1: speed*=0.1
    if 0<weather.rain<1: speed*=0.1
    if 0<weather.fog<1: speed*=0.1
    return speed

