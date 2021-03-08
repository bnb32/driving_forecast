#driving forecast

class Weather:
    pass

class DriveRoute:
    def __init__(self,lats,lons,dists,times):
        self.lats = lats
        self.lons = lons
        self.dists = dists
        self.times = times
        self.speeds = [self.dists[i]/self.times[i] for i in range(len(self.times))]

    def get_weather(self,lat,lon,time):  
        weather=Weather()
        weather.light_rain = False
        weather.light_snow = False
        weather.heavy_rain = False
        weather.heavy_snow = False 
        weather.low_visibility = False 
        return weather

    def get_speed(self,weather,base_speed):
        speed = base_speed
        if weather.heavy_snow: speed*=0.25
        if weather.light_snow: speed*=0.08
        if weather.heavy_rain: speed*=0.1
        if weather.light_rain: speed*=0.08
        if weather.low_visibility: speed*=0.11
        return speed

    def dt(self,i,t):
        weather = self.get_weather(self.lats[i],self.lons[i],t)
        speed = self.get_speed(weather,self.speeds[i])
        return self.dists[i]/speed

    def get_drive_time(self):
        time = 0
        for i in range(len(self.dists)):
            time+=self.dt(i,time)
        return time

