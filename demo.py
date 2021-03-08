from pyroute import DriveRoute

times = [0.2]*20
lats = [1]*21
lons = [1]*21
dists = [1]*20

route = DriveRoute(lats,lons,dists,times)

print(route.get_drive_time())
