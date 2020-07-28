from math import sqrt

def get_int_points_in_circle(center=(0,0), radius=3):
    center_x = center[0]
    center_y = center[1]
    
    coords=[]

    for i in range(-radius, radius+1):
        for j in range(-radius, radius+1):
            if sqrt((i - center_x)**2 + (j-center_y)**2) <= radius:
                coords.append((i,j))

    return coords

print(len(get_int_points_in_circle()))
