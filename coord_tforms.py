import numpy as np

def lat_lon_to_rec(lat, lon):
    return np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)

def rec_to_lat_lon(x, y, z):
    lat = np.arcsin(z)
    lon = np.arctan2(y, x)
    return lat, lon

def rotx(angle):
    rot_matrix = np.array(
        [
            [1, 0, 0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle), np.cos(angle)],
        ]
    )
    return rot_matrix

def roty(angle):
    rot_matrix = np.array(
        [
            [np.cos(angle), 0, np.sin(angle)],
            [0, 1, 0],
            [-np.sin(angle), 0, np.cos(angle)],
        ]
    )
    return rot_matrix

def rotz(angle):
    rot_matrix = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    return rot_matrix

