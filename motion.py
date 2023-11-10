from tkinter import Tk, Canvas
import math
import numpy as np

tk = Tk()
canvas = Canvas(tk, width=500, height=500)
canvas.pack()

class AngleStorage:
    mousex = 250
    mousey = 250
    
    theta = math.radians(45)
    phi = math.radians(45)
    
    sin_th = math.sin(theta)
    cos_th = math.cos(theta)
    
    sin_ph = math.sin(phi)
    cos_ph = math.cos(phi)
    
    def change_perspective(self, evt):
        self.theta = -evt.y/100
        self.phi = -evt.x/100
        
        self.sin_th = math.sin(self.theta)
        self.cos_th = math.cos(self.theta)
        
        self.sin_ph = math.sin(self.phi)
        self.cos_ph = math.cos(self.phi)

MPHI = 1.618
R = 1
M = MPHI*R
positions = np.array([[-R, -M, -M],
                      [-R, -M, M],
                      [-R, M, -M],
                      [-R, M, M],
                      [R, -M, -M],
                      [R, -M, M],
                      [R, M, -M],
                      [R, M, M]])

positions *= 10

def perspective(matrix):
    sin_ph = angles.sin_ph
    cos_ph = angles.cos_ph
    sin_th = angles.sin_th
    cos_th = angles.cos_th
    
    x = matrix[:, 0]
    y = matrix[:, 1]
    z = matrix[:, 2]
    scale = 5
    
    mx = scale * (x - 10)
    my = scale * (y - 10)
    mz = scale * (10 - z)
    
    rx, ry, rz = sin_ph * my + cos_ph * mx, sin_th * mz + cos_th * (cos_ph * my - sin_ph * mx), cos_th * mz + sin_th * (cos_ph * my - sin_ph * mx)
    
    distance = 100 + 400
    
    px = 300*rx/(rz + distance) + 250
    py = 300*ry/(rz + distance) + 250
    
    return np.array([px, py])

angles = AngleStorage()
canvas.bind_all('<Motion>', angles.change_perspective)

def line(index1, index2):
    canvas.create_line(new[0, index1], new[1, index1], new[0, index2], new[1, index2])

while True:
    try:
        tk.update_idletasks()
        tk.update()
        canvas.delete('all')
        global new
        new = perspective(positions)

        line(0, 1)
        line(0, 2)
        line(0, 4)
        line(1, 3)
        line(1, 5)
        line(2, 3)
        line(2, 6)
        line(3, 7)
        line(4, 5)
        line(4, 6)
        line(5, 7)
        line(6, 7)
    except:
        break


