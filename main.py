# File: waterproj/main.py

import contextlib
import numpy as np
import gc
import math
import copy
import time
import tkinter

tk = tkinter.Tk()
canvas = tkinter.Canvas(tk, width=900, height=900, highlightthickness=0, bg="#FFFFFF")
canvas.pack()


# ==================================================================================================== # 
WIDTH = 20                                                                                             # > These are the main variables for the
LENGTH = 20                                                                                            # > system points and values and handy
DEPTH = 20                                                                                             # > items.
                                                                                                       # 
POSITIONS = np.zeros((WIDTH, LENGTH, DEPTH, 4), dtype=np.int8)                                         # > Note: indices are x, y, z and using
POSITIONS[:10, :10, 15:, 0] = 1                                                                        # > bottom-up layering.
                                                                                                       # 
SQRT2 = 1.4142135623730951                                                                             # 
SUPPRESS = contextlib.suppress(IndexError)                                                             # 
                                                                                                       # 
RIVERBED = np.array([[14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], # 
                     [14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0],          # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 14, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0],          # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],            # 
                     [19, 19, 19, 19, 19, 19, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],])          # 
# ==================================================================================================== # 



def determine_choices(x, y, z):
    # Some variables
    point = POSITIONS[x, y, z]
    directions = np.array([[[5,3,1],   # 
                            [7,3,1],   # > This is the main variable that will get
                            [5,3,1]],  # > edited by this function and returned to
                                       # > the refresh() function
                           [[7,3,1],   # 
                            [9,2,1],   # 
                            [7,3,1]],  # 
                                       # 
                           [[5,3,1],   # 
                            [7,3,1],   # 
                            [5,3,1]]]) # 
    
    # This is some initial editing that you
    # can do without having to iterate
    
    # Momentum
    directions[point[1], point[2], point[3]] += 6
    
    # Here comes the actual iteration
    
    i = [-2, -1, 0, 1, 2] # 
    j = [-2, -1, 0, 1, 2] # > These are the items to iterate over
    k = [-2, -1, 0, 1, 2] # 
    for i_x in i:
        for i_y in j:
            for i_z in k:
                
                # Skip the center case
                if i_x == i_y == i_z == 0:
                    continue
                
                # Fetch the point
                try:
                    p = POSITIONS[i_x+x, i_y+y, i_z+z]
                except IndexError:
                    continue
                
                # Skip if the point is nonexistent
                if p[0] == 0:
                    continue
                # Different versions of the variables (saves computation time)
                
                # Projected version
                t_x = i_x//max(abs(i_x), 1) # > This version has the coordinates for
                t_y = i_y//max(abs(i_y), 1) # > the point (roughly) projected inwards
                t_z = i_z//max(abs(i_z), 1) # 
                # "directions" version
                d_x = t_x + 1 # > These are the variables used to refer
                d_y = t_y + 1 # > to a location in "directions"
                d_z = t_z + 1 # 
                _d_x = 1 - t_x # > These are the inverse positions
                _d_y = 1 - t_y # > of the variables above
                _d_z = 1 - t_z # 
                
                if directions[_d_x, _d_y, _d_z]:
                    directions[_d_x, _d_y, _d_z] += directions[d_x, d_y, d_z]
                
                directions[d_x, d_y, d_z] = 0
                
                # This section uses the momentum of the particle
                
                if t_x == p[1] or t_y == p[2] or t_z == p[3]:
                    continue
                
                if (i_x == 2 and p[1] > -1) or (i_x == -2 and p[1] < 1) or \
                    (i_y == 2 and p[2] > -1) or (i_y == -2 and p[2] < 1) or \
                    (i_z == 2 and p[3] > -1) or (i_z == -2 and p[3] < 1):
                   continue
                
                # Pseudo ending position
                p_x = d_x + p[1] # > These are the coordinates for the
                p_y = d_y + p[2] # > projected square point + momentum
                p_z = d_z + p[3] # 
                _p_x = _d_x - p[1] # 
                _p_y = _d_y - p[2] # > These are the inverses of the above
                _p_z = _d_z - p[3] # 
                
                if (p_x == _p_x) and (p_y == _p_y) and (p_z == _p_z):
                    if directions[_d_x, _d_y, _d_z]:
                        directions[_d_x, _d_y, _d_z] += directions[d_x, d_y, d_z]
                    directions[d_x, d_y, d_z] = 0
                else:
                    if directions[_p_x, _p_y, _p_z]:
                        directions[_p_x, _p_y, _p_z] += directions[p_x, p_y, p_z]
                    directions[p_x, p_y, p_z] = 0
                
                # Real ending position
                r_x = i_x + p[1] # > These are the coordinates for the
                r_y = i_y + p[2] # > projected square point + momentum
                r_z = i_z + p[3] # 
                _r_x = _d_x - p[1] # 
                _r_y = _d_y - p[2] # > These are the inverses of the above
                _r_z = _d_z - p[3] # 
                
                if (r_x == _r_x) and (r_y == _r_y) and (r_z == _r_z):
                    if directions[_d_x, _d_y, _d_z]:
                        directions[_d_x, _d_y, _d_z] += directions[d_x, d_y, d_z]
                    directions[d_x, d_y, d_z] = 0
                else:
                    if directions[_r_x, _r_y, _r_z]:
                        directions[_r_x, _r_y, _r_z] += directions[r_x, r_y, r_z]
                    directions[r_x, r_y, r_z] = 0
                
    left = False
    with SUPPRESS:
        if z <= RIVERBED[x-1, y]:
            left = True
    
    right = False
    with SUPPRESS:
        if z <= RIVERBED[x+1, y]:
            right = True
    
    front = False
    with SUPPRESS:
        if z <= RIVERBED[x, y+1]:
            front = True
    
    back = False
    with SUPPRESS:
        if z <= RIVERBED[x, y-1]:
            back = True
    
    down = False
    with SUPPRESS:
        if z <= RIVERBED[x, y] or z == RIVERBED[x, y]+1:
            down = True
    
    # Border cases
    if x == 0 or left:                                   # 
        tmp = (directions[2, :, :] != 0)                 # > This section eliminates the cases where
        directions[2, :, :] += tmp*directions[0, :, :]*3 # > the point is on the border. This has to
        directions[0, :, :] = 0                          # > be done at the end so that it is not
    elif x == WIDTH - 1 or right:                        # > overwritten.
        tmp = (directions[0, :, :] != 0)                 # 
        directions[0, :, :] += tmp*directions[2, :, :]*3 # 
        directions[2, :, :] = 0                          # 
    if y == 0 or back:                                   # 
        tmp = (directions[:, 2, :] != 0)                 # 
        directions[:, 2, :] += tmp*directions[:, 0, :]*3 # 
        directions[:, 0, :] = 0                          # 
    elif y == LENGTH - 1 or front:                       # 
        tmp = (directions[:, 0, :] != 0)                 # 
        directions[:, 0, :] += tmp*directions[:, 2, :]*3 # 
        directions[:, 2, :] = 0                          # 
    if z == 0 or down:                                   # 
        tmp = (directions[:, :, 2] != 0)                 # 
        directions[:, :, 2] += tmp*directions[:, :, 0]*3 # 
        directions[:, :, 0] = 0                          # 
    if z == DEPTH - 1:                                   # 
        tmp = (directions[:, :, 0] != 0)                 # 
        directions[:, :, 0] += tmp*directions[:, :, 2]*3 # 
        directions[:, :, 2] = 0                          # 
    
    return directions

def refresh():
    for x, layer in enumerate(POSITIONS):
        for y, row in enumerate(layer):
            for z, point in enumerate(row):
                
                # Skip if the point is nonexistent
                if point[0] == 0:
                    continue
                
                p1, p2, p3 = choice_analysis(determine_choices(x, y, z))
                p1 = int(p1[0]) - 1
                p2 = int(p2[0]) - 1
                p3 = int(p3[0]) - 1
                POSITIONS[x+p1, y+p2, z+p3] = np.array([1, p1, p2, p3])
                if not (p1 == p2 == p3 == 0):
                    POSITIONS[x, y, z] = 0

def choice_analysis(paths):
    
    # Starting variables
    paths = copy.deepcopy(paths)
    tmpPaths = copy.deepcopy(paths)
    
    # Sets default final direction
    final_direction = (np.array([1]), np.array([1]), np.array([1]))
    
    # 4 iterations (at most) of updating values
    for i in range(4):
        
        # If not, update some values based on opposite coordinates
        i = [0, 1, 2]
        j = [0, 1, 2]
        k = [0, 1, 2]
        for p in i:
            for q in j:
                for r in k:
                    # Skip the center square of the cube
                    if p == q == r == 1:
                        continue
                    # Skip if the point is not "move-to-able"
                    if paths[p, q, r] == 0:
                        continue
                    
                    # Opposite coordinates
                    # This is early evaluation to save computing time if the value gets reused
                    _p = 2 - p
                    _q = 2 - q
                    _r = 2 - r
                    # Coordinate-vs-opposite "sameness" checking
                    # Used to determine whether a position is a center/edge/corner of the cube
                    bp = (p != _p)
                    bq = (q != _q)
                    br = (r != _r)
                    
                    # Subtract opposite point but force the result >= 0
                    tmpPaths[p, q, r] -= min(tmpPaths[p, q, r], paths[_p, _q, _r])
                    
                    # What to add
                    add = 0
                    # The rest of this just manipulates the "add" variable, and, in the end, adds it
                    
                    # Center of a face
                    if bp and not bq and not br:     # 
                        add = np.sum(paths[p, :, :]) # > This section finds all other points
                    elif bq and not bp and not br:   # > on the same face for a given point
                        add = np.sum(paths[:, q, :]) # > at the center of a face
                    elif br and not bp and not bq:   # 
                        add = np.sum(paths[:, :, r]) # 
                    
                    # Edges
                    if bp and bq and not br:                                                  # 
                        add = np.sum(paths[p, q, :])                                          # > This section finds all adjoining
                        add += np.sum(paths[min(1, p):max(1, p)+1, min(1, q):max(1, q)+1, r]) # > centers of faces and corners, as
                    elif bp and br and not bq:                                                # > well as other diagonal edges for
                        add = np.sum(paths[p, :, r])                                          # > a given edge.
                        add += np.sum(paths[min(1, p):max(1, p)+1, q, min(1, r):max(1, r)+1]) # 
                    elif bq and br and not bp:                                                # 
                        add = np.sum(paths[:, q, r])                                          # 
                        add += np.sum(paths[p, min(1, q):max(1, q)+1, min(1, r):max(1, r)+1]) # 
                    
                    # Corners
                    if bp and bq and br:                                                                         # > This section finds all adjacent edges
                        add = np.sum(paths[min(1, p):max(1, p)+1, min(1, q):max(1, q)+1, min(1, r):max(1, r)+1]) # > and centers to a given corner.
                    
                    # Now comes the part where you add (without adding)
                    tmpPaths[p, q, r] = add
        # Ending the iteration
        paths = copy.deepcopy(tmpPaths)
        # Check if max is found
        ismax = (paths==np.amax(paths))
        c = np.count_nonzero(ismax)
        if c == 1:
            final_direction = np.nonzero(ismax)
            break
    # If the max was not found, then the default value will end up getting used
    return final_direction

# ---------------------------- Graphics section ---------------------------- #

# ================================================= # 
class AngleStorage:                                # > This section stores the sin
                                                   # > and cos of the perspective
    sin_th = 0.7717526620201259                    # > angle (the values needed to
    cos_th = -0.6359228165940024                   # > compute the transforms) and
                                                   # > updates it by the position
    sin_ph = 0.9815502530915154                    # > of the mouse
    cos_ph = -0.19120434267030137                  # 
                                                   # 
    theta = 2.26                                   # 
    phi = 1.76318531                               # 
                                                   # 
    two_pi = 2 * math.pi                           # 
                                                   # 
    pressed = False                                # 
    mouseX = 350                                   # 
    mouseY = 350                                   # 
                                                   # 
    def change_perspective(self, evt):             # 
        if not self.pressed:                       # 
            return                                 # 
        self.theta += (evt.y - self.mouseY)/150    # 
        self.theta %= self.two_pi                  # 
        self.phi += (self.mouseX - evt.x)/150      # 
        self.phi %= self.two_pi                    # 
                                                   # 
        self.sin_th = math.sin(self.theta)         # 
        self.cos_th = math.cos(self.theta)         # 
                                                   # 
        self.sin_ph = math.sin(self.phi)           # 
        self.cos_ph = math.cos(self.phi)           # 
                                                   # 
    def press(self, evt):                          # 
        self.pressed = True                        # 
        self.mouseX, self.mouseY = evt.x, evt.y    # 
                                                   # 
    def release(self, evt):                        # 
        self.pressed = False                       # 
                                                   # 
angles = AngleStorage()                            # 
                                                   # 
canvas.bind('<Motion>', angles.change_perspective) # 
canvas.bind('<ButtonPress-1>', angles.press)       # 
canvas.bind('<ButtonRelease-1>', angles.release)   # 
# ================================================ # 

# This is the function that does the main transformation.
def core_graphics(x, y, z, sin_th, cos_th, sin_ph, cos_ph):
    scale = 25
    
    mx = scale * (x - 10)
    my = scale * (y - 10)
    mz = scale * (10 - z)
    
    rx, ry, rz = sin_ph * my + cos_ph * mx, sin_th * mz + cos_th * (cos_ph * my - sin_ph * mx), cos_th * mz + sin_th * (cos_ph * my - sin_ph * mx)
    
    px = 300*rx/(rz + 700) + 450
    py = 300*ry/(rz + 700) + 450
    
    return px, py

# This implements the transformation with multiple data inputs
def integrated_graphics(matrix, bed):
    OVERALL = copy.deepcopy(matrix[:, :, :, 0])
    trig_values = angles.sin_th, angles.cos_th, angles.sin_ph, angles.cos_ph
    
    BED_CONNECTIONS = {}
    for x, row in enumerate(bed):
        for y, point in enumerate(row):
            OVERALL[x, y, point] = 2
            BED_CONNECTIONS[x, y] = []
            with SUPPRESS:
                z_coord = bed[x+1, y]
                BED_CONNECTIONS[x, y].append([x+1, y, z_coord])
            with SUPPRESS:
                z_coord = bed[x, y+1]
                BED_CONNECTIONS[x, y].append([x, y+1, z_coord])
            with SUPPRESS:
                z_coord = bed[x+1, y+1]
                BED_CONNECTIONS[x, y].append([x+1, y+1, z_coord])
    
    OVERALL[[0, 19], :, :][:, [0, 19], :] = 3
    OVERALL[[0, 19], :, :][:, :, [0, 19]] = 3
    OVERALL[:, [0, 19], :][:, :, [0, 19]] = 3
    
    for x, layer in enumerate(OVERALL):
        for y, row in enumerate(layer):
            for z, point in enumerate(row):
                px, py = core_graphics(x, y, z, *trig_values)
                if point == 0:
                    continue
                elif point == 1:
                    canvas.create_oval(px-2, py-2, px+2, py+2, fill="#6495ED", outline="#6495ED")
                elif point == 2:
                    #color = bed_colors[z]
                    #canvas.create_oval(px-2, py-2, px+2, py+2, fill=color, outline=color) # F08080, DC143C
                    vertices = [px, py]
                    for item in BED_CONNECTIONS[x, y]:
                        qx, qy = core_graphics(*item, *trig_values)
                        #canvas.create_line(px, py, qx, qy, width=4, fill=color)
                        vertices.append(qx)
                        vertices.append(qy)
                    with SUPPRESS:
                        canvas.create_polygon(*vertices, fill="#DCDCDC")
                else:
                    color = "#000000"
                    

while True:
    try:
        refresh()
        canvas.delete("all")
        integrated_graphics(POSITIONS, RIVERBED)
        gc.collect()
        tk.update_idletasks()
        tk.update()
    except tkinter.TclError:
        break

for var in vars():
    del var

gc.collect()

