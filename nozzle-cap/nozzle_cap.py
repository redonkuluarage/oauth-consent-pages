"""Parametric yellow nozzle cap (first draft, dims estimated from photos).
Printing orientation: rim down on the bed, tabbed end facing up.
Run: python3 nozzle_cap.py  -> nozzle_cap.stl
All values in mm. EDIT THE PARAMETERS, re-run, re-slice."""
import math
import numpy as np
from manifold3d import Manifold, CrossSection

# ---- parameters (ESTIMATES: verify with calipers/Photo Caliper) ----
H          = 57.0   # overall height, rim to tabbed end face
WALL       = 1.6    # side wall thickness
END_T      = 3.0    # end plate thickness
RIM_OD     = 44.0   # flange OD at open end
RIM_H      = 2.5    # flange height
RIM_T      = 2.0    # flange radial thickness
LOW_OD0    = 41.0   # lower section OD just above flange
STEP_Z     = 31.0   # height (from rim) of the step line
LOW_OD1    = 39.0   # lower section OD at step (top of lower section)
UP_OD0     = 38.0   # upper section OD just above step
UP_OD1     = 35.0   # upper section OD at tabbed end
CLR        = 0.0    # extra radial clearance added to inner features if it's a mating part
CENTER_D   = 8.0
# small holes: (angle deg, radius from center, diameter)
HOLES = [(100,11.0,3.2),(62,11.5,2.2),(20,9.5,2.4),(-25,9.0,2.0),
         (-70,10.0,2.2),(-110,10.5,3.0),(150,11.0,2.0)]
DIMPLES = [(180,11.5,8.0),(0,11.5,8.0)]   # blind circles, 0.4 deep
TAB_ANG = [45,135,225,315]; TAB_R=14.0; TAB_L=8.0; TAB_W=3.0; TAB_H=1.5
RIB_ANG = [-90, 90]; RIB_SPAN=40; RIB_Z=9.0; RIB_T=1.2; RIB_PROJ=1.6; RIB_H=1.4
PLAT_OD = 26.0; PLAT_T=1.2; PLAT_H=6.0     # inner stepped ring under end plate
SEG = 96

def revolve(pts):
    return Manifold.revolve(CrossSection([pts]), SEG)

def cyl(d, h, z=0, x=0, y=0):
    return Manifold.cylinder(h, d/2, d/2, SEG).translate([x, y, z])

# outer shell profile (r,z), solid
r = lambda od: od/2
outer = [(0,0),(r(RIM_OD),0),(r(RIM_OD),RIM_H),(r(LOW_OD0),RIM_H),(r(LOW_OD1),STEP_Z),
         (r(UP_OD0),STEP_Z),(r(UP_OD1),H),(0,H)]
body = revolve(outer)
# hollow bore: wall thickness WALL, leaving END_T plate
inner = [(0,-1),(r(RIM_OD)-RIM_T,-1),(r(RIM_OD)-RIM_T,RIM_H),(r(LOW_OD0)-WALL,RIM_H),
         (r(LOW_OD1)-WALL,STEP_Z),(r(UP_OD0)-WALL,STEP_Z),(r(UP_OD1)-WALL,H-END_T),(0,H-END_T)]
body = body - revolve(inner)

# inner stepped platform
plat = cyl(PLAT_OD, PLAT_H+END_T, H-END_T-PLAT_H) - cyl(PLAT_OD-2*PLAT_T, PLAT_H+END_T+1, H-END_T-PLAT_H-0.5)
body = body + plat

# raised tabs on end face
for a in TAB_ANG:
    t = Manifold.cube([TAB_L, TAB_W, TAB_H+0.5], True).translate([0,0,(TAB_H+0.5)/2])
    t = t.rotate([0,0,0]).translate([TAB_R,0,0]).rotate([0,0,a+90-90])
    # orient long side tangential: cube long axis is X; after translating to +X radial, rotate local 90
    t = Manifold.cube([TAB_W, TAB_L, TAB_H+0.5], True).translate([TAB_R,0,(TAB_H+0.5)/2]).rotate([0,0,a])
    body = body + t.translate([0,0,H-0.5])

# inner snap ribs (arcs on the inner wall)
rin = r(LOW_OD0)-WALL - (RIM_Z:=0)*0
for a in RIB_ANG:
    ro = r(LOW_OD0) - WALL + 0.3  # embed into wall
    ri = ro - RIB_PROJ - 0.3
    arc = (cyl(2*ro, RIB_H, RIB_Z) - cyl(2*ri, RIB_H+1, RIB_Z-0.5))
    wedge = Manifold.cube([60,60,RIB_H+2],False).translate([0,-30,RIB_Z-1])
    # keep only a RIB_SPAN-degree sector
    sec = CrossSection([[(0,0)]+[(60*math.cos(math.radians(RIB_SPAN/2*s)),60*math.sin(math.radians(RIB_SPAN/2*s))) for s in np.linspace(-1,1,12)]]).extrude(RIB_H+2).translate([0,0,RIB_Z-1])
    body = body + (arc ^ sec.rotate([0,0,a]))

# holes through end plate
for ang,rad,d in HOLES:
    x,y = rad*math.cos(math.radians(ang)), rad*math.sin(math.radians(ang))
    body = body - cyl(d, END_T+PLAT_H+2, H-END_T-PLAT_H-1, x, y)
body = body - cyl(CENTER_D, END_T+PLAT_H+2, H-END_T-PLAT_H-1)
for ang,rad,d in DIMPLES:
    x,y = rad*math.cos(math.radians(ang)), rad*math.sin(math.radians(ang))
    body = body - cyl(d, 1, H-0.4, x, y)

mesh = body.to_mesh()
import trimesh
m = trimesh.Trimesh(mesh.vert_properties[:,:3], mesh.tri_verts)
print("watertight:", m.is_watertight, "bbox mm:", m.extents.round(1), "volume cm3:", round(m.volume/1000,2))
m.export(__file__.replace("nozzle_cap.py","nozzle_cap.stl"))
