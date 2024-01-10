# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:11:01 2023

@author: lisa_
"""
#------------------------------------------------------------------------------
# Import Modules
#------------------------------------------------------------------------------
import pyvista as pv
import vtk
import numpy as np 
import math
#------------------------------------------------------------------------------
# Set Variables
#------------------------------------------------------------------------------
## Anatomy (create folder in working directory with same name)
name = 'NiTaStraight'

## Input file names containing line data (from vmtk)
f_centreline = name+str('_centreline.vtk')
f_anastamoses = name+str('_anastamoses.vtk')
anat = name+str('.stl')
## Approximate number of segments into which the expander mesh is to be divided vertically
N = 20
## Number of elements in the expander circumference
N_el = 10

## Radius of bender, expander
r_bender = 1
r_expander = 1

## Distance in mm from origin to cutting plane for conduit (factor to multiply the normal)
n_f = 3

## Define datum axes for visualisation 
datum_x = pv.Spline(pv.PolyData([[-100.0,0.0,0.0], [0.0,0.0,0.0], [100.0,0.0,0.0]]).points)
datum_y = pv.Spline(pv.PolyData([(0.0,-100.0,0.0), (0.0,0.0,0.0), (0.0,100.0,0.0)]).points)
datum_z = pv.Spline(pv.PolyData([(0.0,0.0,-100.0), (0.0,0.0,0.0), (0.0,0.0,100.0)]).points)

#------------------------------------------------------------------------------
# Helpful formula 
#------------------------------------------------------------------------------
def Write_XYZ_to_file(filename, data):
    f1 = open(filename,'w')
    for line in data:
        f1.write(str(line[0])+','+str(line[1])+','+str(line[2])+'\n')

def distance(point1, point2):
    return ((point1[0] - point2[0])**2 +(point1[1] - point2[1])**2 +(point1[2] - point2[2])**2)**0.5

d = distance([-15.65435, 56.29486, 14.68991],[-4.43562, 61.57475, 14.29829])


def M1M2_calculator(data):## Vector between adjacent data points
    M1M2 = [] 
    prev_point = [0,0,0]
    i=0
    length=0
    for point in data:
        x = point[0] - prev_point[0]
        y = point[1] - prev_point[1]
        z = point[2] - prev_point[2]
        mag = distance(point, prev_point)
        length+=mag
        M1M2.append([x,y,z,mag,ids[i]])
        i+=1
        prev_point = point 
    length = length - M1M2[0][3]
    M1M2.pop(0)
    return length, M1M2

def Axis_Average(data, ax):
    total = 0
    for point in data:
        total += point[ax]
    avg = total/len(data)
    return avg

def find_closest_points(set1, set2):
    min_distance = float('inf')
    for point1 in set1:
        for point2 in set2:
            dist = distance(point1, point2)
            if dist < min_distance:
                min_distance = dist
                A = point1
                B = point2
    return A, B, min_distance

def CleanPAOpening(openingData, max_radius, side):
    opening = openingData
    current_max = np.inf
    branch=[]
    furthest_point = [0,0,0]
    pointZ=Axis_Average(opening,2)
    while current_max > max_radius:
        centre = (Axis_Average(opening,0), Axis_Average(opening,1), pointZ)
        newOpening = []
        prev_max = 0
        for point in opening:
            radius = distance(point, centre)
            if radius > prev_max:
                furthest_point = point 
                current_max = radius
                prev_max = radius
        for point in opening:
            if point[0] == furthest_point[0] and point[1] == furthest_point[1] and point[2] == furthest_point[2]:
                branch.append(point)
            else:
                newOpening.append(point)
        opening = newOpening
        if abs(abs(Axis_Average(opening,2))-abs(zavg)) > abs(abs(anasT[side][2])-abs(zavg)) or len(opening) > 0.5*len(openingData):
            pointZ = anasT[side][2]
        else:
            pointZ = Axis_Average(opening,2)
    centreO = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
    centreB = (Axis_Average(branch,0), Axis_Average(branch,1), Axis_Average(branch,2))
    return opening, centreO, branch, centreB

def CleanVCOpening(openingData, max_radius, linept):
    opening = openingData
    current_max = np.inf
    branch = []
    furthest_point = [0,0,0]
    while current_max > max_radius:
        if len(opening) < 0.5*len(openingData):
            centre = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
        elif abs(abs(linept[0])-abs(Axis_Average(opening,0))) < 11 and abs(abs(linept[1])-abs(Axis_Average(opening,1))) < 11:
            centre = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
        else:
            centre = linept
        newOpening = []
        prev_max = 0
        for point in opening:
            radius = distance(point, centre)
            if radius > prev_max:
                furthest_point = point 
                current_max = radius
                prev_max = radius
        for point in opening:
            if point[0] == furthest_point[0] and point[1] == furthest_point[1] and point[2] == furthest_point[2]:
                branch.append(point)
            elif abs(point[2]) < abs(centre[2])-15 or abs(point[2]) > abs(centre[2])+15: 
                branch.append(point)
            else:
                newOpening.append(point)
        opening = newOpening
    centreO = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
    centreB = (Axis_Average(branch,0), Axis_Average(branch,1), Axis_Average(branch,2))
    return opening, centreO, branch, centreB

def CleanBranch(openingData, max_radius):
    opening = openingData
    current_max = np.inf
    branch = []
    furthest_point = [0,0,0]
    while current_max > max_radius:
        centre = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
        newOpening = []
        prev_max = 0
        for point in opening:
            radius = distance(point, centre)
            if radius > prev_max:
                furthest_point = point 
                current_max = radius
                prev_max = radius
        for point in opening:
            if point[0] == furthest_point[0] and point[1] == furthest_point[1] and point[2] == furthest_point[2]:
                continue
            else:
                newOpening.append(point)
        opening = newOpening
    centreO = (Axis_Average(opening,0), Axis_Average(opening,1), Axis_Average(opening,2))
    return opening, centreO

def VesselDiameter(pts):
    tang = np.roll(pts, -1, axis=0)[:-1] - pts[:-1]
    n_points = np.size(tang, 0)
    diameters = []
    
    for i in range(n_points):
        # create a plane to cut using centerline data
        plane = vtk.vtkPlane()
        plane.SetOrigin(pts[i][0], pts[i][1], pts[i][2])
        plane.SetNormal(tang[i][0], tang[i][1], tang[i][2])
    
        # create cutter (effectively a slice through the surface - returns line of intersection
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputDataObject(AlignedAnatomy)
        cutter.Update()
       
        # if the plane intersects multiple regions this isolates the largest one
        cnf = vtk.vtkPolyDataConnectivityFilter()
        cnf.SetInputConnection(cutter.GetOutputPort())
        cnf.SetExtractionModeToLargestRegion()
        cnf.Update()
       
        # computes the perimeter derived diameter of the line of intersecti
        integrator = vtk.vtkIntegrateAttributes()
        integrator.SetInputConnection(cnf.GetOutputPort())
        integrator.Update()
        better = pv.wrap(integrator.GetOutput())
        diameters.append(better.get_array('Length')[0]/math.pi)
    
    return diameters
#------------------------------------------------------------------------------
# Actions
#------------------------------------------------------------------------------

##Open Data
centreline = pv.read(f_centreline)
coords = centreline.points
if coords[0][2] > coords[-1][2]:
    coords = coords[::-1]
ids = np.arange(centreline.n_points, dtype=int)

anastamoses = pv.read(f_anastamoses)
anas_coords = anastamoses.points

##Calculate vectors and vector magnitudes between adjacent points

length, M1M2 = M1M2_calculator(coords)

##Segment line and identify key data points

target_length = length/N
dataPoints = [coords[0]]
seg_length = 0
flength=0
flength2 = 0
for point in M1M2: 
    seg_length+= abs(point[3])
    flength2 += point[3]
    if seg_length > target_length: 
        if abs(seg_length-target_length) < abs(last_length-target_length):
            dataPoints.append(coords[point[4]])
            flength+=seg_length
            seg_length = 0
        else:
            dataPoints.append(coords[prev_point[4]])
            flength+=seg_length
            seg_length = point[3]
    last_length= seg_length
    prev_point = point

#------------------------------------------------------------------------------

## Find z Axis

datamean = np.array(dataPoints).mean(axis=0)

vector = [dataPoints[0][0]-dataPoints[-1][0], dataPoints[0][1]-dataPoints[-1][1], dataPoints[0][2]-dataPoints[-1][2]]
disp = distance(dataPoints[0], dataPoints[-1])
normal = (vector[0]/disp, vector[1]/disp,vector[2]/disp)

start = distance(dataPoints[0], [0,0,0])
end = distance(dataPoints[1], [0,0,0])
linepts = normal*np.mgrid[-start:end:N*1j][:, np.newaxis]
linepts += datamean 

trueZ = []
interval = (start+end)/N
for n in range(0,N):
    point = [0,0,-start+n*interval]
    trueZ.append(point)

lineptsPD = pv.PolyData(linepts)
trueZPD = pv.PolyData(trueZ)

temp, translation_matrix = lineptsPD.align(trueZPD, return_matrix=True)
dataOnZ = pv.PolyData(dataPoints).transform(translation_matrix)

anasT = pv.PolyData(anas_coords).transform(translation_matrix)
anasT_clip = anasT.clip((1,0,0),(15,0,0))
anasT_clip = anasT_clip.clip((-1,0,0),(-15,0,0))

anasT = anasT.points

anas_clip = pv.PolyData(anas_coords).clip((1,0,0),(15,0,0))
anas_clip = anas_clip.clip((-1,0,0),(-15,0,0))

## Open anatomy and extract vessel openings
a = pv.read(anat)

p = pv.Plotter()
p.add_mesh(a, opacity=0.5)
p.add_mesh(pv.PolyData(dataPoints), color='red')
p.add_mesh(anas_clip, color='blue')
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
#p.view_xz()
#p.save_graphic('Graphics/NiTa_raw.svg')
p.show()

a.transform(translation_matrix)

p = pv.Plotter()
p.add_mesh(a, opacity = 0.5)
p.add_mesh(dataOnZ, color='red')
p.add_mesh(anasT_clip, color='blue')
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
#p.view_xz()
#p.save_graphic('Graphics/NiTa_translate.svg')
p.show()


anasT_clip = anasT_clip.points
zavg = Axis_Average(anasT_clip,2)

Edges = a.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)
EdgePoints = Edges.points

## Calculate X Axis and transform centreline data

currentXAxis = [(anasT_clip[0][0], anasT_clip[0][1], zavg), (anasT_clip[-1][0], anasT_clip[-1][1], zavg)]
targetXAxis = [(anasT_clip[0][0],0,0),(anasT_clip[-1][0],0,0)]

NewXAxis, rotation_matrix = pv.PolyData(currentXAxis).align(pv.PolyData(targetXAxis), return_matrix=True)

rotation_matrix[0,3]=0
rotation_matrix[1,3]=0

alignedDataSet = pv.PolyData(dataOnZ).transform(rotation_matrix)

## Transform anatomy 
AlignedAnatomy = a.transform(rotation_matrix)
anasT = pv.PolyData(anasT).transform(rotation_matrix).points

centrelineZ0, origin, Origin_to_Centreline = find_closest_points(alignedDataSet.points, [(0,0,0)])

align_coords = pv.PolyData(coords).transform(translation_matrix)
cl = align_coords.transform(rotation_matrix)

anasT_clip = pv.PolyData(anasT).clip((1,0,0),(15,0,0)).clip((-1,0,0),(-15,0,0))

Zrotate = 0
if (anasT[0][1] + anasT[-1][1])/2 < centrelineZ0[1]:
    alignedDataSet.rotate_vector((0, 0, 1), 180, inplace=True)
    AlignedAnatomy.rotate_vector((0, 0, 1), 180, inplace=True)
    Zrotate = 1
    cl = cl.rotate_vector((0, 0, 1), 180, inplace=True)
    anasT_clip = anasT_clip.rotate_vector((0, 0, 1), 180, inplace=True)

## Save fully aligned anatomy
AlignedAnatomy.save(name+str('/AlignedAnatomy.stl'))
alignedPoints = alignedDataSet.points

##--------------------------------------------------------------------------------------------------------------------

## Extract Edges of Anatomy for pinned BC 
edges = AlignedAnatomy.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)

##Split into quadrants
topOpening = []
baseOpening = []
leftOpening = []
rightOpening = []
for point in edges.points: 
    if point[2] > 30:
        topOpening.append(point)
    elif point[2] < -30: 
        baseOpening.append(point)
    if  -30 < point[2] < 30 :
        if point[0] < 0:
            rightOpening.append(point)
        else: 
            leftOpening.append(point)
            
## Identify top and bottom of VC centreline

top = -1
bottom = 0

##Find right and left side of anastamoses lines
if anasT[0][0] < anasT[-1][0]:
    leftSide=-1
    rightSide=0
else:
    leftSide=0
    rightSide=-1

if Zrotate == 1 :
    if leftSide == 0:
        leftSide = -1
        rightSide = 0
    else:
        leftSide = 0 
        rightSide = -1

clean_top, centreT, branchT, centreBT = CleanVCOpening(topOpening, 15, alignedPoints[top])
clean_base, centreB, branchB, centreBB = CleanVCOpening(baseOpening, 15, alignedPoints[bottom])

clean_leftOpening, centreL, branchL, centreBL = CleanPAOpening(leftOpening, 15, leftSide)
clean_rightOpening, centreR, branchR, centreBR = CleanPAOpening(rightOpening, 15, rightSide)

branchR, centreBR = CleanBranch(branchR, 15)
branchT, centreBT = CleanBranch(branchT, 15)

p = pv.Plotter()
p.add_mesh(AlignedAnatomy, opacity=0.5)
p.add_mesh(pv.PolyData(anasT_clip), color='blue')
p.add_mesh(pv.PolyData(alignedPoints), color='red')
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
p.view_isometric(negative=True)
#p.save_graphic('Graphics/NiTa_Openings.svg')
p.show()


p = pv.Plotter()
p.add_mesh(AlignedAnatomy, opacity=0.5)
p.add_mesh(pv.PolyData(clean_top), color='yellow')
p.add_mesh(pv.PolyData([centreT]), color='yellow')
p.add_mesh(pv.PolyData(clean_base), color='r')
p.add_mesh(pv.PolyData([centreB]), color='r')
p.add_mesh(pv.PolyData(clean_leftOpening), color='black')
p.add_mesh(pv.PolyData([centreL]), color='black')
p.add_mesh(pv.PolyData(clean_rightOpening), color='orange')
p.add_mesh(pv.PolyData([centreR]), color='orange')
p.add_mesh(pv.PolyData(branchT), color='green')
p.add_mesh(pv.PolyData([centreBT]), color='green')
p.add_mesh(pv.PolyData([centreBR]), color='pink')
p.add_mesh(pv.PolyData(branchR), color='pink')
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
p.view_xz(negative=True)
#p.save_graphic('Graphics/NiTa_Openings.svg')
p.show()

Write_XYZ_to_file(name+str('/VesselTop.txt'),clean_top)
Write_XYZ_to_file(name+str('/VesselBase.txt'),clean_base)
Write_XYZ_to_file(name+str('/VesselTopBranch.txt'),branchT)
Write_XYZ_to_file(name+str('/VesselLeft.txt'),clean_leftOpening)
Write_XYZ_to_file(name+str('/VesselRight.txt'),clean_rightOpening)
Write_XYZ_to_file(name+str('/VesselRightBranch.txt'),branchR)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#Find Vessel Radius for expansion

IVCpts = []
SVCpts = []

for point in alignedPoints:
    if 15 < point[2] < 25:
        SVCpts.append(point)
    elif point[2] < -17:
        IVCpts.append(point)

IVC_dia = VesselDiameter(IVCpts)
SVC_dia = VesselDiameter(SVCpts)

IVC_avg_dia = np.average(IVC_dia)
SVC_avg_dia = np.average(SVC_dia)

Write_XYZ_to_file(name+str('/VesselDiameter.txt'),[(SVC_avg_dia,IVC_avg_dia,0)])
#------------------------------------------------------------------------------

##Find Straight Expander Axis
ZAxis = []

up_count = 3
down_count = 3
for point in alignedPoints:
    if point[2] == centrelineZ0[2]:
        continue
    elif point[2] > centrelineZ0[2]:
        up_count +=1
    elif point[2] < centrelineZ0[2]:
        down_count +=1

for n in range(1,down_count+1):
    x = down_count + 1 - n
    point = [0,0,centrelineZ0[2]-x*target_length]
    ZAxis.append(point)

for x in range(0,up_count+1):
    point = [0,0,centrelineZ0[2]+x*target_length]
    ZAxis.append(point)


#------------------------------------------------------------------------------
## Create Straightener Data
length, BigM1M2 = M1M2_calculator(alignedPoints)

end1 = alignedPoints[0]
end2 = alignedPoints[-1]

np1 = (end1[0]-BigM1M2[0][0], end1[1]-BigM1M2[0][1], end1[2]-BigM1M2[0][2])
np2 = (end1[0]-2*BigM1M2[0][0], end1[1]-2*BigM1M2[0][1], end1[2]-2*BigM1M2[0][2])
np3 = (end1[0]-3*BigM1M2[0][0], end1[1]-3*BigM1M2[0][1], end1[2]-3*BigM1M2[0][2])
bottom = np.array([np3,np2, np1])

np4 = (end2[0]+BigM1M2[-1][0], end2[1]+BigM1M2[-1][1], end2[2]+BigM1M2[-1][2])
np5 = (end2[0]+2*BigM1M2[-1][0], end2[1]+2*BigM1M2[-1][1], end2[2]+2*BigM1M2[-1][2])
np6 = (end2[0]+3*BigM1M2[-1][0], end2[1]+3*BigM1M2[-1][1], end2[2]+3*BigM1M2[-1][2])
top = np.array([np4,np5,np6])

straightenerPoints = np.vstack([bottom, alignedPoints, top])

##write straightener to file
spline = pv.Spline(straightenerPoints)
bender = spline.tube(radius=r_bender, n_sides = N_el)

bender.save(name+str('/Bender.stl'))
Write_XYZ_to_file(name+str('/BenderPoints.txt'),straightenerPoints)


##Find Displacements for Bender Set Up
Displacement_Vectors = []
for lineP, axisP in zip(straightenerPoints,ZAxis):
    vector = [axisP[0]-lineP[0],axisP[1]-lineP[1],axisP[2]-lineP[2]]
    Displacement_Vectors.append(vector)

Write_XYZ_to_file(name+str('/BenderDisplacement.txt'),Displacement_Vectors)

##Find Rotations for Bender Set Up 
s_length, s_BigM1M2 = M1M2_calculator(straightenerPoints)
Rotations = []
for v in s_BigM1M2:
    yrot = np.arcsin(v[0]/(v[0]**2+v[1]**2+v[2]**2)**0.5)
    xrot = np.arctan(-v[1]/v[2])
    zrot = 0
    Rotations.append([-xrot,-yrot,zrot])

Rotations.pop(0)
Rotations.append(Rotations[-1])

Write_XYZ_to_file(name+str('/BenderRotations.txt'),Rotations)

#------------------------------------------------------------------------------
##Save conduit line section 

ConduitZ = -IVC_avg_dia
ConduitPoints = []
align_coords = pv.PolyData(coords).transform(translation_matrix)
align_coords = align_coords.transform(rotation_matrix).points

for point in straightenerPoints: 
    if point[2] < ConduitZ:
        ConduitPoints.append(point)

for point in align_coords:
    if ConduitPoints[-1][2] < point[2] < ConduitZ:
        extra_point = point

ConduitPoints.append(extra_point)

Write_XYZ_to_file(name+str('/ConduitPoints.txt'),ConduitPoints)

#----------------------------------------------------------------------------
##Save cutting plane

datamean = np.array(anasT_clip.points).mean(axis=0)
uu, dd, vv = np.linalg.svd(anasT_clip.points - datamean)

start = 15
end = 15
centre_pa = vv[0]*np.mgrid[-start:end:N*1j][:, np.newaxis]
centre_pa += datamean 

theta = np.arctan(-vv[0][2]/vv[0][0])

centre_cv = cl.clip((0,0,1),(0,0,0))
centre_cv_p = centre_cv.clip((0,0,-1),(0,0, -IVC_avg_dia))

datamean = np.array(centre_cv_p.points).mean(axis=0)
uu, dd, vv = np.linalg.svd(centre_cv_p.points - datamean)

start = -15
end = 15
centre_cv = vv[0]*np.mgrid[start:end:8*1j][:, np.newaxis]
centre_cv += datamean 

vv[0][0] = vv[0][2]*np.tan(theta)

start = -15
end = 15
centre_cv_proj = vv[0]*np.mgrid[start:end:8*1j][:, np.newaxis]
centre_cv_proj += datamean 

if vv[0][2] < 0:
    vv[0][0] *= -1
    vv[0][1] *= -1
    vv[0][2] *= -1

centre = (-n_f*vv[0][0], -n_f*vv[0][1], -n_f*vv[0][2])
cut_plane = pv.Plane(center=centre, direction=vv[0], i_size=50, j_size=50, i_resolution=10, j_resolution=10)

Write_XYZ_to_file(name+str('/CutPlane.txt'),[centre,vv[0]])

p = pv.Plotter()
p.add_mesh(pv.Spline(centre_cv), color='r', line_width = 4)
p.add_mesh(pv.Spline(centre_pa), color='b', line_width = 4)
p.add_mesh(pv.Spline(centre_cv_proj), color='green', line_width = 4)
p.add_mesh(cut_plane, color='orange', opacity = 0.5)
p.add_mesh(AlignedAnatomy, opacity = 0.5)
p.add_mesh(pv.Spline(anasT_clip.points), color = 'green', line_width = 2)
p.add_mesh(pv.Spline(centre_cv_p.points), color='green', line_width = 2)
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
p.view_xz()
p.show()


#------------------------------------------------------------------------------
##Create Expander Part
ExpanderPoints = []
for point in ZAxis: 
    if point[2] < 48 and point[2] > -55:
        ExpanderPoints.append(point)

Write_XYZ_to_file(name+str('/ExpanderPoints.txt'),ExpanderPoints)

espline = pv.Spline(pv.PolyData(ExpanderPoints).points)
expander = espline.tube(radius = r_expander, n_sides = N_el)

expander.save(name+str('/Expander.stl'))

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##Visualise Data

p = pv.Plotter()

p.add_mesh(bender, color='r', opacity=0.5)
p.add_mesh(expander, color='b')
p.add_mesh(AlignedAnatomy, opacity=0.5)
p.add_mesh(datum_x, color='black')
p.add_mesh(datum_y, color='black')
p.add_mesh(datum_z, color='black')
p.view_xz()
p.show()


#==============================================================================



