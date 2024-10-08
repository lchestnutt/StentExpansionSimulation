# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 12:43:58 2023

@author: 2459326C
"""

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#imports
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

import sys

from abaqus import *
from abaqusConstants import *
sys.path.insert(8, 
    r'c:/SIMULIA/EstProducts/2020/win_b64/code/python2.7/lib/abaqus_plugins/stlImport')
import stl2inp
import mesh
import __main__
import section 
import regionToolset 
import part
import material 
import assembly 
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#variables
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
## Name of the analysis folder
name = 'ViPuStraight'

myE = 'Expander'
myString = name+str('-Model')
myW = 'Wire'
myS = 'Stent'
myV = 'Vessel'

##File paths to STLs - add abaqus working directory
ExpanderPath = str('C:/AbaqusTemp/')+name+str('/Expander.stl')
BenderPath = str('C:/AbaqusTemp/')+name+str('/Bender.stl')
StentPath = 'C:/AbaqusTemp/ExportedStent.sat'
AnatomyPath = str('C:/AbaqusTemp/')+name+str('/AlignedAnatomy.stl')

## Alter Position of stent from straightening axis
stent_up = 0 ## move stent along Z axis (mm)
stent_rotate = 0 ## rotate stent around Z axis (degrees)

StentXPath = ((5*np.sin(stent_rotate*np.pi/180),5*np.cos(stent_rotate*np.pi/180),36.5+stent_up),(-5*np.sin(stent_rotate*np.pi/180),-5*np.cos(stent_rotate*np.pi/180),36.5+stent_up)) ## Measured from the input file, assuming that the same stent geometry is imported

## Radius of bender, expander (same as preprocessing script)
r_bender = 1
r_expander = 1

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#functions
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

##General

def Import_Stl_file(model,part,path):
    stl2inp.STL2inp(stlfile=path, modelName='temp', mergeNodesTolerance=1E-06)
    mdb.models['temp'].parts.changeKey(fromName='PART-1', toName=part)
    mdb.models[model].Part(part, mdb.models['temp'].parts[part])
    del mdb.models['temp']

def Import_ACIS_file(model, part, path):
    acis = mdb.openAcis(path, scaleFromFile=OFF)
    mdb.models[model].PartFromGeometryFile(name=part, 
        geometryFile=acis, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)

def Read_XYZ_from_file(filename):
    f = open(filename, 'r')
    Data = []
    for point in f:
        p = point.split(',')
        coord = [float(p[0]),float(p[1]),float(p[2])]
        Data.append(coord)
    return Data

def Axis_Average(data, ax):
    total = 0
    for point in data:
        total += point[ax]
    avg = total/len(data)
    return avg

#------------------------------------------------------------------------------
## Datums

def Create_Ref_Point(model,point):
    a = mdb.models[model].rootAssembly
    RP = a.ReferencePoint(point=point)
    r = a.referencePoints

def Create_Set_All_Elements(model,part,set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    p.Set(elements=e, name=set_name)

def Create_Set_All_Cells(model,part,set_name):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    p.Set(cells=c, name=set_name)

def Create_Set_Elements_Sphere(model,part, centre, radius, set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    elements = e.getByBoundingSphere(centre, radius)
    p.Set(elements=elements, name=set_name)

def Create_Set_Node_Ring(model, part, centre, radius, set_name):
    p = mdb.models[model].parts[part]
    n = p.nodes
    nodes = n.getByBoundingSphere(centre, radius)
    p.Set(nodes=nodes, name=set_name)

def Create_Set_Node_Rings(model, part, ring_indices, radius, set_name):
    p = mdb.models[model].parts[part]
    n = p.nodes
    Sets = []
    for i in ring_indices:
        nodes = n.getByBoundingSphere(ExpanderPoints[i], radius=radius)
        nodeSet = p.Set(nodes=nodes, name='temp'+str(i))
        Sets.append(nodeSet)
    p.SetByMerge(name=set_name, sets=Sets[:])
    for i in ring_indices:
        del mdb.models[model].parts[part].sets['temp'+str(i)]

def Create_Set_Node_Ring_Cyl(model, part, centre, radius, set_name):
    p = mdb.models[model].parts[part]
    n = p.nodes
    nodes = n.getByBoundingCylinder((centre[0], centre[1], centre[2]+0.1), (centre[0], centre[1], centre[2]-0.1), radius)
    p.Set(nodes=nodes, name=set_name)

def Create_Set_Nodes(model, part, node_list, set_name):
    p = mdb.models[model].parts[part]
    n = p.nodes
    Sets = []
    for point,x in zip(node_list, range(0,len(node_list))):
        node = n.getByBoundingSphere(point, radius=1e-6)
        nodeSet = p.Set(nodes=(node), name='temp'+str(x))
        Sets.append(nodeSet)
    p.SetByMerge(name=set_name, sets=Sets[:])
    for n in range(0,len(node_list)):
        del mdb.models[model].parts[part].sets['temp'+str(n)]

def Create_Set_Cyl_Elements(model, part, points, radius, set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    Sets = []
    for n in range(0,len(points)-3):
        elements = e.getByBoundingCylinder(points[n], points[n+3], radius=radius+8)
        elementSet = p.Set(elements=(elements), name='temp'+str(n))
        Sets.append(elementSet)
    p.SetByMerge(name=set_name, sets=Sets[:])
    for n in range(0,len(points)-3):
        del mdb.models[model].parts[part].sets['temp'+str(n)]

def Create_Set_Conduit(model, part, points, diameter, normal, start, set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    Sets = []
    end = (start[0]-2*diameter*normal[0], start[1]-2*diameter*normal[1], start[2]-2*diameter*normal[2])
    elements = e.getByBoundingCylinder(start, end, radius=diameter*1.2)
    elementSet = p.Set(elements=(elements), name='temp'+str(0))
    Sets.append(elementSet)
    for n in range(0,len(points)-3):
        elements = e.getByBoundingCylinder(points[n], points[n+3], radius=diameter)
        elementSet = p.Set(elements=(elements), name='temp'+str(n+1))
        Sets.append(elementSet)
    p.SetByMerge(name=set_name, sets=Sets[:])
    for n in range(0,len(points)-3):
        del mdb.models[model].parts[part].sets['temp'+str(n)]

def Create_Set_Boolean(model, part, full_set, cut_set, set_name):
    p = mdb.models[model].parts[part]
    p.SetByBoolean(name=set_name, operation=DIFFERENCE, sets=(p.sets[full_set], p.sets[cut_set]))

def Delete_Element_Sets(model,part,sets):
    p = mdb.models[model].parts[part]
    for group in sets:
        p.deleteElement(elements=[p.sets[group]])

def Create_Datum_Cylind_CSYS_Part(model, part, origin, R, RT, csys_name):
    p = mdb.models[model].parts[part]
    csys = p.DatumCsysByThreePoints(name=csys_name, coordSysType=CYLINDRICAL, origin=origin, point1=R, point2=RT)
    i = p.features[csys_name].id
    return i 

def Create_Datum_Cylind_CSYS(model, origin, R, RT, csys_name):
    a = mdb.models[model].rootAssembly
    csys = a.DatumCsysByThreePoints(name=csys_name, origin=origin, point2=R, coordSysType=CYLINDRICAL, point1=RT)
    i = a.features[csys_name].id
    return i 

def Create_Datum_Axis_Two_Points(model, part, points):
    p = mdb.models[model].parts[part]
    axis = p.DatumAxisByTwoPoint(point1=points[0], point2=points[1])
    i = axis.id
    return i

def Datum_Axis_From_Principle(model, part, axis):
    p = mdb.models[model].parts[part]
    axis = p.DatumAxisByPrincipalAxis(principalAxis=axis)
    i = axis.id
    return i

def Create_Surface_Orphan_Mesh(model, part, side):
    p = mdb.models[model].parts[part]
    if side =='internal':
        side2Elements = p.sets['ALL-Elements'].elements
        p.Surface(side2Elements=side2Elements, name='Internal')
    if side == 'external':
        try:
            side1Elements = p.sets['ALL-Elements'].elements
            p.Surface(side1Elements=side1Elements, name='External')
        except:
            s = p.faces
            side1Faces = s.getSequenceFromMask(mask=('[#ffffffff:54 #3fffffff ]', ), )
            p.Surface(side1Faces=side1Faces, name='External')

def Create_Surface_Stent_Outer(model, name):
    p = mdb.models[model].parts['Stent']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#ffffffff:16 #fffff #0:2 #2 ]', ), )
    p.Surface(side1Faces=side1Faces, name=name)


#------------------------------------------------------------------------------
##Assembly and section Assignment 
def Instance_to_Assembly(model, part, instance):
    a = mdb.models[model].rootAssembly
    p = mdb.models[model].parts[part]
    a.Instance(name=instance, part=p, dependent=ON)
    a.regenerate()

def Rotate_Instance(model, instance, start, direction, angle):
    a = mdb.models[model].rootAssembly
    a.rotate(instanceList=(instance, ), axisPoint=start, axisDirection=direction, angle=angle)

def Edge_Edge_Constraint(model, instance_move, instance_fix, move_axis, fix_axis, Flip_ON_OFF):
    a = mdb.models[model].rootAssembly
    d1 = a.instances[instance_move].datums
    d2 = a.instances[instance_fix].datums
    a.EdgeToEdge(movableAxis=d1[move_axis], fixedAxis=d2[fix_axis], flip=Flip_ON_OFF)

def Create_Material(model, material_name, D, YM, PR, PlasDat, Damp):
    mdb.models[model].Material(name=material_name)
    if D != 'None':
        mdb.models[model].materials[material_name].Density(table=((D, ), ))
    if YM != 'None':
        mdb.models[model].materials[material_name].Elastic(table=((YM, PR), ))
    if PlasDat != 'None':
        mdb.models[model].materials[material_name].Plastic(table=PlasDat)
    if Damp != 'None':
        mdb.models[model].materials[material_name].Damping(alpha=Damp)

def Create_Material_Yeoh(model, material_name, density, table_data, Damp):
    mdb.models[model].Material(name=material_name)
    mdb.models[model].materials[material_name].Density(table=((density,), ))
    mdb.models[model].materials[material_name].Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=YEOH, volumetricResponse=VOLUMETRIC_DATA, table=(table_data, ))
    mdb.models[model].materials[material_name].Damping(alpha=Damp)

def Create_Material_Yeoh_TD(model, material_name, density, table_data, Damp):
    mdb.models[model].Material(name=material_name)
    mdb.models[model].materials[material_name].Density(table=((density,), ))
    mdb.models[model].materials[material_name].Hyperelastic(materialType=ISOTROPIC, testData=ON, type=YEOH, volumetricResponse=VOLUMETRIC_DATA, table=())
    mdb.models[model].materials[material_name].hyperelastic.UniaxialTestData(table=table_data)
    mdb.models[model].materials[material_name].Damping(alpha=Damp)

def Create_Surface_Section(Model, section_name, Density):
    mdb.models[Model].SurfaceSection(name=section_name, useDensity=ON, 
        density=Density)

def Create_Solid_Homogeneous_Section(model,section_name, material_name):
    mdb.models[model].HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

def Create_Shell_Section(model, section_name, thickness, material_name):
    mdb.models[model].HomogeneousShellSection(name=section_name, preIntegrate=OFF,material=material_name, thicknessType=UNIFORM, thickness=thickness, thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)

def Create_Membrane_Section(model, section_name, thickness, material_name):
    mdb.models[model].MembraneSection(name=section_name, material=material_name, thicknessType=UNIFORM, thickness=thickness, thicknessField='', poissonDefinition=DEFAULT)

def Section_Assignment(model, part, set_name, section_name):
    p = mdb.models[model].parts[part]
    region = p.sets[set_name]
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

#------------------------------------------------------------------------------
## BCs and Constraints
def Create_Dynamic_Displacement_Step(model, step_name, prev_step_name, timePeriod, NL_ON_OFF):
    mdb.models[model].ExplicitDynamicsStep(name=step_name, previous=prev_step_name, timePeriod=timePeriod, improvedDtMethod=ON, nlgeom=NL_ON_OFF)

def Create_Smooth_Step_Amplitude(model, amp_name, timeSpan, data_points):
    mdb.models[model].SmoothStepAmplitude(name=amp_name, timeSpan=timeSpan, data=data_points)

def Apply_Displacement_BC_Expand(model, instance, set_name, BC_name, amp_name, CSYS_i, step, u1, u2, u3, ur1, ur2, ur3, deactivate_step):
    a = mdb.models[model].rootAssembly
    region = a.instances[instance].sets[set_name]
    csys = a.datums[CSYS_i]
    mdb.models[model].DisplacementBC(name=BC_name, createStepName=step, region=region, u1=u1, u2=u2, u3=u3, ur1=ur1, ur2=ur2, ur3=ur3, amplitude=amp_name, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=csys)
    if deactivate_step != 'None':
        mdb.models[model].boundaryConditions[BC_name].deactivate(deactivate_step)

def Apply_Displacement_BC_Bend(model, centre, BC_name, step_name, next_step, amp_name, u1, u2, u3, ur1, ur2, ur3):
    a = mdb.models[model].rootAssembly
    r = a.referencePoints
    RPref = r.findAt(centre,)
    region1=regionToolset.Region(referencePoints=(RPref,))
    mdb.models[model].DisplacementBC(name=BC_name, createStepName=step_name, region=region1, u1=u1, u2=u2, u3=u3, ur1=ur1, ur2=ur2, ur3=ur3, amplitude=amp_name, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    mdb.models[model].boundaryConditions[BC_name].deactivate(next_step)

def Apply_Displacement_BC_Pin(model, instance, set_name, CSYS_i, BC_name, step, amp_name, u1, u2, u3, ur1, ur2, ur3):
    a = mdb.models[model].rootAssembly
    region = a.instances[instance].sets[set_name]
    csys = a.datums[CSYS_i]
    mdb.models[model].DisplacementBC(name=BC_name, createStepName=step, region=region, u1=u1, u2=u2, u3=u3, ur1=ur1, ur2=ur2, ur3=ur3, amplitude=amp_name, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=csys)

def Rigid_Ring_Constraint(model, centre, instance, RingSet, constraint_name):
    a = mdb.models[model].rootAssembly
    r = a.referencePoints
    RPref = r.findAt(centre,)
    region1=regionToolset.Region(referencePoints=(RPref,))
    PinRing = a.instances[instance].sets[RingSet]
    mdb.models[model].RigidBody(name=constraint_name, refPointRegion=region1, pinRegion=PinRing)

def Pinned_BC(model, instance, set_name, BC_name, step):
    a = mdb.models[model].rootAssembly
    region = a.instances[instance].sets[set_name]
    mdb.models[model].PinnedBC(name=BC_name, createStepName=step, region=region, localCsys=None)

def ZSymmetry_BC(model, instance, set_name, CSYS_i, BC_name, step):
    a = mdb.models[model].rootAssembly
    region = a.instances[instance].sets[set_name]
    csys = a.datums[CSYS_i]
    mdb.models[model].ZsymmBC(name=BC_name, createStepName=step, region=region, localCsys=csys)

def Create_Symmetry_BC(file, vessel):
    Nodes = Read_XYZ_from_file(name+str(file))
    Centre = (Axis_Average(Nodes,0), Axis_Average(Nodes,1), Axis_Average(Nodes,2))
    Create_Set_Nodes(myString, myV, Nodes, vessel+str('Opening'))
    i = Create_Datum_Cylind_CSYS(myString, Centre, Nodes[0], Nodes[5], str('Cylind_')+vessel)
    ZSymmetry_BC(myString, myV+str('-1'), vessel+str('Opening'), i, vessel+str('-Symmetry'), 'Initial')


#------------------------------------------------------------------------------
## Element Types and Mesh 

def Set_SFM3D3_Elem_Type(model, part, set_name):
    elemType1 = mesh.ElemType(elemCode=SFM3D3, elemLibrary=EXPLICIT)
    p = mdb.models[model].parts[part]
    pickedRegions = p.sets[set_name]
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))

def Set_C3D8R_Elem_Type(model, part, set_name):
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
    p = mdb.models[model].parts[part]
    pickedRegions = p.sets[set_name]
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

def Set_S3R_Elem_Type(model, part, set_name):
    elemType1 = mesh.ElemType(elemCode=S3R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF)
    p = mdb.models[model].parts[part]
    pickedRegions = p.sets[set_name]
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))

def Mesh_Part(model, part, seed_size):
    p = mdb.models[model].parts[part]
    p.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

def Mesh_Stent(model, bias_min, bias_max, globe):
    p = mdb.models[model].parts['Stent']
    e = p.edges
    ## Edge selection extracted from Abaqus RPY file after meshing manually
    pickedEndEdges = e.getSequenceFromMask(mask=(
        '[#14000 #b4740000 #300000a3 #a0000000 #16d03 #d16d0000 #16800016', 
        ' #1450000 #a0000000 #a8001685 #68000016 #1 #1680b4 #68000000', 
        ' #2d0b #db680000 #a00002 #b6d0000 #b800180 #6e000014 #b460001', 
        ' #5000000 #28c000 #2800060 #2d #3 #5168a00 #51400', 
        ' #168a0000 #45140005 #1542a81 #5000 #68000140 #1 #16da0', 
        ' #40005c00 #5a1 #c0000b40 #4a140005 #1680001 #50001400 #68000168', 
        ' #1405001 #16800000 #2800016d #18000a #16d00 #14000a30 #2d0000', 
        ' #a3000 #b8 #28a00028 #6000 #c0002d14 #28000628 #a28000', 
        ' #a0000140 #e80000 #1d000000 #b80000 #40000503 #2d00001 #40000300', 
        ' #1700001 #44050b80 #24489122 #22448912 #49249249 #92492492 #24924924', 
        ' #12244891 #48224489 #89122092 #48912244 #44891224 #24489122 #22448912', 
        ' #12244811 #91024489 #89122448 #48824920 #44891224 #24489122 #22448912', 
        ' #12204891 #91224481 #12209018 #22448809 #12244891 #91224489 #41122448', 
        ' #22448910 #12244891 #44880901 #24489122 #22448912 #12244891 #44104489', 
        ' #24489122 #22448912 #12244891 #91224489 #89122448 #0:2 #3000000 ]', ), )
    p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=pickedEndEdges, minSize=bias_min, 
        maxSize=bias_max, constraint=FINER)
    p.seedPart(size=globe, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()


def Mesh_Stent_Finer(model, bias_min, bias_max, globe):
    p = mdb.models[model].parts['Stent']
    e = p.edges
    ## Edge selection extracted from Abaqus RPY file after meshing manually
    pickedEndEdges = e.getSequenceFromMask(mask=(
        '[#14000 #b4740000 #300000a3 #a0000000 #16d03 #d16d0000 #16800016', 
        ' #50000 #a0000000 #80001400 #68000016 #1 #1680b4 #68000000', 
        ' #2d0b #db680000 #a00002 #b6d0000 #b800180 #6e000014 #b460001', 
        ' #5000000 #28c000 #2800060 #0:2 #5028000 #51400 #168a0000', 
        ' #45140005 #1542a81 #5000 #68000140 #1 #16da0 #40000c00', 
        ' #a1 #c0000b40 #4a140005 #1680001 #50001400 #68000168 #1405001', 
        ' #16800000 #2800016d #18000a #16d00 #14000a30 #2d0000 #0', 
        ' #b8 #28a00028 #6000 #c0002d14 #28000628 #a28000 #a0000140', 
        ' #e80000 #1d000000 #180000 #0 #2d00000 #40000000 #1700001', 
        ' #44050b80 #24089122 #22448812 #49249249 #92492492 #24924924 #12244891', 
        ' #224489 #89100000 #48912244 #44891224 #24089122 #22448812 #12244801', 
        ' #91024089 #89122448 #40000000 #44891224 #24489122 #22448912 #12004891', 
        ' #91224480 #12209008 #22448009 #12244891 #91224489 #22448 #22448910', 
        ' #12244891 #44800900 #24489122 #22448110 #12244891 #44000089 #24489122', 
        ' #22448912 #12244891 #91224489 #89002448 #0:2 #3000000 ]', ), )
    p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=pickedEndEdges, minSize=bias_min, 
        maxSize=bias_max, constraint=FINER)
    p.seedPart(size=globe, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

#------------------------------------------------------------------------------
## Contact 
def Create_Interaction_Property(model, name):
    mdb.models[model].ContactProperty(name)
    mdb.models[model].interactionProperties[name].TangentialBehavior(formulation=FRICTIONLESS)

def Create_Interaction_Excl(model, step, int_name, excluded_instances, int_prop, deactivate_step):
    mdb.models[model].ContactExp(name=int_name, createStepName=step)
    mdb.models[model].interactions[int_name].includedPairs.setValuesInStep(stepName=step, useAllstar=ON)
    for instance in excluded_instances:
        try: 
            r12=mdb.models[model].rootAssembly.instances[instance].surfaces['External']
            r22=mdb.models[model].rootAssembly.instances[instance].surfaces['Internal']
            mdb.models[model].interactions[int_name].excludedPairs.setValuesInStep(stepName=step, addPairs=((ALLSTAR, r12), (ALLSTAR, r22)))
        except:
            r12=mdb.models[model].rootAssembly.instances[instance].surfaces['External']
            mdb.models[model].interactions[int_name].excludedPairs.setValuesInStep(stepName=step, addPairs=((ALLSTAR, r12), ))
    mdb.models[model].interactions[int_name].contactPropertyAssignments.appendInStep(stepName=step, assignments=((GLOBAL, SELF, int_prop), ))
    if deactivate_step != 'None':
        mdb.models[model].interactions[int_name].deactivate(deactivate_step)

def Create_Interaction_Incl(model, step, int_name, included_instances, included_surfaces, int_prop, deactivate_step):
    mdb.models[model].ContactExp(name=int_name, createStepName=step)
    mdb.models[model].interactions[int_name].includedPairs.setValuesInStep(stepName=step, useAllstar=ON)
    for instance, surface in zip(included_instances, included_surfaces):
            r12=mdb.models[model].rootAssembly.instances[instance[0]].surfaces[surface[0]]
            r22=mdb.models[model].rootAssembly.instances[instance[1]].surfaces[surface[1]]
            mdb.models[model].interactions[int_name].includedPairs.setValuesInStep(stepName=step, useAllstar=OFF, addPairs=((r12, r22), ))
    mdb.models[model].interactions[int_name].contactPropertyAssignments.appendInStep(stepName=step, assignments=((GLOBAL, SELF, int_prop), ))
    if deactivate_step != 'None':
        mdb.models[model].interactions[int_name].deactivate(deactivate_step)

#------------------------------------------------------------------------------
## Other

def Set_Output_Interval(model, numIntervals):
    mdb.models[model].fieldOutputRequests['F-Output-1'].setValues(numIntervals=numIntervals)

def Mass_Scale(model, start_step, factor):
    mdb.models[model].steps[start_step].setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, factor, 0.0, None, 0, 0, 0.0, 0.0, 0, None), ), improvedDtMethod=ON)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#action
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
mdb.Model(name=myString)

#Import Expander, Anatomy and bender from stl Mesh
Import_Stl_file(myString,myE,ExpanderPath)

Import_Stl_file(myString,myW,BenderPath)
Import_Stl_file(myString,myV, AnatomyPath)

#Import Stent from exported part file 
Import_ACIS_file(myString,myS,StentPath)

## Mesh Stent 
Mesh_Stent_Finer(myString, 0.1, 0.6, 0.08)

Create_Set_All_Elements(myString, myE, 'ALL-Elements')
Create_Set_All_Elements(myString, myW, 'ALL-Elements')
Create_Set_All_Elements(myString, myV, 'ALL-Elements')
Create_Set_All_Cells(myString, myS, 'ALL-Cells')

#Access ref points for the expander
ExpanderPoints = Read_XYZ_from_file(name+str('/ExpanderPoints.txt'))

#Access ref points for the bender
BenderPoints = Read_XYZ_from_file(name+str('/BenderPoints.txt')) 

##Add instances to assembly
Instance_to_Assembly(myString, myE, myE+str('-1'))
Instance_to_Assembly(myString, myE, myE+str('-2'))
Instance_to_Assembly(myString, myW, myW+str('-1'))
Instance_to_Assembly(myString, myS, myS+str('-1'))
Rotate_Instance(myString, myS+str('-1'), (1,0,0), (-2,0,0), 180)
Instance_to_Assembly(myString, myV, myV+str('-1'))

## Create Ref Geo for instance alignment 
move_axis = Create_Datum_Axis_Two_Points(myString, myS, StentXPath)
fix_axis = Datum_Axis_From_Principle(myString, myW, XAXIS)
Edge_Edge_Constraint(myString, myS+str('-1'), myW+str('-1'), move_axis, fix_axis, OFF)

##Create Stent 
Create_Material(myString, 'Steel', 7.5e-09, 193000.0, 0.3, ((200.0, 0.0), (410.0, 0.02), (530.0, 0.1), (600.0, 0.2), (630, 0.35), (640, 0.4)), 1)

Create_Solid_Homogeneous_Section(myString, 'Steel-Section', 'Steel' )
Section_Assignment(myString, myS, 'ALL-Cells', 'Steel-Section')

Create_Surface_Section(myString, 'Surface-Section', 7.5E-09)
Section_Assignment(myString, myE, 'ALL-Elements', 'Surface-Section')
Section_Assignment(myString, myW, 'ALL-Elements', 'Surface-Section')


##Create Conduit/Tissue sets
diameter = Read_XYZ_from_file(name+str('/VesselDiameter.txt'))
expander_dispSVC = np.ceil(1.2*diameter[0][0])/2-r_expander
expander_dispIVC = np.ceil(1.1*diameter[0][1])/2-r_expander
wall_thick = diameter[0][0]*0.05

ConduitPoints = Read_XYZ_from_file(name+str('/ConduitPoints.txt'))

CutPlane = Read_XYZ_from_file(name+str('/CutPlane.txt'))

Create_Set_Conduit(myString,myV,ConduitPoints, diameter[0][1], CutPlane[1], CutPlane[0], 'Conduit')

Create_Set_Boolean(myString, myV, 'ALL-Elements', 'Conduit', 'Tissue')

## Assign Conduit/Tissue Materials 
Create_Material_Yeoh_TD(myString, 'Tissue', 1.06e-09, ((0.0, 0.0), (0.0547, 0.1), (0.228, 0.2), (0.319, 0.3), (1.12, 0.4), (3.05, 0.5), (4.11, 0.6), (4.95, 0.7)), 1)
Create_Material_Yeoh_TD(myString, 'Conduit', 5.54e-010, ((0.0,0.0), (0.453,0.101),(0.680,0.164),(1.36,0.228),(2.57,0.29),(3.85,0.345),(6.88,0.441),(7.78,0.472)), 1)

Create_Shell_Section(myString,'Tissue-Section', wall_thick, 'Tissue')
Section_Assignment(myString, myV, 'Tissue', 'Tissue-Section')

Create_Shell_Section(myString,'Conduit-Section', wall_thick, 'Conduit')
Section_Assignment(myString, myV, 'Conduit', 'Conduit-Section')

##Create Steps for bending and expansion 
Create_Dynamic_Displacement_Step(myString, 'Bend', 'Initial', 0.25, ON)
Create_Dynamic_Displacement_Step(myString, 'Expand1', 'Bend', 0.25, ON)
Create_Dynamic_Displacement_Step(myString, 'Expand2', 'Expand1', 2, ON)
Create_Dynamic_Displacement_Step(myString, 'Recoil', 'Expand2', 2, ON)
Create_Dynamic_Displacement_Step(myString, 'Rest', 'Recoil', 2, ON)

## Create Amplitude Data 
bend_amp_data = ((0.0, 0.0), (0.24, 1.0))
Create_Smooth_Step_Amplitude(myString, 'Smooth-Bend', STEP, bend_amp_data)

expand1_amp_data = ((0,0,), (0.23,1))
Create_Smooth_Step_Amplitude(myString, 'Smooth-Expand1', STEP, expand1_amp_data)

i = Create_Datum_Cylind_CSYS(myString, (0,0,0), (1,0,0), (0,1,0), 'Cylind')

check_2 = (0.45*diameter[0][0]-r_expander)/expander_dispSVC
expand_amp_data = ((0.,0.), (1.2,check_2), (2,1.))

Create_Smooth_Step_Amplitude(myString, 'Smooth-Expand2', STEP, expand_amp_data)

recoil_amp_data = ((0,1.), (0.8,check_2), (2,0))
Create_Smooth_Step_Amplitude(myString, 'Balloon-Recoil', STEP, recoil_amp_data)

#Create Bending BCs and constraints
BenderDisplacements = Read_XYZ_from_file(name+str('/BenderDisplacement.txt'))
BenderRotations = Read_XYZ_from_file(name+str('/BenderRotations.txt'))
BenderRotations.append([0,0,0])

for d, rot, n in zip(BenderDisplacements, BenderRotations, range(0,len(BenderPoints))): 
    Create_Ref_Point(myString, BenderPoints[n])
    Create_Set_Node_Ring(myString, myW, BenderPoints[n], r_bender+0.5, 'NodeRing-'+str(n+1))
    Rigid_Ring_Constraint(myString, BenderPoints[n], myW+str('-1'), 'NodeRing-'+str(n+1),'RingConstraint'+str(n+1))
    Apply_Displacement_BC_Bend(myString, BenderPoints[n], 'BendBC-'+str(n+1), 'Bend', 'Expand1', 'Smooth-Bend', d[0], d[1], d[2], rot[0], rot[1], rot[2])

## Create and Merge Sets

SVC_indices = []
IVC_indices = []
for n in range(1,len(ExpanderPoints)-1):
    if ExpanderPoints[n][2] > 0:
        SVC_indices.append(n)
    else: 
        IVC_indices.append(n)

Create_Set_Node_Rings(myString, myE, [0, len(ExpanderPoints)-1], r_expander+0.2, 'Expander-Ends')
Create_Set_Node_Rings(myString, myE, SVC_indices, r_expander+0.2, 'SVC-Expander')
Create_Set_Node_Rings(myString, myE, IVC_indices, r_expander+0.2, 'IVC-Expander')

## Create Expand1 BCs
Apply_Displacement_BC_Expand(myString, myE+str('-1'), 'SVC-Expander', 'WireExpandBC-SVC','Smooth-Expand1', i, 'Expand1', 5-r_expander, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-1'), 'IVC-Expander', 'WireExpandBC-IVC','Smooth-Expand1', i, 'Expand1', 5-r_expander, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-1'), 'Expander-Ends', 'WireExpandBC-Ends','Smooth-Expand1', i, 'Expand1', 0.0, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')

## Create Expand2 BCs
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'SVC-Expander', 'ExpandBC-SVC','Smooth-Expand2', i, 'Expand2', expander_dispSVC, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'IVC-Expander', 'ExpandBC-IVC','Smooth-Expand2', i, 'Expand2', expander_dispIVC, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'Expander-Ends', 'ExpandBC-Ends','Smooth-Expand2', i, 'Expand2', 0.0, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')

## Create Recoil BCs
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'SVC-Expander', 'RecoilBC-SVC','Balloon-Recoil', i, 'Recoil', expander_dispSVC, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'IVC-Expander', 'RecoilBC-IVC','Balloon-Recoil', i, 'Recoil', expander_dispIVC, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')
Apply_Displacement_BC_Expand(myString, myE+str('-2'), 'Expander-Ends', 'RecoilBC-Ends','Balloon-Recoil', i, 'Recoil', 0.0, 0.0, 0.0, UNSET, UNSET, UNSET, 'None')

## Pin Stent in Recoil 
Create_Set_Nodes(myString, myS, [[-3.64658,-1.510496,36.5],[3.646625,-1.510388,36.499985],[3.64662,1.5104,36.499985]], 'Recoil-Pins')
Apply_Displacement_BC_Pin(myString, myS+str('-1'), 'Recoil-Pins', i, 'Recoil-Pins', 'Recoil', UNSET, UNSET, UNSET, 0, UNSET, UNSET, UNSET)

## Symmetry at PA Openings
Create_Symmetry_BC('/VesselTopBranch.txt', 'SVCBranch')
Create_Symmetry_BC('/VesselLeft.txt', 'LPA')
Create_Symmetry_BC('/VesselRight.txt', 'RPA')
Create_Symmetry_BC('/VesselRightBranch.txt', 'RPABranch')
Create_Symmetry_BC('/VesselTop.txt', 'SVC')
Create_Symmetry_BC('/VesselBase.txt', 'IVC')


##Set Element Types
Set_SFM3D3_Elem_Type(myString, myE, 'ALL-Elements')
Set_SFM3D3_Elem_Type(myString, myE, 'ALL-Elements')
Set_SFM3D3_Elem_Type(myString, myW, 'ALL-Elements')
Set_C3D8R_Elem_Type(myString, myS, 'ALL-Cells')
Set_S3R_Elem_Type(myString, myV, 'ALL-Elements')

## Create Contact Interactions 
Create_Interaction_Property(myString, 'Int-Prop-1')

Create_Surface_Orphan_Mesh(myString, myE, 'internal')
Create_Surface_Orphan_Mesh(myString, myE, 'external')
Create_Surface_Orphan_Mesh(myString, myW, 'internal')
Create_Surface_Orphan_Mesh(myString, myW, 'external')
Create_Surface_Orphan_Mesh(myString, myV, 'internal')
Create_Surface_Orphan_Mesh(myString, myV, 'external')
Create_Surface_Orphan_Mesh(myString, myS, 'external')

Create_Interaction_Incl(myString, 'Bend', 'Contact-Bending', [[myV+str('-1'), myW+str('-1')]], [['Internal', 'External']], 'Int-Prop-1', 'Expand1')
Create_Interaction_Incl(myString, 'Expand1', 'Contact-Expansion1', [[myV+str('-1'), myE+str('-1')]], [['Internal', 'External']], 'Int-Prop-1', 'Expand2')
Create_Interaction_Excl(myString, 'Expand2', 'Contact-Expansion2', [myW+str('-1'), myE+str('-1')], 'Int-Prop-1', 'None')

## Mass Scaling 
Mass_Scale(myString, 'Bend', 10000)

##Output Interval 
Set_Output_Interval(myString, 50)

