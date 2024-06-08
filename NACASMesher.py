 #   This program is free software: you can redistribute it and/or modify
 #   it under the terms of the GNU General Public License as published by
 #   the Free Software Foundation, either version 3 of the License, or
 #   (at your option) any later version.

 #  This program is distributed in the hope that it will be useful,
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #   GNU General Public License for more details.

 #   You should have received a copy of the GNU General Public License
 #   along with this program.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np
import math
import matplotlib.pyplot as plt
import csv

#=======================================================
# NACA CURVE GENERATION FUNCTIONS
#=======================================================
def dyc_dx( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    ((2.0 * m) / np.power(p,2)) * (p - x / c),
                    ((2.0 * m) /np.power(1-p,2)) * (p - x / c ))

def thickness( x, t, c ):
    A = 0.2969
    B = -0.1260
    C = -0.3516
    D = 0.2843
    E = -0.1036 #Closed circular trailing edge
    #E = -0.1015

    AxMinus1 =  A * (np.sqrt(x/c))
    Bx1 =  B * (x/c)
    Cx2 =  C * np.power(x/c,2)
    Dx3 =  D * np.power(x/c,3)
    Ex4 =  E * np.power(x/c,4)

    return 5*t*c*(AxMinus1 + Bx1 + Cx2 + Dx3 + Ex4)

def camber_line( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    m * (x / np.power(p,2)) * (2.0 * p - (x / c)),
                    m * ((c - x) / np.power(1-p,2)) * (1.0 + (x / c) - 2.0 * p ))

def naca4(x, m, p, t, c=1):
    if(m == 0 and p == 0):
        yt = thickness(x, t, c)
        return ((x, yt), 
                (x, -yt))
    else:   
        dyc_dx_var = dyc_dx(x, m, p, c)
        th = np.arctan(dyc_dx_var)

        yt = thickness(x, t, c)
        yc = camber_line(x, m, p, c)  

        return ((x - yt*np.sin(th), yc + yt*np.cos(th)), 
                (x + yt*np.sin(th), yc - yt*np.cos(th)))

#=======================================================
# CONFIGURATION
#=======================================================

#NACA CONFIGURATION

#naca0012 
m = 0.0
p = 0.0
t = 0.12
c = 1.0



#MESH CONFIGURATION
vector_option = True

parameter_vector = [4000, 12, 20, 0, 1, 1, 0.001, 0.003, 0.035, 0.4, 0.3, 0.000005, 1.31, 1, 1, 1, 17, 31, 74, 21, 20, 0.25]

if(vector_option is True):
    number_of_points = parameter_vector[0] 
    DistanceTo_Inlet = parameter_vector[1]
    DistanceTo_Outlet = parameter_vector[2]
    Angle_Of_Response = parameter_vector[3]
    Z_Thickness = parameter_vector[4]
    Mesh_Scale = parameter_vector[5]
    LeadingEdge_CellSize = parameter_vector[6]
    TrailingEdge_CellSize = parameter_vector[7]
    Middle_CellSize = parameter_vector[8]
    SeparatingPoint_Position = parameter_vector[9]
    BL_Thickness = parameter_vector[10]
    FirstLayer_Thickness = parameter_vector[11]
    Expansion_Ratio = parameter_vector[12]
    Inlet_MaxCellSize = parameter_vector[13]
    Outlet_MaxCellSize = parameter_vector[14]
    Inlet_X_Outlet_MaxCellSize = parameter_vector[15]
    MeshNum_OnBL = parameter_vector[16]
    MeshNum_OutBL = parameter_vector[17]
    MeshNum_AtTail = parameter_vector[18]
    MeshNum_Leading = parameter_vector[19]
    MeshNum_Trailing = parameter_vector[20]
    Inlet_ExpansionRatio_2 = parameter_vector[21]
else:
    number_of_points = 4000

    DistanceTo_Inlet = 12 # Distance to inlet (x chord length) - DEFAULT: 12
    DistanceTo_Outlet = 20 # Distance to outlet (x chord length) - DEFAULT: 20
    Angle_Of_Response = 0 # Angle of response (degree) - DEFAULT: 0
    Z_Thickness = 1 # Depth in Z direction - DEFAULT: 1
    Mesh_Scale = 1 # Mesh Scale - DEFAULT: 1

    LeadingEdge_CellSize = 0.01 # Cell size at leading edge - DEFAULT: 0.01
    TrailingEdge_CellSize = 0.03 # Cell size at trailing edge - DEFAULT: 0.03
    Middle_CellSize = 0.035 # Cell size in middle - DEFAULT: 0.035
    SeparatingPoint_Position = 0.4 # Separating point position (% from leading point) - DEFAULT: 0.4

    BL_Thickness = 0.05 # Boundary layer thickness - DEFAULT: 0.5
    FirstLayer_Thickness = 5e-7 # Fist layer thickness - DEFAULT: 0.005
    Expansion_Ratio = 2.3 # Expansion ratio - DEFAULT: 1.2
    #0.2; 5e-7; 2 | 1.279989e-02
    #0.1; 5e-7; 2 | 1.19049999903e-02
    #0.05; 5e-7; 2 | 1.171929e-02

    Inlet_MaxCellSize = 2 # Inlet Max Cell Size - DEFAULT: 1
    Outlet_MaxCellSize = 200 # Outlet Max Cell Size - DEFAULT: 1
    Inlet_X_Outlet_MaxCellSize = 50 # Inlet_X_Outlet_MaxCellSize - DEFAULT: 1

    #Secondary
    MeshNum_OnBL = 17 # Parameter - DEFAULT: 17
    MeshNum_OutBL = 31 # Parameter - DEFAULT: 31
    MeshNum_AtTail = 220 # Parameter - DEFAULT: 74
    MeshNum_Leading = 400 # Parameter - DEFAULT: 21
    MeshNum_Trailing = 200 # Parameter - DEFAULT: 20

    Inlet_ExpansionRatio_2 = 0.25 # Parameter - DEFAULT: 0.25

#=======================================================
# NACA POINTS CREATION FINAL FUNCTIONS
#=======================================================

number_of_points = round(number_of_points/2)
#number_of_points -= 1

x = np.linspace(1,0,number_of_points)

points = np.round(naca4(x, m, p, t, c),6)

X = np.concatenate((points[0][0],points[1][0][::-1]))
Y = np.concatenate((points[0][1],points[1][1][::-1]))

df = np.column_stack((X,Y))

pair = [0.000000, -0.000000];
#print(Coords)
index_of_midpoint = np.where((df == pair).all(axis=1))[0][0]
#print(index_of_midpoint)

#=======================================================
# MESH GENERATION FUNCTIONS AND FILE WRITTING
#=======================================================

Division_Point = SeparatingPoint_Position 
ExpansionRatio_Leading = Middle_CellSize/LeadingEdge_CellSize
ExpansionRatio_Trailing = Middle_CellSize/TrailingEdge_CellSize
Inlet_ExpansionRatio_1 = TrailingEdge_CellSize/Inlet_MaxCellSize

LastLayer_Thickness = FirstLayer_Thickness*Expansion_Ratio**MeshNum_OnBL
ExpansionRatio_5 = Expansion_Ratio**MeshNum_OnBL
ExpansionRatio_6 = Inlet_MaxCellSize/LastLayer_Thickness
ExpansionRatio_7 = Outlet_MaxCellSize/TrailingEdge_CellSize
ExpansionRatio_AtOutlet = Inlet_X_Outlet_MaxCellSize*ExpansionRatio_6/Outlet_MaxCellSize*(MeshNum_OutBL+MeshNum_OnBL)/MeshNum_OutBL


ExpansionRatio1_a = LastLayer_Thickness/DistanceTo_Inlet*(ExpansionRatio_6-1)+1
ExpansionRatio1_b = math.log(ExpansionRatio_6)/math.log(ExpansionRatio1_a)

ExpansionRatio2_a = TrailingEdge_CellSize/DistanceTo_Outlet*(ExpansionRatio_7-1)+1
ExpansionRatio2_b = math.log(ExpansionRatio_7)/math.log(ExpansionRatio2_a)

ExpensionRatio3_a = LeadingEdge_CellSize/SeparatingPoint_Position*(ExpansionRatio_Leading-1)+1
ExpensionRatio3_b = math.log(ExpansionRatio_Leading)/math.log(ExpensionRatio3_a)

ExpensionRatio4_a = TrailingEdge_CellSize/(1-Middle_CellSize)*(ExpansionRatio_Trailing-1)+1
ExpensionRatio4_b = math.log(ExpansionRatio_Trailing)/math.log(ExpensionRatio4_a)


# Writes blockMeshDict file into system directory of current folder
f = open('system/blockMeshDict' ,'w')
f.write('/*--------------------------------*- C++ -*----------------------------------*\ \n')
f.write('  =========                 | \n')
f.write('  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox \n')
f.write('   \\    /   O peration     | Website:  https://openfoam.org \n')
f.write('    \\  /    A nd           | Version:  6 \n')
f.write('     \\/     M anipulation  | \n')
f.write('\*---------------------------------------------------------------------------*/ \n')
f.write('FoamFile \n')
f.write('{ \n')
f.write('    version     2.0; \n')
f.write('    format      ascii; \n')
f.write('    class       dictionary; \n')
f.write('    object      blockMeshDict; \n')
f.write('} \n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
f.write('convertToMeters '+str(Mesh_Scale) +'; \n')
f.write('\n')

f.write('geometry\n')
f.write('{\n')
f.write('}\n')
f.write('\n')

f.write('vertices\n')
f.write('(\n')
f.write(' 	('+ str(0)						+ 	' '	+ 	str(0)													+	' '	+ 	str(0) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+ 	str(0)													+	' '	+ 	str(0) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+ 	str(DistanceTo_Inlet)									+	' '	+	str(0) +') \n')
f.write(' 	('+ str(-DistanceTo_Inlet+1)	+   ' '	+ 	str(0)													+	' '	+	str(0) +') \n')
f.write(' 	('+ str(0)						+ 	' '	+ 	str(0)													+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+ 	str(0)													+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+ 	str(DistanceTo_Inlet)									+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(-DistanceTo_Inlet+1)	+ 	' '	+	str(0)													+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(math.sin(math.pi/180*Angle_Of_Response)*(DistanceTo_Outlet+1))	+	' '	+	str(0) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(DistanceTo_Inlet)									+	' '	+	str(0) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(math.sin(math.pi/180*Angle_Of_Response)*(DistanceTo_Outlet+1))	+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(DistanceTo_Inlet)									+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+	str(-DistanceTo_Inlet)									+	' '	+	str(0) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+	str(-DistanceTo_Inlet)									+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(-DistanceTo_Inlet)									+	' '	+	str(0) +') \n')
f.write(' 	('+ str(DistanceTo_Outlet+1)	+ 	' '	+	str(-DistanceTo_Inlet)									+	' '	+	str(Z_Thickness) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+	str(0)													+	' '	+	str(0) +') \n')
f.write(' 	('+ str(1)						+ 	' '	+	str(0)													+	' '	+	str(Z_Thickness) +') \n')

f.write('  \n')
f.write('); \n')
f.write('  \n')
f.write('\n')

f.write('blocks\n')
f.write('(\n')
f.write('	 hex (0	1	2	3	4	5	6	7)	('+str(MeshNum_Leading+MeshNum_Trailing)+' '+ str(MeshNum_OnBL+MeshNum_OutBL) + ' 1	)	//block \n')
f.write('    edgeGrading \n')
f.write('    ( \n')

f.write('    //	x-direction	expansion ratio \n')
f.write('    ( \n')
f.write('    	(' + str(SeparatingPoint_Position) + ' ' + str(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(ExpansionRatio_Leading) +' ) \n')
f.write('    	(' + str(1-Division_Point) + ' ' + str(1-MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(1/ExpansionRatio_Leading) + ') \n')
f.write('    ) \n')
f.write('    ' + str(Inlet_ExpansionRatio_1) + ' ' + str(Inlet_ExpansionRatio_1) + ' \n')

f.write('    ( \n')
f.write('    	(' + str(Division_Point) + ' ' + str(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(ExpansionRatio_Leading) + ')\n')
f.write('    	(' + str(1-Division_Point) + ' ' + str(1-MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(1/ExpansionRatio_Trailing) + ') \n')
f.write('    ) \n')

f.write('    //	y-direction	expansion ratio \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-(BL_Thickness/DistanceTo_Inlet)) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-(BL_Thickness/DistanceTo_Inlet)) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-(BL_Thickness/DistanceTo_Inlet)) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('     \n')
f.write('    //	z-direction	expansion ratio \n')
f.write('    1	1	1	1 \n')
f.write('    ) \n')	
f.write('    \n')
f.write('    hex	(1	8	9	2	5	10	11	6)	(	'+ str(MeshNum_AtTail) + ' ' + str(MeshNum_OnBL + MeshNum_OutBL) + ' ' + '1)	//block	2 \n') 	
f.write('    edgeGrading \n')
f.write('    ( \n')
f.write('    //	x-direction	expansion ratio \n')
f.write('    	' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) +' \n')
f.write('    //	y-direction	expansion ratio \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('    '+  str(ExpansionRatio_AtOutlet) + '  ' + str(ExpansionRatio_AtOutlet) + ' \n')
f.write('    ( \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(ExpansionRatio_5) + ') \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(ExpansionRatio_6) + ') \n')
f.write('    ) \n')
f.write('     \n')
f.write('    //	z-direction	expansion ratio \n')
f.write('    1	1	1	1 \n')
f.write('    ) \n')	
f.write('     \n')		
f.write('    hex	(3	12	16	0	7	13	17	4)	(	'+str(MeshNum_Leading+MeshNum_Trailing)+' '+ str(MeshNum_OnBL+MeshNum_OutBL) + ' 1	)	//block	2 \n') 	
f.write('    edgeGrading \n')
f.write('    ( \n')
f.write('    //	x-direction	expansion ratio \n')
f.write('    	' + str(Inlet_ExpansionRatio_1) + ' \n')
f.write('    ( \n')
f.write('    	(' + str(Division_Point) + ' ' + str(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(ExpansionRatio_Leading) + ') \n')
f.write('    	(' + str(1-Division_Point) + ' ' + str(1-(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing))) + ' ' + str(1/ExpansionRatio_Trailing) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(Division_Point) + ' ' + str(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing)) + ' ' + str(ExpansionRatio_Leading) + ') \n')
f.write('    	(' + str(1-Division_Point) + ' ' + str(1-(MeshNum_Leading/(MeshNum_Leading+MeshNum_Trailing))) + ' ' + str(1/ExpansionRatio_Trailing) + ') \n')
f.write('    ) \n')
f.write('    	' + str(Inlet_ExpansionRatio_1) + ' \n')
f.write('    //	y-direction	expansion ratio \n')
f.write('    ( \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) + ') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) + ') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) + ') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) + ') \n')
f.write('    ) \n')
f.write('    ( \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) + ') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) + ') \n')
f.write('    ) \n')
f.write('     \n')
f.write('    //	z-direction	expansion ratio \n')
f.write('    1	1	1	1 \n')	
f.write('    ) \n')		
f.write('     \n')
f.write('     \n')
f.write('     \n')
f.write('     \n')
f.write('     \n')
f.write('    hex	(12	14	8	16	13	15	10	17)	(	' + str(MeshNum_AtTail) + ' ' + str(MeshNum_OnBL+MeshNum_OutBL)  + ' 1)	//block	4 \n')	
f.write('    edgeGrading \n')
f.write('    ( \n')
f.write('    //	x-direction	expansion ratio \n')
f.write('    	' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) + ' ' + str(ExpansionRatio_7) +' \n')
f.write('    //	y-direction	expansion ratio \n')
f.write('    ( \n')
f.write('    	(' + str(1-BL_Thickness/DistanceTo_Inlet) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) + ') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) + ') \n')
f.write('    ) \n')
f.write('    '+  str(1/ExpansionRatio_AtOutlet) + '  ' + str(1/ExpansionRatio_AtOutlet) + ' \n')
f.write('    ( \n')
f.write('    	(' + str(1-(BL_Thickness/DistanceTo_Inlet)) + ' ' + str(1-(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL))) + ' ' + str(1/ExpansionRatio_6) +') \n')
f.write('    	(' + str(BL_Thickness/DistanceTo_Inlet) + ' ' + str(MeshNum_OnBL/(MeshNum_OutBL+MeshNum_OnBL)) + ' ' + str(1/ExpansionRatio_5) +  ') \n')
f.write('    ) \n')
f.write('     \n')
f.write('    //	z-direction	expansion ratio \n')
f.write('    1	1	1	1 \n')
f.write('    ) \n')		
f.write(');\n\n')
f.write('     \n')
f.write('edges\n')
f.write('(\n')
f.write('   arc	3 2	(' + str(-DistanceTo_Inlet*math.sin(math.pi/4)+1) + ' ' + str(DistanceTo_Inlet*math.sin(math.pi/4)) + '  0) \n')
f.write('   arc	7 6	(' + str(-DistanceTo_Inlet*math.sin(math.pi/4)+1) + ' ' + str(DistanceTo_Inlet*math.sin(math.pi/4)) +' '+  str(Z_Thickness) +') \n')
f.write('     \n')
f.write('	spline	1	0   \n')
f.write('	(\n')
for i in range(0, index_of_midpoint+1):
#for i in range(10):
	linestring='    ('+str(X[i])+' '+str(Y[i]) +' '+ str(0)+')\n'
	f.write(linestring)
#f.write('	(  0	0	0)	   \n')
f.write('	)\n')
f.write('	\n')
f.write('	\n')
f.write('	spline	5	4	   \n')
f.write('	(\n')
for i in range(0, index_of_midpoint+1):
#for i in range(10):
	linestring='    ('+str(X[i])+' '+str(Y[i]) +' '+ str(Z_Thickness)+')\n'
	f.write(linestring)
f.write('	)\n')
f.write('	\n')
f.write('	\n')
f.write('	arc	3 12	(' + str(-DistanceTo_Inlet*math.sin(math.pi/4)+1) + ' ' + str(-DistanceTo_Inlet*math.sin(math.pi/4)) + '  0) \n')
f.write('	arc	7 13	(' + str(-DistanceTo_Inlet*math.sin(math.pi/4)+1) + ' ' + str(-DistanceTo_Inlet*math.sin(math.pi/4)) +' '+  str(Z_Thickness) +') \n')
f.write('	\n')
f.write('	spline	0	16   \n')
f.write('	(\n')
#for i in range(10):
for i in range(index_of_midpoint+2, len(X)):
	linestring='    ('+str(X[i])+' '+str(Y[i])+' '+  str(0) + ')\n'
	f.write(linestring)
f.write('	)\n')
f.write('	\n')
f.write('	\n')
f.write('	spline	4	17	  \n')
f.write('	(\n')
for i in range(index_of_midpoint+2, len(X)):
#for i in range(10):
	linestring='    ('+str(X[i])+' '+str(Y[i])+' '+  str(Z_Thickness) + ')\n'
	f.write(linestring)

f.write('	)\n')
f.write('	\n')
f.write('	\n')
f.write(');\n')
f.write('	\n')
f.write('faces\n')
f.write('(\n')
f.write('	\n')
f.write(');\n\n')
f.write('faces\n')
f.write('(\n')
f.write('	\n')
f.write(');\n')
f.write('	\n')
f.write('	\n')
f.write('defaultPatch')
f.write('{\n')
f.write('	name frontAndBack;\n')
f.write('	type  empty;\n')
f.write('	\n')
f.write('}\n')
f.write('	\n')
f.write('boundary\n')
f.write('(\n')
f.write('    inlet\n')
f.write('    {\n')
f.write('        type patch;\n')
f.write('        faces\n')
f.write('        (\n')
f.write('          (9 2 6 11)  	\n')
f.write('          (2 3 7 6  )  \n')
f.write('          (3 12 13 7)  \n')
f.write('          (12 15 14 13) \n')
f.write('         );\n')
f.write('     }\n')
f.write('    outlet\n')
f.write('    {\n')
f.write('        type patch;\n')
f.write('        faces\n')
f.write('        (\n')
f.write('          (8 9 10 11)  	 \n')
f.write('          (15 8 10 14  ) 	 \n')
f.write('         );\n')
f.write('     }\n')
f.write('    walls\n')
f.write('    {\n')
f.write('        type wall;\n')
f.write('        faces\n')
f.write('        (\n')
f.write('          (0 1 5 4)  	\n')
f.write('          (0 4 17 16  ) \n')
f.write('         );\n')
f.write('     }\n')
f.write('    interface1\n')
f.write('    {\n')
f.write('        type patch;\n')
f.write('        faces\n')
f.write('        (\n')
f.write('          (1 8 10 5)  	\n')
f.write('         );\n')
f.write('     }\n')
f.write('    interface2\n')
f.write('    {\n')
f.write('        type patch;\n')
f.write('        faces\n')
f.write('        (\n')
f.write('          (16 17 10 8)   	\n')
f.write('         );\n')
f.write('     }\n')
f.write(');\n')
f.write('mergePatchPairs\n')
f.write('(\n')
f.write('( interface1 interface2 )\n')
f.write(');\n')

f.close()
