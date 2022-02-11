#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
\geometricManipulations\__init__.py: This module has been build to have an standard class to interact with pointcloud data, lines, planes and statisticall information
"""
__author__          = "Robert Rijnbeek"
__version__         = "1.0.1"
__maintainer__      = "Robert Rijnbeek"
__email__           = "robert270384@gmail.com"
__status__          = "Development"

__creation_date__   = '11/02/2022'
__last_update__     = '11/02/2022'

# =============== IMPORTS ===============

import numpy as np
import time
import random

from math import cos, sin, radians, sqrt
from scipy.spatial import distance
from scipy.stats import binned_statistic_dd

# ===============  CODE  ===============

class GeometricManipulation:
    def __init__(self):
        pass

# ===== P O I N T    D I S T A N C E =======

    def Absolute3DDistance(self,P1, P2):
        """
        INPUT:
            - P1 (TUPLE, LIST OR ARRAY) Lenght 3
            - P1 (TUPLE, LIST OR ARRAY) Lenght 3
        OUTPUT:
            - FLOAT: Distance between P1 an p2
        """
        return sqrt((P1[0]-P2[0])**2+(P1[1]-P2[1])**2+(P1[2]-P2[2])**2)

    def ClosestPoint(self,point, point_list):
        """
        INPUT:
            - point (TUPLE, LIST OR ARRAY) Lenght 3
            - point_list List of point (TUPLE, LIST OR ARRAY with  Lenght 3)
        OUTPUT:
            - point of pointlist that is nearest to point (first argument)
        """
        distances=distance.cdist([point], point_list)

        closest_index = distances.argmin()
        #closest_index = distance.cdist([point], point_list).argmin()
        return point_list[closest_index], np.amin(distances)

# ==== T R A N S F O R M A T I O N   M A T R I X =====

    def TransformationMatrixByVector(self,INIT_VECTOR, FINAL_VECTOR):
        """
        INPUT:
            - INIT_VECTOR: (TUPLE, LIST OR ARRAY) Lenght 3 That is the direction of the vector.
            - FINAL_VECTOR: (TUPLE, LIST OR ARRAY) Lenght 3 That is the direction THAT WILL BE CONVERTED THE INIT_VECTOR
        OUTPUT:
            MATRIX 4 x 4 : That represent the transformation matrix
        """
        try:
            f=INIT_VECTOR
            t=FINAL_VECTOR
            c = np.dot(f, t)
            if abs(1-c) > 0.000000001 :
                v = np.cross(f, t)
                u = v/np.linalg.norm(v)
                h = (1 - c)/(1 - c**2)
                vx, vy, vz = v
                rot =[[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy,0],
                    [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx,0],
                    [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2,0],
                    [0,0,0,1]]
                return rot
            else:
                print("for resolution reasons  the transformation matrix is the base transformation matrix"  )
                rot =[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
                return rot
        except:
            print("error creating transformation matriz")
            return False
    
    def Transformation3x3MatrixByVector(self,INIT_VECTOR, FINAL_VECTOR):
        """
        INPUT:
            - INIT_VECTOR: (TUPLE, LIST OR ARRAY) Lenght 3 That is the direction of the vector.
            - FINAL_VECTOR: (TUPLE, LIST OR ARRAY) Lenght 3 That is the direction THAT WILL BE CONVERTED THE INIT_VECTOR
        OUTPUT:
            MATRIX 4 x 4 : That represent the transformation matrix
        """
        try:
            f=INIT_VECTOR
            t=FINAL_VECTOR
            c = np.dot(f, t)
            if abs(1-c) > 0.000000001 :
                v = np.cross(f, t)
                u = v/np.linalg.norm(v)
                h = (1 - c)/(1 - c**2)
                vx, vy, vz = v
                rot =[[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
                    [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
                    [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]]
                return np.array(rot)
            else:
                print("for resolution reasons  the transformation matrix is the base transformation matrix"  )
                rot =[[1,0,0],[0,1,0],[0,0,1]]
                return np.array(rot)
        except:
            print("error creating transformation matriz")
            return False

    def trig(self,angle):
        """
        INPUT: angle
        OUTPUT: The cos and the sin of the angle
        """
        r = radians(angle)
        return cos(r), sin(r)

    def TransformationMatrix(self,rotation=(0,0,0), translation=(0,0,0)): # redefine function name
        """
        INPUT:
            - rotation (TUPLE, LIST OR ARRAY) Lenght 3 That represent the rotation angle components.
            - translation: (TUPLE, LIST OR ARRAY) Lenght 3 That represent the traslation
        OUTPUT:
            MATRIX 4 x 4 : That represent the transformation matrix that include the rotation and the traslation
        """
        xC, xS = self.trig(rotation[0])
        yC, yS = self.trig(rotation[1])
        zC, zS = self.trig(rotation[2])
        dX = translation[0]
        dY = translation[1]
        dZ = translation[2]
        return [[yC*xC, -zC*xS+zS*yS*xC, zS*xS+zC*yS*xC, dX],
                [yC*xS, zC*xC+zS*yS*xS, -zS*xC+zC*yS*xS, dY],
                [-yS, zS*yC, zC*yC, dZ],
                [0, 0, 0, 1]]


    def TransformationMatrix3x3(self,rotation=(0,0,0)): # redefine function name
        """
        INPUT:
            - rotation (TUPLE, LIST OR ARRAY) Lenght 3 That represent the rotation angle components.
            - translation: (TUPLE, LIST OR ARRAY) Lenght 3 That represent the traslation
        OUTPUT:
            MATRIX 4 x 4 : That represent the transformation matrix that include the rotation and the traslation
        """
        xC, xS = self.trig(rotation[0])
        yC, yS = self.trig(rotation[1])
        zC, zS = self.trig(rotation[2])
        return np.array([[yC*xC, -zC*xS+zS*yS*xC, zS*xS+zC*yS*xC],
                [yC*xS, zC*xC+zS*yS*xS, -zS*xC+zC*yS*xS],
                [-yS, zS*yC, zC*yC]])

    def CombinationOfTransformationMatrix(self,MATRIX1,MATRIX2): #see if we must use also for matrix 3*3
        """
        INPUT:
            - MATRIX1: MATRIX 4 x 4 : That represent the transformation matrix 1
            - MATRIX1: MATRIX 4 x 4 : That represent the transformation matrix 1
        OUTPUT: 
            - MATRIX 4 x 4 That repressent the transformation matrix that is combination of both MATRIX1 and MATRIX 2
        """
        try:
            return [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*MATRIX2)] for X_row in MATRIX1]
        except:
            print("problems creating new transformation Matrix with: MATRIX1="+str(MATRIX1) + " and MATRIX2=" + str(MATRIX2) )
            return False

    def CombinationOfTransformationMatrix3x3(self,MATRIX1,MATRIX2): #see if we must use also for matrix 3*3
        """
        INPUT:
            - MATRIX1: MATRIX 3 x 3 : That represent the transformation matrix 1
            - MATRIX1: MATRIX 3 x 3 : That represent the transformation matrix 1
        OUTPUT: 
            - MATRIX 3 x 3 That repressent the transformation matrix that is combination of both MATRIX1 and MATRIX 2
        """
        try:
            return np.dot(MATRIX1,MATRIX2)
        except:
            print("problems creating new transformation Matrix with: MATRIX1="+str(MATRIX1) + " and MATRIX2=" + str(MATRIX2) )
            return False


# ===== P O I N T    T R A N S F O R M A T I O N =======
    
    def PointTransform(self,point=(0,0,0), matrix=(0,0,0)):
        """
        INPUT: 
            - point: (TUPLE, LIST OR ARRAY) Lenght 3 That represent a single point
            - vector: Matrix 4 x 4 that represent the transformation matrix
        OUTPUT:
            - point: (TUPLE, LIST OR ARRAY) of Lenght 3 that represent the transformation of 'point'
        """
        p = [0,0,0]
        for r in range(3):
            p[r] += matrix[r][3]
            for c in range(3):
                p[r] += point[c] * matrix[r][c]
        return p

    def PointListTransform(self,pointlist,matrix):
        """
        INPUT: 
            - pointlist: List of (TUPLE, LIST OR ARRAY) Lenght 3 That represent a list of points
            - matrix: Matrix 4 x 4 that represent the transformation matrix
        OUTPUT:
            - pointlist: (TUPLE, LIST OR ARRAY) of Lenght 3 that represent the transformation of the 'pointlist'
        """
        return list(map(lambda x: self.PointTransform(x, matrix),pointlist))
    
    def PointListTransform3x3(self,pointlist,matrix):
        """
        INPUT: 
            - pointlist: ARRAY Lenght 3 That represent a list of points
            - matrix: array Matrix 3 x 3 that represent the transformation matrix
        OUTPUT:
            - pointlist:  ARRAY of Lenght 3 that represent the transformation of the 'pointlist'
        """
        return np.dot(pointlist, matrix.T)


# ==== P O I N T    S U B S A M P L E =======

    def RandomSubsampleOfPointList(self,POINTLIST,SAMPLE_COUNT):
        """
        INPUT:
            - POINTLIST: pointlist: List of (TUPLE, LIST OR ARRAY) Lenght 3 That represent a list of points
            - SAMPLE_COUNT (Integer): Count of points that will be the subset of points that wil, be selected
        OUTPUT: 
            - Pointlist:(TUPLE, LIST OR ARRAY) of Lenght 3 that represent the subset of the pointlist
        """
        try:
            if len(POINTLIST)>=SAMPLE_COUNT:
                return random.sample(POINTLIST, SAMPLE_COUNT)
            else:
                print("SAMPLE_COUNT is larger than lenght of pointlist")
                return False
        except:
            print("Error with 'RandomSubsampleOfPointList' function")
            return False

    """
    def DensitySubsambpleOfPointList(self,POINTLIST,DENSITY_OF_POINT_LIST)
        ...

    """
        

# ===== P L A N E     F I T ======

    def MakePlaneBy3Points(self, POINT_LIST):
        """
        INPUT:
            - pointlist: List of (TUPLE, LIST OR ARRAY) Lenght 3 That represent a list of 3 points
        OUTPUT:
            List of lenght 4 that represent the 4 coeficients [a, b, c, d] of the plane: ax+by+cz+d=0
        """
        [p1,p2,p3]=POINT_LIST
        v1=p3-p1
        v2=p2-p1
        cp=np.cross(v1,v2)
        [a,b,c] = cp
        d=np.dot(cp, p3)*(-1)
        return [a,b,c,d]/np.linalg.norm(cp)

    def MakePlaneBy3Points1(self, POINT_LIST):
        """
        INPUT:
            - pointlist: List of (TUPLE, LIST OR ARRAY) Lenght 3 That represent a list of 3 points
        OUTPUT:
            List of lenght 4 that represent the 4 coeficients [a, b, c, d] of the plane: ax+by+cz+d=0
        """
        print(POINT_LIST[:,0])
        return 
    
    def MakePlaneFromPointCloud(self,POINT_LIST):
        """
        INPUT:
            - pointlist: List of (TUPLE, LIST OR ARRAY) Lenght 3 That represent a list of points
        OUTPUT:
            List of lenght 4 that represent the 4 coeficients [a, b, c, d] of the plane: ax+by+cz+d=0
        """
        try:
            if len(POINT_LIST)>1:
                xyz=np.array(POINT_LIST)
                centroid         = xyz.mean(axis = 0)
                xyzT             = np.transpose(xyz)
                xyzR             = xyz - centroid                         #points relative to centroid
                #xyzRT            = np.transpose(xyzR)
                u, sigma, v       = np.linalg.svd(xyzR)
                normal            = v[2]                              
                normal            = normal / np.linalg.norm(normal)  
                d = -(normal[0] * centroid[0] + normal[1] * centroid[1] + normal[2] * centroid[2])
                return [normal[0],normal[1],normal[2],d]
            else:
                print("PointList must be have at least more than 2 points")
                return False
        except:
            print("no plane fit possible")
            return False

    def MakePlaneByNormalAndPoint(self,NORMAL,POINT): #Need to see if 'd' is correct maybe it is -1 for the other planes
        normal=np.array(NORMAL)
        normal=normal / np.linalg.norm(normal) 
        d=-(normal[0]*POINT[0]+normal[1]*POINT[1]+normal[2]*POINT[2])
        return [normal[0],normal[1],normal[2],d]

    def PlaneTransform(self,PLANE_COEF,MATRIX):
        [a,b,c,d]=PLANE_COEF
        normal=[a,b,c]
        point_of_plane=self.GetRandomPointFromPlane(PLANE_COEF)
        new_normal=self.PointTransform(point=normal, matrix=MATRIX)
        new_point=self.PointTransform(point=point_of_plane, matrix=MATRIX)
        plane_coef=self.MakePlaneByNormalAndPoint(new_normal,new_point)
        return plane_coef

# ===== P L A N E / P O I N T S     M A N I P U L A T I O N S =======

    def GetRandomPointFromPlane(self,PLANE_COEF): #NEW 13-6-2019
        [x1,y1,z1] = np.random.rand(3)
        [a,b,c,d]=PLANE_COEF
        if c is not 0:
            z=(-d-a*x1-b*y1)/c
            return [x1,y1,z]
        elif b is not 0:
            y=(-d-a*x1-c*z1)/b
            return [x1,y,z1]
        else:
            x = (-d-b*y1-c*z1)/a
            return [x,y1,z1]
    
    def DistanceToPlane(self,POINT,PLANE_COEF):
        """
        INPUT:
            - POINT: (TUPLE, LIST OR ARRAY) Lenght 3 That represent a single point
            - PLANE_COEF: List of lenght 4 that represent the 4 coeficients [a, b, c, d] of the plane: ax+by+cz+d=0
        Output:
            - Float number that represent the distance from point to plane
        """  
        x1, y1 ,z1 = POINT
        [a,b,c,d]= PLANE_COEF 
        d = abs((a * x1 + b * y1 + c * z1 + d))  
        e = sqrt(a * a + b * b + c * c)
        return d/e 

    

    def PointListDistanceToPlane1(self,POINT_LIST,PLANE_COEF):
        """
        INPUT:
            - POINT_LIST: List of (TUPLE, LIST OR ARRAY) of  Lenght 3 That represent a list of points
            - PLANE_COEF: List of lenght 4 that represent the 4 coeficients [a, b, c, d] of the plane: ax+by+cz+d=0
        Output:
            - list of Float numbers that represent the distance from the points to plane
        """  
        return list(map(lambda x : self.DistanceToPlane(x,PLANE_COEF),POINT_LIST))

    def PointListDistanceToPlane(self, point,coef):
        d=coef[3]
        coef=coef[0:3]
        return (np.dot(point,coef)+d)/np.linalg.norm(coef)

    def StandardDesviation(self,DISTANCE_LIST):
        return sqrt(sum(list(map(lambda x : x**2,DISTANCE_LIST)))/(len(DISTANCE_LIST)-1))


# ===== 2   P L A N E    I N T E R S C T I O N  =======

    def plane_intersect(self,a, b):
        """
        a, b   4-tuples/lists
            Ax + By +Cz + D = 0
            A,B,C,D in order  
        output: 1 point and 1 vector
        """
        a_vec, b_vec = np.array(a[:3]), np.array(b[:3])
        aXb_vec = np.cross(a_vec, b_vec)
        A = np.array([a_vec, b_vec, aXb_vec])
        #print(np.linalg.det(A))
        if np.linalg.det(A) < 0.01: # changed 
            #print("paralell planes")
            return False
        else:
            d = np.array([-a[3], -b[3], 0.]).reshape(3,1)
            # could add np.linalg.det(A) == 0 test to prevent linalg.solve throwing error
            p_inter = np.linalg.solve(A, d).T
            #return p_inter[0], (p_inter + aXb_vec)[0]
            return p_inter[0], aXb_vec


# ===== P L A N E  /  L I N E    I N T E R S E C T I O N ======

    def LinePlaneCollision(self,PLANE_COEF, rayDirection, rayPoint, epsilon=1e-6):
        """
        PLANE_COEF (A,B,C,D) => [1,3,5,6]   for a plane   A*x + B*y + C*z +D =0  =>  x + 3*y + 5*z + 6 = 0
        rayDirection => vector (X,Y,Z) 
        rayPoint => Point (Px,Py,Pz)
        """
        try:
            planeNormal=np.array(PLANE_COEF[:3])
            planePoint=np.array([0,0,(-PLANE_COEF[0]-PLANE_COEF[1]-PLANE_COEF[3])/PLANE_COEF[2]])
            rayPoint=np.array(rayPoint)
            rayDirection=np.array(rayDirection)
            ndotu = planeNormal.dot(rayDirection)
            if abs(ndotu) < epsilon:
                raise RuntimeError("no intersection or line is within plane")
            w = rayPoint - planePoint
            si = -planeNormal.dot(w) / ndotu
            Psi = w + si * rayDirection + planePoint
            return Psi
        except:
            print("No intersection with plane and ray")
            return False

# ===== P O I N T / L I N E   D I S T A N C E ======

    def PointLineDistanceCalcultations(self,POINT1,POINT_LINE,VECTOR_LINE):
        """
        INPUT: 
            POINT1: One point (x,y,z)
            POINT_LINE: One point frome line
            VECTOR_LINE: The vector from line
        OUTPUT:
            NO ERRORS: Tuple with 2 components: (distance, P)
                - distance: nearest distance between point and line
                - P: Point of line that has the nearest distance to selected point
            ERRORS: Boolean =>False Happens when norm of vecto is near 0 or when somthing get whrong in the calculations
        """
        try:
            A=np.array(POINT1)
            B=np.array(POINT_LINE)
            C=B+np.array(VECTOR_LINE)
            coef=np.linalg.norm(C - B)
            if coef>0.00001:
                d=(C - B) / np.linalg.norm(C - B)
                v= A-B
                t=np.dot(d,v)
                P=B + t*d
                distance=np.linalg.norm(P-A)
                return distance, P
            else:
                print("Vector line is to near (0,0,0)")
                return False
        except:
            print("Problems finding distance between Line and point")
            return False

    def PointLineDistance(self,POINT,POINT_LINE,VECTOR_LINE):
        """
        INPUT: 
            POINT:          One point (x,y,z)
            POINT_LINE:     One point frome line
            VECTOR_LINE:    The vector from line
        OUTPUT:
            NO ERRORS:  distance between point and line (FLOAT Number)
            ERRORS:     Boolean =>False Happens when norm of vecto is near 0 or when somthing get whrong in the calculations
        """
        try:
            A=np.array(POINT)
            B=np.array(POINT_LINE)
            C=B+np.array(VECTOR_LINE)
            BA=A-B
            BC=B-C
            return np.linalg.norm(np.cross(BA,BC))/np.linalg.norm(BC)
        except: 
            print("Problems finding distance between Line and point")
            return False
    

    def PointListToLineDistanceCalculation(self,POINT_LIST, POINT_LINE,VECTOR_LINE):
        """

        """
        A=np.array(POINT_LIST)
        B=np.array(POINT_LINE)
        C=B+np.array(VECTOR_LINE)
        BA=A-B
        BC=B-C
        CROSS=np.linalg.norm(np.cross(BA,BC),axis=1)/np.linalg.norm(BC)
        return CROSS

# ====== R A Y  &   R A Y   D I S T A N C E ===============

    def Two3DRayCrossCalcultation(self,POINT1,VECTOR1,POINT2, VECTOR2):
        """
        INPUT: 
            POINT1: One point frome line 1
            VECTOR1: The vector of line 1
            POINT2: One point frome line 2
            VECTOR2: The vector of line 2
        OUTPUT:
            NO ERRORS: Tuple with 4 components: (midpoint, distance, pp1, pp2)
                - midpoint: is the point hows is the nearest between line 1 and 2
                - distance: nearest distance between line 1 and 2
                - pp1: The nearest Point of line 1 to line 2
                - pp2: The nearest Point of line 2 to line 1
            ERRORS: Boolean =>False Happens when Line 1 and line 2 are paralell or when somthing wron in the calculations
        """
        try:
            #L1
            V1=np.array(VECTOR1)
            P1=np.array(POINT1)
            #L2
            V2=np.array(VECTOR2)
            P2=np.array(POINT2)
            #Perpendicular Vector
            V3=np.cross(V1,V2)
            #print(np.linalg.norm(V3))
            if np.linalg.norm(V3)> 0.01:
                # CALCULATIONS
                a=np.array(np.transpose([V1,-V2,V3]))
                b=np.array(P2-P1)
                coef = np.linalg.solve(a, b)
                #Output generating
                pp1= P1 + coef[0]*V1
                pp2=P2 + coef[1]*V2
                distance= np.linalg.norm(pp1-pp2)
                midpoint = (pp1 + pp2)*0.5

                return midpoint, distance, pp1, pp2
            else:
                #print("vector1 and vector2 are parallel")
                return False
        except:
            print("Problems calculating nearest point cros fronm two lines")
            return False

# ===== P O I N T    P L A N E     F I L T E R ========

    def PointListDistanceFilter(self,POINT_LIST,PLANE_COEF,MIN_DISTANCE):
        return list(filter(lambda x : self.DistanceToPlane(x,PLANE_COEF)<MIN_DISTANCE,POINT_LIST))


# ===== P O L Y L I N E   M A N I P U L A T I O N S  =======

    #unit normal vector of plane defined by points a, b, and c
    def unit_normal(self,a, b, c):
        x = np.linalg.det([[1,a[1],a[2]],
            [1,b[1],b[2]],
            [1,c[1],c[2]]])
        y = np.linalg.det([[a[0],1,a[2]],
            [b[0],1,b[2]],
            [c[0],1,c[2]]])
        z = np.linalg.det([[a[0],a[1],1],
            [b[0],b[1],1],
            [c[0],c[1],1]])
        magnitude = (x**2 + y**2 + z**2)**.5
        return (x/magnitude, y/magnitude, z/magnitude)

    #area of polygon 
    def poly_area(self,poly):
        if len(poly) < 3: # not a plane - no area
            return 0
        total = [0, 0, 0]
        N = len(poly)
        for i in range(N):
            vi1 = poly[i]
            vi2 = poly[(i+1) % N]
            prod = np.cross(vi1, vi2)
            total[0] += prod[0]
            total[1] += prod[1]
            total[2] += prod[2]
        result = np.dot(total, self.unit_normal(poly[0], poly[1], poly[2]))
        return abs(result/2)

        
# ===== C L U S T E R I N G   O F   P O I N T C L O U D =========

    def XYZClusterFilter(self, POINT_LIST, X_INTERVAL=None, Y_INTERVAL=None, Z_INTERVAL = None, Count = True, Inbox = True,Outbox = True):
        """
        INPUT:
            POINT_LIST => Array point list POINT_LIST=[[1,3,6],[5,7,9],........]
            X_INTERVAL, Y_INTERVAL, Z_INTERVAL   => Array of lenght 2 or as option None
            InBox, Outbox => Optional boolean input default for both is True
        OUTPUT: 4 VAR TUPLE
            0: Index of list <ARRAY>
            1: Count of indexes <INTEGER> or None
            2: selected point of point cloud <ARRAY> MATRIX 3 x Number of points
            3: unselected points of point cloud <ARRAY> Matrix 3 x Number of points
        """
        try:
            pts1=np.array(POINT_LIST)
            colx=0
            coly=1
            colz=2
            if X_INTERVAL==None:colx=None
            if Y_INTERVAL==None:coly=None
            if Z_INTERVAL==None:colz=None
            columns=list(filter(None.__ne__,[colx,coly,colz]))
            pts=pts1[:,columns]
            ll = np.array(list(filter(None.__ne__,[None if X_INTERVAL==None else X_INTERVAL[0], None if Y_INTERVAL==None else Y_INTERVAL[0],None if Z_INTERVAL==None else Z_INTERVAL[0]])))
            ur = np.array(list(filter(None.__ne__,[None if X_INTERVAL==None else X_INTERVAL[1], None if Y_INTERVAL==None else Y_INTERVAL[1],None if Z_INTERVAL==None else Z_INTERVAL[1]])))
            inidx = np.all(np.logical_and(ll <= pts, pts <= ur), axis=1)
            if Count:
                Count=len(np.where(inidx)[0])
            else:
                Count=None
            if Inbox:
                Inbox = pts1[inidx]
            else:
                Inbox = None
            if Outbox:
                Outbox = pts1[np.logical_not(inidx)]
            else:
                Outbox = None
            return np.where(inidx)[0], Count, Inbox, Outbox
        except:
            print("problems with clustering the PointList")
            return False

# ===== C R E A T E    I N D E X O B J E C T ======

    def CreateIndexObject(self,POINT_LIST,MIN_INTERVAL,MAX_INTERVAL,PARTITIONS, INDEX_LIST = None ,NEARBY_CLUSTERS=True,TIMEOUT=60):
        """
        INPUT:
            POINT_LIST => Array point list POINT_LIST=[[1,3,6],[5,7,9],........]
            MIN_INTERVAL, MAX_INTERVAL => Array of lenght 3 
            PARTITIONS => two options:
                            - Interger => Number of partition on x,y,z axis
                            - Array of length 3 => number of partition on eath of the x,y,z axis 
            INDEX_LIST =>   - None by default
                            - List or array if indexes of the PointList
            NEARBY_CLUSTERS => boolean: - True (default) if output object calculate the nearby clusters of eatch cluster 
                                        - False if output object does not calculate the nearby clusters of eatch cluster 
            TIMEOUT => Number (integer or float). Number of seconds that the function have to generate the index object before it interuped it. default 60 seconds. Seconds 
        OUTPUT:
            Without errors: INDEX OBJECT => To see axample of object ixwcute the function: self.IndexObjectExample()
            with errors: Boolean False

        """
        if isinstance(INDEX_LIST,type(None)):
            index_option=0
        elif isinstance(INDEX_LIST,(list,np.ndarray)):
            if isinstance(INDEX_LIST,list):
                INDEX_LIST=np.array(INDEX_LIST)
            if len(POINT_LIST)==len(INDEX_LIST):
                index_option=1
            else:
                print("Lenght of index list is not the same than lenght of point list")
                return False
        else:
            print("INDEX_LIST is not a list or array")
            return False
        
        
        range_list= np.transpose([MIN_INTERVAL,MAX_INTERVAL]).tolist()
        
        hst, edges, bincounts=binned_statistic_dd(POINT_LIST, values = 1, statistic ='count', bins=np.array(PARTITIONS, dtype=float),range=np.array(range_list, dtype=float),expand_binnumbers=False)
        
        edges_x ,edges_y ,edges_z = edges
        edges_x=np.insert(edges_x,0,0).tolist()
        edges_y=np.insert(edges_y,0,0).tolist()
        edges_z=np.insert(edges_z,0,0).tolist()
        len_edges_x=len(edges_x)
        len_edges_y=len(edges_y)
        len_edges_z=len(edges_z)

        if isinstance(PARTITIONS,int):
            delta_x=(-range_list[0][0]+range_list[0][1])/PARTITIONS
            delta_y=(-range_list[1][0]+range_list[1][1])/PARTITIONS
            delta_z=(-range_list[2][0]+range_list[2][1])/PARTITIONS
        elif isinstance(PARTITIONS,list):
            if len(PARTITIONS)==3:
                delta_x=(-range_list[0][0]+range_list[0][1])/PARTITIONS[0]
                delta_y=(-range_list[1][0]+range_list[1][1])/PARTITIONS[1]
                delta_z=(-range_list[2][0]+range_list[2][1])/PARTITIONS[2]
            else:
                print("Lenght of the partition list is not 3")
                return False
        else:
            print("'PARTITIONS' is not a number or a list ")
            return False
        
        index=0
        ID=0
        index_object=[]
        for ind_x,row_x in enumerate(edges_x):
            for ind_y,row_y in enumerate(edges_y):
                for ind_z,row_z in enumerate(edges_z):
                    if ((ind_x==0) or (ind_y==0) or (ind_z==0) or (ind_x==(len_edges_x-1)) or (ind_y==(len_edges_y-1)) or (ind_z==(len_edges_z-1))):
                        index +=1
                    else:
                        point_index=np.where(bincounts==index)[0]
                        if index_option==0:
                            original_point_index=point_index
                        else:
                            original_point_index=INDEX_LIST[point_index]
                        count=len(point_index)
                        x_interval={"min":row_x,"max":row_x+delta_x}
                        y_interval={"min":row_y,"max":row_y+delta_y}
                        z_interval={"min":row_z,"max":row_z+delta_z}
                        Range={"x_interval":x_interval,"y_interval":y_interval,"z_interval":z_interval}
                        row={"id":ID,"x_id":ind_x-1,"y_id":ind_y-1,"z_id":ind_z-1,"range":Range,"point_index":point_index,"original_point_index":original_point_index,"point_count":count,"nearby":None}
                        index_object.append(row)
                        ID +=1
                        index +=1
        #print(row_object)
        if NEARBY_CLUSTERS :
            index=0
            for row in index_object:
                x_id=row["x_id"]
                y_id=row["y_id"]
                z_id=row["z_id"]
                id_list=[]
                for newrow in index_object:  
                    new_x_id=newrow["x_id"]
                    new_y_id=newrow["y_id"]
                    new_z_id=newrow["z_id"]
                    if (new_x_id==(x_id+1) or new_x_id==(x_id-1) or new_x_id==(x_id)) and (new_y_id==(y_id+1) or new_y_id==(y_id-1) or new_y_id==(y_id) ) and (new_z_id==(z_id+1) or new_z_id==(z_id-1) or new_z_id==(z_id)):
                        id_list.append(newrow["id"])
                index_object[index]["nearby"]=id_list
                index += 1
        return index_object

    def CreateIndexObject1(self,POINT_LIST,MIN_INTERVAL,MAX_INTERVAL,PARTITIONS, INDEX_LIST = None ,NEARBY_CLUSTERS=True,TIMEOUT=60):
        """
        INPUT:
            POINT_LIST => Array point list POINT_LIST=[[1,3,6],[5,7,9],........]
            MIN_INTERVAL, MAX_INTERVAL => Array of lenght 3 
            PARTITIONS => two options:
                            - Interger => Number of partition on x,y,z axis
                            - Array of length 3 => number of partition on eath of the x,y,z axis 
            INDEX_LIST =>   - None by default
                            - List or array if indexes of the PointList
            NEARBY_CLUSTERS => boolean: - True (default) if output object calculate the nearby clusters of eatch cluster 
                                        - False if output object does not calculate the nearby clusters of eatch cluster 
            TIMEOUT => Number (integer or float). Number of seconds that the function have to generate the index object before it interuped it. default 60 seconds. Seconds 
        OUTPUT:
            Without errors: INDEX OBJECT => To see axample of object ixwcute the function: self.IndexObjectExample()
            with errors: Boolean False
        """
        try:
            if isinstance(INDEX_LIST,type(None)):
                index_option=0
            elif isinstance(INDEX_LIST,(list,np.ndarray)):
                if isinstance(INDEX_LIST,list):
                    INDEX_LIST=np.array(INDEX_LIST)
                if len(POINT_LIST)==len(INDEX_LIST):
                    index_option=1
                else:
                    print("Lenght of index list is not the same than lenght of point list")
                    return False
            else:
                print("INDEX_LIST is not a list or array")
                return False

            timeout = time.time() + TIMEOUT
            index_object=[]
            ID=0
            [min_x,min_y,min_z]=MIN_INTERVAL
            [max_x,max_y,max_z]=MAX_INTERVAL

            if isinstance(PARTITIONS,int):           
                PARTITIONS=PARTITIONS+1
                x_intervals=np.linspace(min_x,max_x,PARTITIONS)
                y_intervals=np.linspace(min_y,max_y,PARTITIONS)
                z_intervals=np.linspace(min_z,max_z,PARTITIONS)
            elif isinstance(PARTITIONS,list):
                if len(PARTITIONS)==3:
                    x_intervals=np.linspace(min_x,max_x,PARTITIONS[0]+1)
                    y_intervals=np.linspace(min_y,max_y,PARTITIONS[1]+1)
                    z_intervals=np.linspace(min_z,max_z,PARTITIONS[2]+1)
                else:
                    print("Lenght of the partition list is not 3")
                    return False
            else:
                print("'PARTITIONS' is not a number or a list ")
                return False

            len_int_x=len(x_intervals)
            len_int_y=len(y_intervals)
            len_int_z=len(z_intervals)

            for i in range(0,len_int_x-1):
                for j in range(0,len_int_y-1):
                    for k in range(0,len_int_z-1):
                        x_min=x_intervals[i]
                        y_min=y_intervals[j]
                        z_min=z_intervals[k]
                        x_max=x_intervals[i+1]
                        y_max=y_intervals[j+1]
                        z_max=z_intervals[k+1]
                        x_interval={"min":x_min,"max":x_max}
                        y_interval={"min":y_min,"max":y_max}
                        z_interval={"min":z_min,"max":z_max}
                        Range={"x_interval":x_interval,"y_interval":y_interval,"z_interval":z_interval}
                        index, Count, Inbox, Outbox = self.XYZClusterFilter( POINT_LIST, X_INTERVAL=[x_min,x_max], Y_INTERVAL=[y_min,y_max], Z_INTERVAL = [z_min,z_max], Count = True, Inbox = False,Outbox = False)
                        if index_option==0:
                            original_point_index=index
                        else:
                            original_point_index=INDEX_LIST[index]
                        row={"id":ID,"x_id":i,"y_id":j,"z_id":k,"range":Range,"point_index":index,"original_point_index":original_point_index,"point_count":Count,"nearby":None}
                        index_object.append(row)
                        ID +=1
                        if time.time() > timeout:
                            print("timeout error")
                            return False
            
            if NEARBY_CLUSTERS :
                index=0
                for row in index_object:
                    x_id=row["x_id"]
                    y_id=row["y_id"]
                    z_id=row["z_id"]
                    id_list=[]
                    for newrow in index_object:  
                        new_x_id=newrow["x_id"]
                        new_y_id=newrow["y_id"]
                        new_z_id=newrow["z_id"]
                        if (new_x_id==(x_id+1) or new_x_id==(x_id-1) or new_x_id==(x_id)) and (new_y_id==(y_id+1) or new_y_id==(y_id-1) or new_y_id==(y_id) ) and (new_z_id==(z_id+1) or new_z_id==(z_id-1) or new_z_id==(z_id)):
                            id_list.append(newrow["id"])
                    index_object[index]["nearby"]=id_list
                    index += 1
            return index_object
        except:
            print("Error making the index object")
            return False

    def IndexObjectExample(self):
        return [{'range': {'x_interval': {'max': 1.0, 'min': 0.0}, 'y_interval': {'max': 1.0, 'min': 0.0}, 'z_interval': {'max': 0.3333333333333333, 'min': 0.0}}, 'id': 0, 'point_cound': 5, 'nearby': [0, 1], 'point_index': [ 1,  3, 11, 12, 18], 'z_id': 0, 'y_id': 0, 'x_id': 0}, {'range': {'x_interval': {'max': 1.0, 'min': 0.0}, 'y_interval': {'max': 1.0, 'min': 0.0}, 'z_interval': {'max': 0.6666666666666666, 'min': 0.3333333333333333}}, 'id': 1, 'point_cound': 6, 'nearby': [0, 1, 2], 'point_index': [ 0,  6,  8, 13, 14, 17], 'z_id': 1, 'y_id': 0, 'x_id': 0}, {'range': {'x_interval': {'max': 1.0, 'min': 0.0}, 'y_interval': {'max': 1.0, 'min': 0.0}, 'z_interval': {'max': 1.0, 'min': 0.6666666666666666}}, 'id': 2, 'point_cound': 9, 'nearby': [1, 2], 'point_index': [ 2,  4,  5,  7,  9, 10, 15, 16, 19], 'z_id': 2, 'y_id': 0, 'x_id': 0}]

# ===== I N D E X O B J E C T    O P E R A T I O N S ======

    def IndexObjectRowInRangeQ(self,POINT,INDEX_OBJECT_ROW):
        Range=INDEX_OBJECT_ROW["range"]
        x_interval=Range["x_interval"]
        y_interval=Range["y_interval"]
        z_interval=Range["z_interval"]
        if ((x_interval["max"]>=POINT[0]>=x_interval["min"]) and (y_interval["max"]>=POINT[1]>=y_interval["min"]) and (z_interval["max"]>=POINT[2]>=z_interval["min"])):
            return True
        else:
            return False



    def GetIndexRowFromCoordenates(self, COORDENATE, INDEX_OBJECT):
        return list(filter(lambda x : self.IndexObjectRowInRangeQ(COORDENATE,x),INDEX_OBJECT))

    def GetIndexRowFromID(self,ID,INDEX_OBJECT):
        return list(filter(lambda x : x["id"]==ID,INDEX_OBJECT))



    def GetIDFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["id"]
    
    def GetX_IDFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["x_id"]

    def GetY_IDFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["y_id"]
    
    def GetZ_IDFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["z_id"]

    def GetXYZ_ID_ListFromIndexRow(self, INDEX_OBJECT_ROW):
        return [INDEX_OBJECT_ROW["x_id"],INDEX_OBJECT_ROW["y_id"],INDEX_OBJECT_ROW["z_id"]]

    def GetPointIndexFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["point_index"]

    def GetOriginalPointIndexFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["original_point_index"]

    def GetNearbyIDListFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["nearby"]
    
    def GetPointCountFromIndexRow(self, INDEX_OBJECT_ROW):
        return INDEX_OBJECT_ROW["point_count"]

    def GetX_IntervalFromIndexRow(self,INDEX_OBJECT_ROW):
        Range=INDEX_OBJECT_ROW["range"]
        interval=Range["x_interval"]
        return [interval["min"],interval["max"]]
    
    def GetY_IntervalFromIndexRow(self,INDEX_OBJECT_ROW):
        Range=INDEX_OBJECT_ROW["range"]
        interval=Range["y_interval"]
        return [interval["min"],interval["max"]]
    
    def GetZ_IntervalFromIndexRow(self,INDEX_OBJECT_ROW):
        Range=INDEX_OBJECT_ROW["range"]
        interval=Range["z_interval"]
        return [interval["min"],interval["max"]]

    def GetIntervalFromIndexRow(self,INDEX_OBJECT_ROW):
        Range=INDEX_OBJECT_ROW["range"]
        x_interval=Range["x_interval"]
        y_interval=Range["y_interval"]
        z_interval=Range["z_interval"]
        return np.transpose([[x_interval["min"],x_interval["max"]],[y_interval["min"],y_interval["max"]],[z_interval["min"],z_interval["max"]]])

    def GetPointListFromIndexRow(self, POINT_LIST, INDEX_OBJECT_ROW):
        point_list = POINT_LIST[INDEX_OBJECT_ROW["point_index"]]
        return point_list

    def GetPointListByIDFromIndexObject(self, ID, POINT_LIST, INDEX_OBJECT,NEARBY_POINTS=False):
        try:
            index_row = self.GetIndexRowFromID(ID,INDEX_OBJECT)[0]
            point_list = self.GetPointListFromIndexRow(POINT_LIST, index_row)
            if NEARBY_POINTS:
                index_list=self.GetNearbyIDListFromIndexRow(index_row)
                nearby_points=[point_list]
                for nearby_ID in index_list:
                    if nearby_ID==ID:
                        continue
                    else:
                        #nearby_index_row=self.GetIndexRowFromID(nearby_ID,INDEX_OBJECT)[0]
                        #nearby_point_list = self.GetPointListFromIndexRow(POINT_LIST, nearby_index_row)
                        nearby_point_list=self.GetPointListByIDFromIndexObject(nearby_ID, POINT_LIST, INDEX_OBJECT,NEARBY_POINTS=False)
                        nearby_points.append(nearby_point_list)
                return np.concatenate(nearby_points)
            else:
                return point_list
        except:
            print("ERROR get point list")
            return False
        

if __name__ == '__main__':
    pass