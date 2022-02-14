import sys
sys.path.append('.')
from geometricManipulation import GeometricManipulation

if __name__ == '__main__':
    
    pr = GeometricManipulation()
    result = pr.Absolute3DDistance([1,4,6.7],[2,2,2])
    print(result)

    a=pr.CombinationOfTransformationMatrix3x3([[1,2,3],[5,6,7],[8,9,10]],[[1,0,0],[1,1,0],[0,0,1]])
    print(a)
