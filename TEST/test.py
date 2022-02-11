import sys
sys.path.append('.')
from geometricManipulation import GeometricManipulation

if __name__ == '__main__':
    
    pr = GeometricManipulation()
    result = pr.Absolute3DDistance([1,4,6.7],[2,2,2])
    print(result)