from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math, random


class VDW:
    '''
    Calculate vdw surface any molecule
    
    vdw radii from www.ccdc.cam.ac.uk
    
    '''
    vdw_radii = {'Ac': 2.0, 'Al': 2.0, 'Am': 2.0, 'Sb': 2.0, 'Ar': 1.88, 'As': 1.85, 'At': 2.0, 'Ba': 2.0, 'Bk': 2.0, 'Be': 2.0, 'Bi': 2.0, 'Bh': 2.0, 'B': 2.0, 'Br': 1.85,
       'Cd': 1.58, 'Cs': 2.0, 'Ca': 2.0, 'Cf': 2.0, 'C': 1.7, 'Ce': 2.0, 'Cl': 1.75, 'Cr': 2.0, 'Co': 2.0, 'Cu': 1.4, 'Cm': 2.0, 'Ds': 2.0, 'Db': 2.0, 'Dy': 2.0,
       'Es': 2.0, 'Er': 2.0, 'Eu': 2.0, 'Fm': 2.0, 'F': 1.47, 'Fr': 2.0, 'Gd': 2.0, 'Ga': 1.87, 'Ge': 2.0, 'Au': 1.66, 'Hf': 2.0, 'Hs': 2.0, 'He': 1.4, 'Ho': 2.0,
       'H': 1.2, 'In': 1.93, 'I': 1.98, 'Ir': 2.0, 'Fe': 2.0, 'Kr': 2.02, 'La': 2.0, 'Lr': 2.0, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2.0, 'Mg': 1.73, 'Mn': 2.0, 'Mt': 2.0,
       'Md': 2.0, 'Hg': 1.55, 'Mo': 2.0, 'Nd': 2.0, 'Ne': 1.54, 'Np': 2.0, 'Ni': 1.63, 'Nb': 2.0, 'N': 1.55, 'No': 2.0, 'Os': 2.0, 'O': 1.52, 'Pd': 1.63, 'P': 1.8,
       'Pt': 1.72, 'Pu': 2.0, 'Po': 2.0, 'K': 2.75, 'Pr': 2.0, 'Pm': 2.0, 'Pa': 2.0, 'Ra': 2.0, 'Rn': 2.0, 'Re': 2.0, 'Rh': 2.0, 'Rb': 2.0, 'Ru': 2.0, 'Rf': 2.0,
       'Sm': 2.0, 'Sc': 2.0, 'Sg': 2.0, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27, 'Sr': 2.0, 'S': 1.8, 'Ta': 2.0, 'Tc': 2.0, 'Te': 2.06, 'Tb': 2.0, 'Tl': 1.96,
       'Th': 2.0, 'Tm': 2.0, 'Sn': 2.17, 'Ti': 2.0, 'W': 2.0, 'U': 1.86, 'V': 2.0, 'Xe': 2.16, 'Yb': 2.0, 'Y': 2.0, 'Zn': 1.39, 'Zr': 2.0}

    def __init__(self,coordinates,points = 80):
        '''
        Args:
            coordinates = [[atom_type,x,y,z],..]
            points = how much points per one atom (80 by default)
    
        '''
        
        self.coordinates = coordinates
        self.number_of_atoms = len(coordinates)
        self.number_of_points = self.number_of_atoms * points

        self.add_radii()


    def add_radii(self):
        '''
        Returns:
            coordinates = [[radius for atom_type,x,y,z],..]
    
        '''
        
        for atom in self.coordinates:
            for k,v in self.vdw_radii.items():
                if atom[0] == k:
                        atom[0] = v
        pass

        
    @classmethod
    def get_from_xyz(cls,file_name):
        '''
        Returns:
            coordinates = [[atom_type,x,y,z],..] from xyz. file
    
        '''        
        coordinates = []
        with open(file_name,'r') as inp:
            lines = inp.readlines()
            for line in lines:
                line = line.split()
                if len(line) == 4:
                    coordinates.append([line[0],float(line[1]),float(line[2]),float(line[3])])
        return cls(coordinates)


    @classmethod
    def get_from_cif(cls,file_name):
        pass



    def generate_point(self, xyz_point):
        '''
        generate point on surface of sphere

        '''
        
        radius,x,y,z = xyz_point
        
        theta = random.random()*2*math.pi
        phi = random.random()*math.pi

        sin_phi = math.sin(phi)
        x_new = radius*math.cos(theta)*sin_phi + x
        y_new = radius*math.sin(theta)*sin_phi + y
        z_new = radius*math.cos(phi) + z
        
        return [x_new,y_new,z_new]



    def lenght(self,one_point, other_point):
        '''
        Returns:
            distance between two points
    
        '''
        
        vector = list(map(lambda i,j: j-i, one_point, other_point))
        return (sum(map(lambda i: i**2, vector)))**(0.5)



    def create_vdw_surface(self):
        '''
        Returns:
            points on vdw surface = [[x,y,z],..]
    
        '''
        points = list()
        while len(points) != self.number_of_points:
            
            atom = random.choice(self.coordinates)
            point = self.generate_point(atom)
            l = []
            for i in self.coordinates:
                if i != atom:
                    diff = self.lenght(point,i[1:]) - i[0]
                    l.append(diff)
            if all(i >= 0 for i in l):
                    points.append(point)
                    
        return points



    def show_surface(self, points):
        '''
        Returns:
            3d plot of surface
    
        '''
        x = [i[0] for i in points]
        y = [i[1] for i in points]
        z = [i[2] for i in points]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.axis('equal')

        plt.show()


    def save_xyz(self,points, nameoffile = 'pointcloud.xyz'):
        '''
        Returns:
            save xyz coordinates of points of vdw surface to XYZ. file
    
        '''
        
        points = [[str(i) for i in j] for j in points]
        
        with open(nameoffile,'r') as f: #delete content if file exists
            pass
        
        with open(nameoffile,'a') as f:
            for i in points:
                f.write('   '.join(i)+'\n')

