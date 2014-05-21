import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
import cartopy.mpl.geoaxes as geo_ax
 
import iris
from iris.analysis.interpolate import regrid
from iris.cube import Cube

class VectorField(object):
    """
    A collection of vector components which together represent a single
    vector field.

    """


class VectorField2D(VectorField):
    """
    A collection of vector components which together represent a single,
    two dimensional vector field.

    """

    def __init__(self, i_comp=None, j_comp=None, magnitude=None,
                 direction=None):
        self.components = []
        if isinstance(i_comp, Cube):
            self.i_component = i_comp
            self.components.append(self.i_component)
        else:
            if i_comp is not None:
                raise TypeError('VectorField components must be Iris Cubes\n'
                                'not {}'.format(type(i_comp)))
        if isinstance(j_comp, Cube):
            self.j_component = j_comp
            self.components.append(self.j_component)
        else:
            if j_comp is not None:
                raise TypeError('VectorField components must be Iris Cubes\n'
                                'not {}'.format(type(j_comp)))
        if hasattr(self, 'i_component') and hasattr(self, 'j_component'):
            self._regrid_i_j()
        if isinstance(magnitude, Cube):
            self.magnitude = mag
            self.components.append(self.magnitude)
        else:
            if magnitude is not None:
                raise TypeError('VectorField components must be Iris Cubes\n'
                                'not {}'.format(type(magnitude)))
            else:
                self.magnitude = self._magnitude()
        if isinstance(direction, Cube):
            self.direction = direction
            self.components.append(self.direction)
        else:
            if direction is not None:
                raise TypeError('VectorField components must be Iris Cubes\n'
                                'not {}'.format(type(direction)))
            else:
                self.direction = self._direction()
        # catch multiple crs inconsistency error
        _crs = self.crs()
        
        
    def crs(self):
        """
        Returns the horizontal geospatial coordinate reference system the
        vector field dimension coordiantes are defined with respect to.

        Returns a ValueError if all components do not define the same dimension
        coordinate CRSs.
        """
        crs_set = []
        for cube in self.components:
            for dim_coord in cube.dim_coords:
                if dim_coord.coord_system:
                    if dim_coord.coord_system not in crs_set:
                        crs_set.append(dim_coord.coord_system)
        if len(crs_set) == 0:
            res = None
        elif len(crs_set) != 1:
            raise ValueError('Vector Field components are defined with respect'
                             ' to different coordinate reference systems')
        else:
            res = crs_set.pop()
        return res

    def _components_colocate(self):
        """
        Returns True if all components colocate, else returns False.
        
        """
        result = True
        for cube in self.components:
           cube_result = cube.coords() == self.components[0].coords()
           if not cube_result:
               result = False
        return result

    def _c_stagger_common_grid(self):
        """
        return a grid cube which may be used to resample non-colocated
        components

        factored for Arakawa C grid only
        """
        x = self.i_component.coord(axis='x') != self.j_component.coord(axis='x')
        y = self.i_component.coord(axis='y') != self.j_component.coord(axis='y')
        if x and y:
            grid_cube = self.i_component.copy()
            cdim = grid_cube.coord_dims(grid_cube.coord(axis='x'))
            if len(cdim) == 1:
                cdim = cdim[0]
            else:
                raise ValueError('X coordiante spans more than 1 dimension')
            grid_cube.remove_coord(grid_cube.coord(axis='x'))
            grid_cube.add_dim_coord(self.j_component.coord(axis='x').copy()
                                    , cdim)
            res = grid_cube
        else:
            res = None
        return res
            

    def _regrid_i_j(self):
        if self._components_colocate():
            self.coloc_i_component = self.i_component
            self.coloc_j_component = self.j_component
        grid_cube = self._c_stagger_common_grid()
        if grid_cube and not hasattr(self, 'coloc_i_component') and \
            not hasattr(self, 'coloc_j_component'):
            self.coloc_i_component = regrid(self.i_component, grid_cube)
            self.coloc_j_component = regrid(self.j_component, grid_cube)
            
        

    def _magnitude(self):
        """
        returns magnitude if available, or calculates it using components and
        linear interpolation if required
        
        """
        if hasattr(self, 'magnitude'):
            res = self.magnitude
        else:
            res = (self.coloc_i_component ** 2 + self.coloc_j_component **2)**0.5
        return res

    def _direction(self):
        if hasattr(self, 'direction'):
            res = self.direction
        else:
            res = self.coloc_i_component.copy()
            res.standard_name = None
            res.long_name = 'direction (from North)'
            res.units = iris.unit.Unit('radians')
            res.data = np.arctan(self.coloc_i_component.data / self.coloc_j_component.data)
        return res

    def streamplot(self, projection=None):
        """
        produce a streamplot of the vector field

        """
        if projection is None:
            projection = self.crs().as_cartopy_projection()

        ax = plt.axes(projection=projection)
        ax.set_extent([self.coloc_i_component.coord(axis='x').points.min(),
                       self.coloc_i_component.coord(axis='x').points.max(),
                       self.coloc_i_component.coord(axis='y').points.min(),
                       self.coloc_i_component.coord(axis='y').points.max()],
                self.crs().as_cartopy_projection())
        ax.coastlines()

        ax.streamplot(self.coloc_i_component.coord(axis='x').points,
                      self.coloc_i_component.coord(axis='y').points,
                      self.coloc_i_component.data,
                      self.coloc_j_component.data,
                      transform=self.crs().as_cartopy_projection(),
                      linewidth=2, density=2, color=self.magnitude.data)
        plt.show()


def main():
    fpath = '/data/local/dataZoo/FF/ukV/ukvtuv.T+0'
    xwind = iris.load_cube(fpath, 'x_wind')
    ywind = iris.load_cube(fpath, 'y_wind')
    wind = VectorField2D(i_comp=xwind, j_comp=ywind)
    wind.streamplot(ccrs.PlateCarree())
    # bug in vector transformation code??
    # wind.streamplot()
    

if __name__ == '__main__':
    main()
    
    
