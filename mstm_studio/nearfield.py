
from mstm_studio.mstm_spectrum import SPR, Material, SpheresOverlapError
import os   # file path operations
import datetime
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


class NearField(SPR):
    '''
    Calculate field distribution map at fixed wavelength
    '''

    def __init__(self, wavelength, mstm_path='~/bin/mstm.x',
                 environment_material='Air', temp_dir=None):
        super().__init__([wavelength], mstm_path,
                         environment_material, temp_dir)
        self.set_incident_field(fixed=True, azimuth_angle=0.0,
                                polar_angle=0.0, polarization_angle=0.0)
        self.set_plane()
        self.paramDict['calculate_near_field'] = 1  # do it

    def set_plane(self, plane='zx', hmin=-10., hmax=10.,
                  vmin=-10., vmax=10., step=1., offset=0.):
        '''
        Determine the plane and grid for near field computation.

        plane: one of 'yz'|'zx'|'xy'
        hmin, hmax, vmin, vmax: horizontal and vertical sizes
        step:   size of the grid grain
        offset: shift of the plane
        '''
        self.hmin = hmin
        self.hmax = hmax
        self.vmin = vmin
        self.vmax = vmax
        self.step = step
        self.offset = offset
        self.nh = int(np.round((hmax - hmin) / step)) + 1
        self.nv = int(np.round((vmax - vmin) / step)) + 1
        print('Field computation grid: %ix%i' % (self.nh, self.nv))
        if plane == 'zy':
            plane = 'yz'
        if plane == 'xz':
            plane = 'zx'
        if plane == 'yx':
            plane = 'xy'
        # 1: y - z plane;  2: z - x plane;  3: x - y
        if plane == 'yz':
            self.paramDict['near_field_plane_coord'] = 1
        elif plane == 'zx':
            self.paramDict['near_field_plane_coord'] = 2
        elif plane == 'xy':
            self.paramDict['near_field_plane_coord'] = 3
        else:
            raise Exception('Wrong plane specification! \n %s' % plane)
        # TODO: extend to calculate complex E vector = 1,
        #       complex E and H vectors = 2
        self.paramDict['near_field_output_data'] = 0  # |E|^2
        self.paramDict['near_field_output_file'] = 'nf-temp.dat'

        k = 2.0 * np.pi / self.wavelengths[0]
        self.paramDict['spacial_step_size'] = k * self.step
        self.paramDict['near_field_plane_position'] = k * self.offset
        self.paramDict['near_field_plane_vertices'] = [k * self.hmin,
                                                       k * self.vmin,
                                                       k * self.hmax,
                                                       k * self.vmax]

    def _read_output(self, tmpdir):
        ''' read nearfield spatial distribution
            from file specified in `near_field_output_file`
            parameter
        '''
        fn = os.path.join(tmpdir,
                          self.paramDict['near_field_output_file'])
        with open(fn) as fout:
            fout.readline()  # skip 1st
            nsph = int(fout.readline().strip())  # no. of spheres in plane
        data = np.loadtxt(fn, skiprows=2 + nsph)
        self.field = np.reshape(data[:, 2], [self.nh, self.nv])
        return self.field

    def _get_grid_hv(self):
        h = np.arange(self.hmin, self.hmax + self.step / 2., self.step)
        v = np.arange(self.vmin, self.vmax + self.step / 2., self.step)
        return h, v

    def write(self, filename):
        ''' save field data to text file'''
        h, v = self._get_grid_hv()
        with open(filename, 'w') as fout:
            if self.paramDict['near_field_plane_coord'] == 1:
                fout.write('# Y[nm]\tZ[nm]\t|E|^2\r\n')
            elif self.paramDict['near_field_plane_coord'] == 2:
                fout.write('# Z[nm]\tX[nm]\t|E|^2\r\n')
            elif self.paramDict['near_field_plane_coord'] == 3:
                fout.write('# X[nm]\tY[nm]\t|E|^2\r\n')

            for i, x in enumerate(h):
                for j, y in enumerate(v):
                    fout.write('%.4f\t%.4f\t%.8f\r\n' % (x, y,
                                                         self.field[i, j]))

    def plot(self, fig=None, axs=None, caxs=None):
        '''
        Show 2D field distribution

        Parameters:

        fig:
            matplotlib figure
        axs:
            matplotlib axes
        caxs:
            matplotlib axes for colorbar

        Returns:
            filled/created fig and axs objects
        '''
        x, y = self._get_grid_hv()
        xx, yy = np.meshgrid(x, y)
        xx = np.transpose(xx)
        yy = np.transpose(yy)
        zz = self.field
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        im = axs.pcolormesh(xx, yy, zz, cmap='hot', shading='auto')
        if caxs is None:
            caxs = fig.add_axes([0.9, 0.1, 0.05, 0.8])  # left, bottom, width, height
        else:
            caxs.clear()
        fig.colorbar(im, cax=caxs, orientation='vertical')
        if self.paramDict['near_field_plane_coord'] == 1:
            axs.set_xlabel('Y, nm')
            axs.set_ylabel('Z, nm')
        elif self.paramDict['near_field_plane_coord'] == 2:
            axs.set_xlabel('Z, nm')
            axs.set_ylabel('X, nm')
        elif self.paramDict['near_field_plane_coord'] == 3:
            axs.set_xlabel('X, nm')
            axs.set_ylabel('Y, nm')
        axs.set_aspect('equal')
        if flag:
            plt.show()
        return fig, axs


if __name__ == '__main__':
    from mstm_studio.mstm_spectrum import Material, ExplicitSpheres
    mat1 = Material(os.path.join('nk', 'etaSilver.txt'))
    nf = NearField(wavelength=340, temp_dir='./temp/')
    nf.environment_material = 'glass'
    nf.set_plane(plane='xz', hmin=-10, hmax=20, vmin=-15, vmax=15, step=0.25)

    spheres = ExplicitSpheres(2, [0, 0, 0, 5, 0, 0, 11, 3],
                              mat_filename=2*[mat1])
    nf.set_spheres(spheres)
    nf.simulate()
    nf.plot()
    nf.write('nearfield.dat')

    print('See you!')
