
from mstm_studio.mstm_spectrum import SPR, Material
import os   # file path operations
import datetime
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
try:
    from miepython.field import eh_near_cartesian
except ImportError:
    print('Mie theory is disabled. Please install `miepython` package')
    pass


class NearFieldMie(object):
    '''
    Calculate |E|^2 field distribution on a 2D map
    for single sphere within Mie theory
    '''
    def calculate(self,
                  wavelength=550,
                  material=1.5+0.1j,
                  environment_material=1.,
                  radius=40,
                  plane='zx', hmin=-10., hmax=10.,
                  vmin=-10., vmax=10., step=1.,
                  include_incident=True):
        '''
        '''
        self.plane = plane
        if isinstance(material, Material):
            material = material.get_nk([wavelength])[0]
        self.h = np.arange(hmin, hmax + step / 2., step)
        self.v = np.arange(vmin, vmax + step / 2., step)
        plane = plane.upper()
        if plane == 'ZX':
            Z, X = np.meshgrid(self.h, self.v, indexing='xy')
            Y = np.zeros_like(X)
        elif plane == 'XZ':
            X, Z = np.meshgrid(self.h, self.v, indexing='xy')
            Y = np.zeros_like(X)
        elif plane == 'ZY':
            Z, Y = np.meshgrid(self.h, self.v, indexing='xy')
            X = np.zeros_like(Y)
        elif plane == 'YZ':
            Y, Z = np.meshgrid(self.h, self.v, indexing='xy')
            X = np.zeros_like(Y)
        elif plane == 'XY':
            X, Y = np.meshgrid(self.h, self.v, indexing='xy')
            Z = np.zeros_like(Y)
        elif plane == 'YX':
            Y, X = np.meshgrid(self.h, self.v, indexing='xy')
            Z = np.zeros_like(Y)
        else:
            raise Exception(f'Wrong plane: {plane}')

        self.E_xyz, self.H_xyz = eh_near_cartesian(
            lambda0=wavelength,  # Vacuum wavelength
            d_sphere=2*radius,  # Sphere diameter
            m_sphere=np.conj(material),  # Sphere refractive index
            n_env=environment_material,  # Refractive index of the surrounding medium
            x=X, y=Y, z=Z,  # cartesian coordinates where to calculate
            include_incident=include_incident  # Include incident field for points outside sphere
        )
        return np.sum(np.abs(self.E_xyz)**2, axis=0)


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
    # ~ from matplotlib.patches import Circle
    wl = 240
    # ~ matsph = Material(os.path.join('nk', 'etaSilver.txt'))
    matsph = 0.5 + 0.1j
    matrix = 1.5
    hmin, hmax, vmin, vmax, step = -20, 20, -20, 20, 0.25
    # ~ hmin, hmax, vmin, vmax, step = -60, 60, -60, 60, 30
    a = 10

    # mie python
    nf = NearFieldMie()
    e2 = nf.calculate(wavelength=wl,
                      material=matsph,
                      environment_material=matrix,
                      radius=a,
                      plane='xz', hmin=hmin, hmax=hmax,
                      vmin=vmin, vmax=vmax, step=step,
                      include_incident=True)
    fig, ax = plt.subplots(1, 1)  #, figsize=(5, 6))
    im = ax.pcolormesh(nf.h, nf.v, e2, cmap='hot', shading='auto')
    ax.set_aspect('equal')
    # ~ ax.add_patch(Circle((0.0, 0.0), a,
                 # ~ fill=False, color='white', lw=1.2))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    caxs = fig.add_axes([0.9, 0.1, 0.05, 0.8])  # left, bottom, width, height
    fig.colorbar(im, cax=caxs, orientation='vertical')
    plt.tight_layout()
    plt.savefig('nf_miepython.png')
    plt.show()

    # mstm v.3
    nf = NearField(wavelength=wl)
    nf.environment_material = matrix
    nf.set_incident_field(fixed=True, azimuth_angle=0.0,
                          polar_angle=0.0, polarization_angle=0.0)
    nf.set_plane(plane='xz', hmin=hmin, hmax=hmax,
                 vmin=vmin, vmax=vmax, step=step)

    spheres = ExplicitSpheres(1, [0, 0, 0, a],
                              mat_filename=Material(matsph))
    # ~ spheres = ExplicitSpheres(2, [-6, 0, 0, 5, 6, 0, 0, 5],
                          # ~ mat_filename=2*[Material('nk/etaGold.txt')])
    nf.set_spheres(spheres)
    nf.simulate()
    fig, ax = plt.subplots(1, 1, figsize=(5, 6))
    nf.plot(fig=fig, axs=ax)
    plt.tight_layout()
    plt.savefig('nf_mstmv3.png')
    plt.show()
    # ~ nf.write('nearfield.dat')

    print('See you!')
