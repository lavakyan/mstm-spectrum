
from mstm_studio.mstm_spectrum import SPR, Material, SpheresOverlapError
import subprocess
import os   # to delete files after calc.
import sys  # to check whether running on Linux or Windows
import tempfile
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

    def __init__(self, wavelength):
        super().__init__([wavelength])
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

        self.hmin = hmin
        self.hmax = hmax
        self.vmin = vmin
        self.vmax = vmax
        self.step = step
        self.offset = offset
        self.nh = int(np.round((hmax - hmin) / step)) + 1
        self.nv = int(np.round((vmax - vmin) / step)) + 1
        print('Field computation grid: %ix%i' % (self.nh, self.nv))

    def simulate(self):
        if self.paramDict['number_spheres'] == 0:  # np spheres
            return self.wavelengths, np.zeros_like(self.wavelengths)
        if self.spheres.check_overlap():
            raise SpheresOverlapError('Spheres overlapping!')
        if isinstance(self.environment_material, Material):
            material = self.environment_material
        else:
            print(self.environment_material)
            material = Material(self.environment_material)
        with tempfile.TemporaryDirectory() as tmpdir:  # COPIED FROM PARENT!
            print('Using temporary directory: %s' % tmpdir)
            outFID = open(os.path.join(tmpdir, 'scriptParams.inp'), 'w')
            outFID.write('begin_comment\n')
            outFID.write('**********************************\n')
            outFID.write('  MSTM input for NearField calculation\n')
            outFID.write('  Generated by python script\n')
            outFID.write('  %s\n' %
                         datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
            outFID.write('**********************************\n')
            outFID.write('end_comment\n')
            for key in self.paramDict.keys():
                if key not in self.local_keys:
                    outFID.write(key + '\n')
                    if isinstance(self.paramDict[key], str):
                        svalue = self.paramDict[key]
                    else:
                        if isinstance(self.paramDict[key], list):
                            svalue = '  '.join(map(str, self.paramDict[key]))
                        else:
                            svalue = str(self.paramDict[key])
                        # replace exponent symbol
                        svalue = svalue.replace('e', 'd', 1)
                    outFID.write('%s \n' % svalue)

            lam = self.wavelengths[0]
            k = 2.0 * 3.14159 / lam
            outFID.write('begin_comment\n')
            outFID.write('**********************************\n')
            outFID.write('  Wavelength  %.3f \n' % lam)
            outFID.write('**********************************\n')
            outFID.write('end_comment\n')
            outFID.write('output_file\n')
            outFID.write('mstm_l%.0f.out\n' % (lam * 1000))
            outFID.write('length_scale_factor\n')
            outFID.write('  %.6f\n' % k)
            outFID.write('medium_real_ref_index\n')
            outFID.write('  %f\n' % material.get_n(lam))
            outFID.write('medium_imag_ref_index\n')
            outFID.write('  %f\n' % material.get_k(lam))

            outFID.write('spacial_step_size\n')
            outFID.write('  %f\n' % (k * self.step))
            outFID.write('near_field_plane_position\n')
            outFID.write('  %f\n' % (k * self.offset))
            outFID.write('near_field_plane_vertices\n')
            outFID.write('  %f  %f  %f  %f\n' % (k * self.hmin, k * self.vmin,
                                                 k * self.hmax, k * self.vmax))

            outFID.write('sphere_sizes_and_positions\n')

            for i in range(len(self.spheres)):
                a = self.spheres.a[i]
                if a > 0:  # consider only positive radii
                    x = self.spheres.x[i]
                    y = self.spheres.y[i]
                    z = self.spheres.z[i]
                    n = self.spheres.materials[i].get_n(lam)
                    k = self.spheres.materials[i].get_k(lam)
                    outFID.write('  %.4f  %.4f  %.4f  %.4f  %.3f  %.3f \n' %
                                 (a, x, y, z, n, k))

            outFID.write('end_of_options\n')
            outFID.close()

            # run the binary
            if sys.platform == 'win32':
                si = subprocess.STARTUPINFO()
                si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                subprocess.call('%s scriptParams.inp > NUL' % self.command,
                                shell=True, startupinfo=si, cwd=tmpdir)
            else:
                subprocess.call('%s scriptParams.inp > /dev/null' %
                                self.command, shell=True, cwd=tmpdir)

            # parse the simulation results
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
        ''' save field data '''
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

    def plot(self):
        '''
        Show 2D field distribution
        '''
        x, y = self._get_grid_hv()
        xx, yy = np.meshgrid(x, y)
        xx = np.transpose(xx)
        yy = np.transpose(yy)
        zz = self.field
        plt.pcolormesh(xx, yy, zz, cmap='hot', shading='auto')
        plt.colorbar()
        if self.paramDict['near_field_plane_coord'] == 1:
            plt.xlabel('Y, nm')
            plt.ylabel('Z, nm')
        elif self.paramDict['near_field_plane_coord'] == 2:
            plt.xlabel('Z, nm')
            plt.ylabel('X, nm')
        elif self.paramDict['near_field_plane_coord'] == 3:
            plt.xlabel('X, nm')
            plt.ylabel('Y, nm')
        plt.gca().set_aspect('equal')
        plt.show()
        return plt


if __name__ == '__main__':
    from mstm_studio.mstm_spectrum import Material, ExplicitSpheres
    mat1 = Material(os.path.join('nk', 'etaSilver.txt'))
    nf = NearField(wavelength=340)
    nf.environment_material = 'glass'
    nf.set_plane(plane='xz', hmin=-10, hmax=20, vmin=-15, vmax=15, step=0.25)

    spheres = ExplicitSpheres(2, [0, 0, 0, 5, 0, 0, 11, 3],
                              mat_filename=2*[mat1])
    nf.set_spheres(spheres)
    nf.simulate()
    nf.plot()
    nf.write('nearfield.dat')

    print('See you!')
