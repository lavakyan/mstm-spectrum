
from mstm_studio.mstm_spectrum import SPR, Material, SpheresOverlapError
# from mstm_studio.nearfield import NearField
import os   # file path operations
import datetime
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

## auxilary functions to write inp file
def _float2str(f):
    return str(f).replace('e', 'd', 1)

def _complex2str(c):
    return f'({_float2str(c.real)},{_float2str(c.imag)})'

def _bool2str(b):
    if b:
        return '.true.'
    else:
        return '.false.'

def _value2str(v):
    if isinstance(v, str):
        return v
    elif isinstance(v, bool):
        return _bool2str(v)
    elif isinstance(v, complex):
        return _complex2str(v)
    else:  # float, int
        return _float2str(v)

def _list2str(vs):
    result = ''
    for v in vs:
        result = f'{result} {_value2str(v)},'
    return result
## end of auxilary functions

class SPR_v4(SPR):
    '''
    Class for calculation of surface plasmin resonance (SPR),
    running MSTM external code.
    The MSTM executable should be set in MSTM_BIN environment
    variable. Default is ~/bin/mstm_v4.x (not yet)
    '''

    paramDict = {
      'number_spheres': 0,
      'length_scale_factor': 1.0,          # 2π/λ[nm]
      'ref_index_scale_factor': 1.0+0.0j,  # multiplier for spheres
      'number_plane_boundaries': 0,        # layered environment
      'layer_ref_index': 1.0+0.0j,         # refractive indeces of layeres
      'layer_thickness': '',       # thiknesses of layers.
                                   # 0th layer is below 0 by Z.
                                   # others layers with this thicknesses
      #  'medium_chiral_factor': 0.0+0.0j,
      'periodic_lattice': False,        # new in ver.4
      'cell_width': [20, 20],           # periodic 2D lattice parameters

      'mie_epsilon': 1.0E-12,           # Convergence criterion for determining the number of orders
                                        # in the Mie expansions. Negative value - number of orders.
      'translation_epsilon': 1.0E-8,    # Convergence criterion for estimating the maximum order of the cluster T matrix
      'solution_epsilon': 1.0E-8,       # Precision of linear equation system solution
      'max_iterations': 5000,           # with account of all iterations
      'translation_epsilon': 1E-6,      # Error criterion for determining truncation degree when expanding fields
      't_matrix_convergence_epsilon': 1.0E-6,
      'max_t_matrix_order':  100,       # up to 120
      #  'plane_wave_epsilon': 1E-3,       # Precision of expansion of incedent field (both for palne and gaussian waves)
      #  'iterations_per_correction': 20,  # ignored for big 'near_field_translation_distance'
      'calculate_scattering_matrix': True,
      #  'near_field_translation_distance': 1.0E6,  # can be big real, small real or negative. TWEAK FOR PERFORMANCE
      'random_orientation': True,
      # 'incidence_average': True,          # alternative way to `random_orientation` (Monte-Carlo)
      # 'number_incident_directions': 100,  # ^^
      # 'azimuthal_average' : True,
      'incident_beta_deg': 0,           # parameters for fixed orientation
      'incident_alpha_deg': 0,
      # 'incident_frame': True,           # scattering matrix at θ = 0 corresponds to the incident direction, and
                                          # θ = β, ϕ = 180◦ would point in the z direction in the sphere coordinate system
      'incident_frame': False,          # scattering matrix with respect to the sphere coordinate system, so that θ = 0 would point
                                        # in the z direction, and θ = β, ϕ = α would point in the incident direction
      'scattering_map_model': 0,        # 0 - prints the scattering matrix at discrete values of θ over a circle
                                        # 1 - prints full 2D scattering matrix
      'normalize_s11': True,
      'gaussian_beam_constant': 0,      # CB = 1/(k ω0). CB = 0 - plane wave
      'gaussian_beam_focal_point': [0.0, 0.0, 0.0],  # does not alters results for plane wave and random orientations
      # 'write_sphere_data': False,     # removed from mstm2023

      'output_file': 'test.dat',            # should change for each run

      'calculate_near_field': False,   # no near field calculations
    }
    # keys that require setup at every wavelength
    local_keys = ['output_file', 'length_scale_factor',
                  'number_plane_boundaries', 'layer_ref_index',
                  't_matrix_file']

    _search_path_win =['mstm2023.exe'] + SPR._search_path_win
    _search_path_nix = ['~/bin/mstm2023.x', './mstm2023.x'] + SPR._search_path_nix

    def _write_input(self, tmpdir):
        '''
        Writes file scriptParams.inp in specified directory.
        Input:
        tmpdir -- temporaty directory to store input file
        '''
        print('Using temporary directory: %s' % tmpdir)
        outFID = open(os.path.join(tmpdir, 'scriptParams.inp'), 'w')
        outFID.write('!**********************************\n')
        outFID.write('!  MSTM input for SPR calculation\n')
        outFID.write('!  Generated by python script\n')
        outFID.write('!  %s\n' %
                     datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
        outFID.write('!**********************************\n')
        for key in self.paramDict.keys():
            if key not in self.local_keys:
                outFID.write(key + '\n')
                value = self.paramDict[key]
                if isinstance(value, list):
                    svalue = _list2str(value)
                else:
                    svalue = _value2str(value)
                outFID.write('%s \n' % svalue)

        for wl in self.wavelengths:
            k0 = 2.0 * np.pi / wl
            outFID.write('!**********************************\n')
            outFID.write('!  Wavelength  %.3f \n' % wl)
            outFID.write('!**********************************\n')
            outFID.write('output_file\n')
            outFID.write(f'mstm_l{wl*1000:.0f}.out\n')
            outFID.write('append_output_file\n')
            outFID.write('  .false.\n')
            outFID.write('length_scale_factor\n')
            outFID.write('  %.6f\n' % k0)
            outFID.write('layer_ref_index\n')
            outFID.write(f'  {_complex2str(self._environment_material.get_nk(wl))}\n')

            outFID.write('sphere_data\n')
            for i in range(len(self.spheres)):
                a = np.abs(self.spheres.a[i])  # non-negative value
                x = self.spheres.x[i]
                y = self.spheres.y[i]
                z = self.spheres.z[i]
                self.spheres.materials[i].D = 2 * a
                n = self.spheres.materials[i].get_n(wl)
                k = self.spheres.materials[i].get_k(wl)
                outFID.write('  %.4f, %.4f, %.4f, %.4f, (%.3f, %.3f) \n' %
                             (x, y, z, a, n, k))
                             # ~ (x*k0, y*k0, z*k0, a*k0, n, k))
            outFID.write('end_of_sphere_data\n')
            if wl < self.wavelengths[-1]:  # not the last
                outFID.write('new_run\n')

        outFID.write('end_of_options\n')
        outFID.close()
        return

    def _read_output(self, tmpdir):
        if self.paramDict['periodic_lattice']:
            self.reflectance = []
            self.absorptance = []
            self.transmittance = []
            self.reflectance_par = []
            self.absorptance_par = []
            self.transmittance_par = []
            self.reflectance_ort = []
            self.absorptance_ort = []
            self.transmittance_ort = []
            for wl in self.wavelengths:
                fnl = os.path.join(tmpdir, f'mstm_l{wl*1000:.0f}.out')
                with open(fnl, 'r') as fout:
                    while True:
                        line = fout.readline()
                        if not line:
                            raise Exception(f'Unexpected end of file: mstm_l{wl*1000:.0f}.out')
                        if 'scattering by periodic lattice ' in line:
                            break
                        elif 'unit cell reflectance, absorptance, transmittance (unpol, par, perp)' in line:
                            values = map(float,
                                         fout.readline().strip().split())
                            values = list(values)
                            self.reflectance.append(float(values[0]))
                            self.absorptance.append(float(values[1]))
                            self.transmittance.append(float(values[2]))
                            self.reflectance_par.append(float(values[3]))
                            self.absorptance_par.append(float(values[4]))
                            self.transmittance_par.append(float(values[5]))
                            self.reflectance_ort.append(float(values[6]))
                            self.absorptance_ort.append(float(values[7]))
                            self.transmittance_ort.append(float(values[8]))
                os.remove(fnl)
            self.reflectance = np.array(self.reflectance)
            self.absorptance = np.array(self.absorptance)
            self.transmittance = np.array(self.transmittance)
            self.reflectance_par = np.array(self.reflectance_par)
            self.absorptance_par = np.array(self.absorptance_par)
            self.transmittance_par = np.array(self.transmittance_par)
            self.reflectance_ort = np.array(self.reflectance_ort)
            self.absorptance_ort = np.array(self.absorptance_ort)
            self.transmittance_ort = np.array(self.transmittance_ort)
            # ext c.s. = -lnT / (n * L) = -lnT * cell_x*cell_y / N
            # n = N / (cell_x*cell_y*L) - concentration
            # ext. eff = ext c.s. / (sph.area)
            # sph.area is hard to calculate.
            # For now - just a sum of c.s. of all spheres
            cell_x, cell_y = self.paramDict['cell_width']
            coef = cell_x * cell_y / (np.pi * np.sum(self.spheres.a**2) * len(self.spheres))
            self.extinction = -coef * np.log(self.transmittance)
            self.extinction_par = -coef * np.log(self.transmittance_par)
            self.extinction_ort = -coef * np.log(self.transmittance_ort)
            return self.wavelengths, self.extinction
        elif self.paramDict['random_orientation']:  # random orient.
            self.extinction = []
            self.absorbtion = []
            self.scattering = []
            for wl in self.wavelengths:
                fnl = os.path.join(tmpdir, f'mstm_l{wl*1000:.0f}.out')
                with open(fnl, 'r') as fout:
                    while True:
                        line = fout.readline()
                        if not line:
                            raise Exception(f'Unexpected end of file: mstm_l{wl*1000:.0f}.out')
                        if 'total scattering' in line: # next line after
                            break                      # required
                        elif 'total extinction, absorption, scattering efficiencies' in line:
                            values = map(float,
                                         fout.readline().strip().split())
                            values = list(values)
                            self.extinction.append(float(values[0]))
                            self.absorbtion.append(float(values[1]))
                            self.scattering.append(float(values[2]))
                os.remove(fnl)
            self.extinction = np.array(self.extinction)
            self.absorbtion = np.array(self.absorbtion)
            self.scattering = np.array(self.scattering)
            return self.wavelengths, self.extinction
        else:  # fixed orientation, no periodicity
            self.extinction = []  # unploraized
            self.absorbtion = []
            self.scattering = []
            self.extinction_par = []  # parallel polarization (\hat \alpha)
            self.absorbtion_par = []
            self.scattering_par = []
            self.extinction_ort = []  # perpendicular polarization (\hat \beta)
            self.absorbtion_ort = []
            self.scattering_ort = []
            for wl in self.wavelengths:
                fnl = os.path.join(tmpdir, f'mstm_l{wl*1000:.0f}.out')
                with open(fnl, 'r') as fout:
                    while True:
                        line = fout.readline()
                        if not line:
                            raise Exception(f'Unexpected end of file: mstm_l{wl*1000:.0f}.out')
                        if 'down and up hemispherical scattering efficiencies' in line:
                            break
                        elif 'total extinction, absorption, scattering efficiencies' in line:
                            # total extinction, absorption, scattering efficiencies (unpol, par, perp incidence)
                            values = map(float,
                                         fout.readline().strip().split())
                            values = list(values)
                            self.extinction.append(float(values[0]))
                            self.absorbtion.append(float(values[1]))
                            self.scattering.append(float(values[2]))
                            self.extinction_par.append(float(values[3]))
                            self.absorbtion_par.append(float(values[4]))
                            self.scattering_par.append(float(values[5]))
                            self.extinction_ort.append(float(values[6]))
                            self.absorbtion_ort.append(float(values[7]))
                            self.scattering_ort.append(float(values[8]))
                os.remove(fnl)
            self.extinction = np.array(self.extinction)
            self.absorbtion = np.array(self.absorbtion)
            self.scattering = np.array(self.scattering)
            self.extinction_par = np.array(self.extinction_par)
            self.absorbtion_par = np.array(self.absorbtion_par)
            self.scattering_par = np.array(self.scattering_par)
            self.extinction_ort = np.array(self.extinction_ort)
            self.absorbtion_ort = np.array(self.absorbtion_ort)
            self.scattering_ort = np.array(self.scattering_ort)
            return self.wavelengths, self.extinction

    def plot(self):
        '''
        Plot results with matplotlib.pyplot
        '''
        if self.paramDict['periodic_lattice']:  # PBC
            if self.paramDict['random_orientation']:  # random
                plt.plot(self.wavelengths, self.transmittance, 'r-', label='T')
                plt.plot(self.wavelengths, self.absorptance, 'g-', label='A')
                plt.plot(self.wavelengths, self.reflectance, 'b-', label='R')
            else:
                plt.plot(self.wavelengths, self.transmittance_par, 'r-', label='T par')
                plt.plot(self.wavelengths, self.absorptance_par, 'g-', label='A par')
                plt.plot(self.wavelengths, self.reflectance_par, 'b-', label='R par')
                plt.plot(self.wavelengths, self.transmittance_ort, 'r--', label='T ort')
                plt.plot(self.wavelengths, self.absorptance_ort, 'g--', label='A ort')
                plt.plot(self.wavelengths, self.reflectance_ort, 'b--', label='R ort')
        elif self.paramDict['random_orientation']:  # random, no PBC
            plt.plot(self.wavelengths, self.extinction, 'r-', label='extinction')
            plt.plot(self.wavelengths, self.absorbtion, 'g-', label='absorbtion')
            plt.plot(self.wavelengths, self.scattering, 'b-', label='scattering')
        else:
            plt.plot(self.wavelengths, self.extinction_par, 'r-',  label='extinction par.')
            plt.plot(self.wavelengths, self.extinction_ort, 'r--', label='extinction ort.')
            plt.plot(self.wavelengths, self.absorbtion_par, 'g-',  label='absorbtion par.')
            plt.plot(self.wavelengths, self.absorbtion_ort, 'g--', label='absorbtion ort.')
            plt.plot(self.wavelengths, self.scattering_par, 'b-',  label='scattering par.')
            plt.plot(self.wavelengths, self.scattering_ort, 'b--', label='scattering ort.')
        plt.legend()
        plt.show()
        return plt

    def write(self, filename):
        '''
        Save results to file
        '''
        if self.paramDict['random_orientation']:  # random
            fout = open(filename, 'w')
            fout.write('#Wavel.\tExtinct.\n')
            for i in range(len(self.wavelengths)):
                fout.write('%.4f\t%.8f\r\n' % (self.wavelengths[i],
                                               self.extinction[i]))
            fout.close()
        else:   # fixed
            fout = open(filename, 'w')
            fout.write('#Wavel.\tExt_par\tExt_ort\n')
            for i in range(len(self.wavelengths)):
                fout.write('%.4f\t%.8f\t%.8f\r\n' % (self.wavelengths[i],
                           self.extinction_par[i], self.extinction_ort[i]))
            fout.close()

    def set_incident_field(self, fixed=False, beta_angle=0.0,
                           alpha_angle=0.0, polarization_angle=0.0):
        '''
        Set incident wave orientation and polarization

        Parameters:

            fixed: bool
                True  - fixed orientation and polarized light
                False - average over all orientations and polarizations

            beta_angle:  float (degrees)
                polar angle (from axis Z)

            alpha_angle: float (degrees)
                azimutal angle (from axis X)

            polarization_angle: float (degrees)
                ?? not used in MSTM v.4 ??
        '''
        self.paramDict['random_orientation'] = not fixed
        if fixed:
            self.paramDict['incident_beta_deg'] = beta_angle
            self.paramDict['incident_alpha_deg'] = alpha_angle
            # self.paramDict['polarization_angle_deg'] = polarization_angle

    def set_boundary(self, pbc=False, cell_x=100., cell_y=100.):
        '''
        Set peroidic boundary conditions (PBC).
        Periodicity can be only in XY plane.
        Spheres should not intercept the PBC box.

        Parameters:

            pbc: bool
                Use PBC

            cell_x: float
                the X size of PBC box

            cell_y: float
                the Y size of PBC box
        '''
        self.paramDict['periodic_lattice'] = pbc
        self.paramDict['cell_width'] = [cell_x, cell_y]
        # TODO: sanity checks?
        if self.paramDict['random_orientation']:
            print('Switching to fixed orientation')
            self.set_incident_field(True)

    def set_layers(self, mats=[], depths=[]):
        '''
        Layers in Z direction:
        z < 0 -- governed by `environment_material`
        0 < z < depth[0] -- mats[0]
        depth[0] < z < depth[1] -- mats[1]
        etc.

        Defaul is no layers.

        mats: list of Materials
            materials of layers

        depths: list of float
            the size of layers
        '''

        # TODO
        pass


class NearField_v4(SPR_v4):
    '''
    Calculate field distribution map at fixed wavelength
    using MSTM v.4 code (significantly reworked).

    Note: The precision of internal field is sensitive
    to the `near_field_expansion_spacing` and
    `near_field_expansion_order` parameters. Please tweak
    them if obtained too high values.

    Example of usage:

    .. code-block:: python

        from mstm_studio.mstm_spectrum import Material, ExplicitSpheres
        wl = 240  # wavelength
        # grid:
        hmin, hmax, vmin, vmax, step = -20, 20, -20, 20, 0.25
        a = 10  # sphere radius
        nf = NearField_v4(wavelength=wl,
                          incident_default_mode=True)
        nf.environment_material = 1.5  # matrix n
        spheres = ExplicitSpheres(2, [-6, 0, 0, 5, 6, 0, 0, 5],
                                  mat_filename=2*[Material(0.2 + 0.8j)])
        nf.set_spheres(spheres)
        nf.simulate()  # do calc

        # plot results
        fig, ax = plt.subplots(1, 1, figsize=(5, 6))
        nf.plot(fig=fig, axs=ax, mode='par')
        plt.tight_layout()
        plt.savefig('nf_mstm4.png')
        plt.show()

    '''
    def __init__(self, wavelength, mstm_path='~/bin/mstmswd.x',
                 environment_material='Air', temp_dir=None,
                 incident_default_mode=True):
        # ~ print(mstm_path)
        super().__init__([wavelength], mstm_path,
                         environment_material, temp_dir)
        # ~ print(self.command)
        self.paramDict['calculate_near_field'] = True  # do nearfield
        self.set_incident_field(fixed=True,
                                beta_angle=0.0,
                                alpha_angle=0.0)
        if incident_default_mode:
            # controls where the incident field appears in the calculation results:
            # = 1 - is the standard model,
            # where the external field is = scattered + incident and the field
            # inside the spheres is calculated from the internal field expansions;
            # != 1 has the external field due solely to the scattered field,
            # and the field inside the particles is now internal-incident.
            self.paramDict['near_field_calculation_model'] = 1
        else:
            self.paramDict['near_field_calculation_model'] = 0
        self.paramDict['store_surface_vector'] = True  # unless doubt
        self.paramDict['near_field_expansion_spacing'] = 1  # default 5
        self.paramDict['near_field_expansion_order'] = 10  # def 10. higher - more accurate, but slower
        self.paramDict['near_field_output_file'] = 'nf-temp.dat'
        self.set_plane()

    def set_plane(self, plane='XY', hmin=-10., hmax=10.,
                  vmin=-10., vmax=10., step=1., offset=0.):
        '''
        Determine the plane and grid for near field computation.

        plane: one of 'yz'|'zx'|'xy'
        hmin, hmax, vmin, vmax: horizontal and vertical sizes
        step:   size of the grid grain
        offset: shift of the plane
        '''
        hmin += -step / 100.  # add small value in
        vmin += -step / 100.  # attempt to diminish
        hmax += step / 100.   # rounding problems
        vmax += step / 100.
        self.hmin = hmin
        self.hmax = hmax
        self.vmin = vmin
        self.vmax = vmax
        self.step = step
        self.offset = offset
        self.nh = int(np.round((hmax - hmin) / step))  # MSTM has its
        self.nv = int(np.round((vmax - vmin) / step))  # own roundings..
        if plane.upper() in ['YX', 'XY']:
            self.plane = 'XY'
            min_border=[hmin, vmin, offset]
            max_border=[hmax, vmax, offset]
        elif plane.upper() in ['XZ', 'ZX']:
            self.plane = 'XZ'
            min_border=[hmin, offset, vmin]
            max_border=[hmax, offset, vmax]
        elif plane.upper() in ['YZ', 'ZY']:
            self.plane = 'YZ'
            min_border=[offset, hmin, vmin]
            max_border=[offset, hmax, vmax]
        else:
            raise Exception(f'Wrong plane {plane}')

        k0 = 2 * np.pi / self.wavelengths[0]
        self.paramDict['near_field_minimum_border'] = [k0 * v for v in min_border]
        self.paramDict['near_field_maximum_border'] = [k0 * v for v in max_border]
        self.paramDict['near_field_step_size'] = k0 * step
        print('Estimated field grid: %ix%i' % (self.nh, self.nv))
        return

    def _read_output(self, tmpdir):
        '''
        read nearfield spatial distribution
        from file specified in `near_field_output_file`

        Stored internal values of fields E, H, both
        could be par(∥) or ort (⊥) towards incidence,
        all projected on x, y, z directions
        with complex values. 24 items in total

        Returns:
            2d array of |E|^2 - the square of total electric field

            The complex values of field compontents are stored
            as internal variables:

                self.Epar_x, self.Epar_y, self.Epar_z,
                self.Hpar_x, self.Hpar_y, self.Hpar_z,
                self.Eort_x, self.Eort_y, self.Eort_z,
                self.Hort_x, self.Hort_y, self.Hort_z

        '''
        fn = os.path.join(tmpdir,
                          self.paramDict['near_field_output_file'])
        with open(fn) as fout:
            fout.readline()  # skip " run number:"
            fout.readline()  # skip value
            nsph = int(fout.readline().strip())  # no. of spheres in plane
            for _ in range(nsph):
                fout.readline()  # skip
            nbou = int(fout.readline().strip())  # no. of layer boundaries
            for _ in range(nbou):
                fout.readline()  # skip
            fout.readline()  # skip echo of near_field_minimum_border
            fout.readline()  # skip echo of near_field_maximum_border
            s = fout.readline().strip()
            nx, ny, nz = [int(v) for v in s.split()]
            print('dimensions from mstm')
            print(nx, ny, nz)
        nskip = 2 + 1 + nsph + 1 + nbou + 2 + 1
        data = np.loadtxt(fn, skiprows=nskip)
        # The next N lines have 27 columns associated with each calculation
        # point: x, y, z, E∥ , H∥ , E⊥ , H⊥ ;
        # each vector field has 6 columns: Re Ex , Im Ex , Re Ey , and so on,
        # and ∥, ⊥ correspond to the parallel and perpendicular incident
        # polarization states.
        # code in mstm-scatprops-26.f90:
        # write(outputunit,'(27es12.4)') rpos(:),earray(:,1,ix,iy),harray(:,1,ix,iy), &
        #                earray(:,2,ix,iy),harray(:,2,ix,iy)
        # confirms my undertanding. But - numberical comparison with v3 is bad.
        print(data.shape)
        if self.plane == 'XY':
            self.nh = nx
            self.nv = ny
        elif self.plane == 'XZ':
            self.nh = nx
            self.nv = nz
        elif self.plane == 'YZ':
            self.nh = ny
            self.nv = nz

        self.Epar_x = np.reshape(data[:, 3], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:, 4], [self.nv, self.nh])
        self.Epar_y = np.reshape(data[:, 5], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:, 6], [self.nv, self.nh])
        self.Epar_z = np.reshape(data[:, 7], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:, 8], [self.nv, self.nh])
        self.Hpar_x = np.reshape(data[:, 9], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,10], [self.nv, self.nh])
        self.Hpar_y = np.reshape(data[:,11], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,12], [self.nv, self.nh])
        self.Hpar_z = np.reshape(data[:,13], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,14], [self.nv, self.nh])
        self.Eort_x = np.reshape(data[:,15], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,16], [self.nv, self.nh])
        self.Eort_y = np.reshape(data[:,17], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,18], [self.nv, self.nh])
        self.Eort_z = np.reshape(data[:,19], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,20], [self.nv, self.nh])
        self.Hort_x = np.reshape(data[:,21], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,22], [self.nv, self.nh])
        self.Hort_y = np.reshape(data[:,23], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,24], [self.nv, self.nh])
        self.Hort_z = np.reshape(data[:,25], [self.nv, self.nh]) + \
                 1j * np.reshape(data[:,26], [self.nv, self.nh])

        self.Epar_xyz = np.array(
                    [self.Epar_x,
                     self.Epar_y,
                     self.Epar_z])
        self.Eort_xyz = np.array(
                    [self.Eort_x,
                     self.Eort_y,
                     self.Eort_z])
        self.field = np.sum(np.abs(self.Epar_xyz)**2, axis=0)
        return self.field

    def _get_grid_hv(self):
        h = np.linspace(self.hmin, self.hmax, self.nh)
        v = np.linspace(self.vmin, self.vmax, self.nv)
        return h, v

    def write(self, filename):
        ''' save field data to text file'''
        h, v = self._get_grid_hv()
        with open(filename, 'w') as fout:
            fout.write(f'# {self.plane[0]}[nm]\t{self.plane[1]}[nm]\t|E|^2\r\n')
            for i, x in enumerate(h):
                for j, y in enumerate(v):
                    fout.write('%.4f\t%.4f\t%.8f\r\n' % (x, y,
                                                         self.field[j, i]))

    def plot(self, fig=None, axs=None, caxs=None, mode='par'):
        '''
        Show 2D field distribution

        Parameters:

        fig:
            matplotlib figure
        axs:
            matplotlib axes
        caxs:
            matplotlib axes for colorbar
        mode: 'total' | 'par' | 'ort'
           plot averaged field or parallel or orthogonal
           components

        Returns:
            filled/created fig and axs objects
        '''
        x, y = self._get_grid_hv()
        xx, yy = np.meshgrid(x, y)
        zz = self.field
        if mode == 'par':
            zz = np.abs(self.Epar_x)**2 + \
                 np.abs(self.Epar_y)**2 + \
                 np.abs(self.Epar_z)**2
        elif mode == 'ort':
            zz = np.abs(self.Eort_x)**2 + \
                 np.abs(self.Eort_y)**2 + \
                 np.abs(self.Eort_z)**2
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
        axs.set_xlabel(f'{self.plane[0]}, nm')
        axs.set_ylabel(f'{self.plane[1]}, nm')
        axs.set_aspect('equal')
        if flag:
            plt.show()
        return fig, axs


if __name__ == '__main__':
    from mstm_studio.mstm_spectrum import Material, ExplicitSpheres

    if True:
        mat1 = Material(os.path.join('nk', 'etaGold.txt'))
        mat2 = Material(os.path.join('nk', 'etaSilver.txt'))
        wls = np.linspace(300, 800, 100)
        print('old SPR')
        spr = SPR(wls)
        spr.environment_material = 'air'
        spheres = ExplicitSpheres(2, [-20, 0, 0, 10, 10, 0, 0, 12],
                                  mat_filename=[mat1, mat2])
        spr.set_spheres(spheres)
        # ~ spr.set_incident_field(fixed=False)
        # ~ spr.set_incident_field(fixed=True, azimuth_angle=90.0, polar_angle=90.0,
                               # ~ polarization_angle=45.0)
        spr.simulate()
        # input()
        # spr.write('test.dat')
        spr.plot()

        print('new SPR')
        spr = SPR_v4(wls)  #, temp_dir='./temp/')
        spr.environment_material = 'air'
        spr.set_spheres(spheres)
        spr.set_incident_field(fixed=False)
        # ~ spr.set_incident_field(fixed=True, beta_angle=90.0, alpha_angle=90.0)
        spr.simulate()
        spr.write('test.dat')
        spr.plot()

    if True:
        wl = 240
        matsph = 0.5 + 0.1j
        matrix = 1.5
        hmin, hmax, vmin, vmax, step = -20, 20, -20, 20, 0.25
        a = 10

        nf = NearField_v4(wavelength=wl,
                          incident_default_mode=True)
        nf.environment_material = matrix
        spheres = ExplicitSpheres(1, [0, 0, 0, a],
                                  mat_filename=Material(matsph))
        nf.set_plane(plane='xz', hmin=hmin, hmax=hmax,
                     vmin=vmin, vmax=vmax, step=step)
        # ~ spheres = ExplicitSpheres(2, [-6, 0, 0, 5, 6, 0, 0, 5],
                                  # ~ mat_filename=2*[Material('nk/etaGold.txt')])
        # ~ nf.set_plane(plane='YZ', hmin=-20., hmax=20.,
                      # ~ vmin=-15., vmax=15., step=0.5, offset=0.)
        nf.set_spheres(spheres)
        nf.simulate()

        fig, ax = plt.subplots(1, 1, figsize=(5, 6))
        nf.plot(fig=fig, axs=ax, mode='par')
        plt.tight_layout()
        plt.savefig('nf_mstm4.png')
        plt.show()
        # ~ nf.plot(mode='ort')
        # ~ nf.write('nearfield.dat')

    print('See you!')
