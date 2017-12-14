#! /usr/bin/env python
#
# Support module generated by PAGE version 4.9
# In conjunction with Tcl version 8.6
#    Dec 03, 2017 12:51:25 AM
#    Dec 03, 2017 09:18:10 PM
#    Dec 03, 2017 11:35:42 PM
# manually edited afterwards

import sys
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#from matplotlib.backend_bases import key_press_handler
from itertools import cycle
import numpy as np
from scipy import interpolate
from mstm_spectrum import (Material, SingleSphere, Background,
                           LinearBackground, LorentzBackground, SPR,
                           LogNormalSpheres)
from fit_spheres_optic import (Fitter, FixConstraint, EqualityConstraint,
                               ConcentricConstraint)
#import threading
#import time
import copy
try:
    from Tkinter import Frame, Label, Entry, Toplevel, Spinbox, StringVar
    from tkColorChooser import askcolor
except ImportError:
    from tkinter import Frame, Label, Entry, Toplevel, Spinbox, StringVar
    from tkinter.colorchooser import askcolor
import tkFileDialog, tkSimpleDialog, tkMessageBox
try:
    import ttk
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    py3 = 1



def btConstraintsClick(event=None):
    global w, spheres
    if spheres is None:
        nspheres = 0
    else:
        nspheres = len(spheres)
    w.constr_win_app.show_window(nspheres)


fitter = None

def btStartFitClick(event=None):
    global w, spheres, fitter

    #~ if spheres is None:
        #~ tkMessageBox.showwarning('Warnign', 'No spheres to fit')
    if (fitter is not None) and fitter.isAlive():
        tkMessageBox.showwarning('Warning', 'Fitting is already running')
        return
    if (fitter is not None):
        fitter = create_fitter(get_wavelengths(), w.edExpFileName.get(),
                               w.cbBkgMethod.get())
    fitter.set_scale(float(w.edSpecScale.get()))
    fitter.set_background(w.cbBkgMethod.get(),
                          initial_values=get_bkg_params())
    fitter.set_spheres(copy.deepcopy(spheres))
    if spheres is not None:  # fit only background
        fitter.set_matrix(get_matrix_material())
    else:
        print('WARNING: setting zero material!')
        fitter.set_matrix(Material(0))
    # set constraints
    fitter.add_constraint(w.constr_win_app.get_constraints_list())
    s = fitter.report_freedom()
    print(s)
    s += fitter.report_result('Initial parameters')
    print(s)
    if tkMessageBox.askokcancel('Continue?', s):
        fitter.start()

def btStopFitClick(event=None):
    global fitter
    fitter.stop()
    #~ fitting_thread._Thread_stop()

def btLoadExpClick(event=None):
    global w, fitter
    ftypes = [('Text files', '*.txt'), ('Dat files', '*.dat'),
              ('Exp. files', '*.exp'), ('All files', '*')]
    fn = tkFileDialog.askopenfilename(filetypes=ftypes)
    w.edExpFileName.delete(0, 'end')
    w.edExpFileName.insert(0, fn)
    if fn:
        wls = get_wavelengths()
        try:
            fitter = create_fitter(wls, fn, w.cbBkgMethod.get())
        except Exception as err:
            tkMessageBox.showerror('Error', str(err))
            return
        btPlotExpClick(event)

def btPlotExpClick(event=None):
    global fitter
    axs = w.TPanedwindow3_p1.axs
    axs.clear()
    axs.plot(fitter.wls, fitter.exp, 'ro', label='Exp.')
    x, y = load_spec(w.model_fn)
    axs.plot(x, y, 'b-', label='Model')
    axs.set_ylabel('Intensity')
    axs.set_xlabel('Wavelength, nm')
    axs.legend()
    w.TPanedwindow3_p1.canvas.draw()

def btCalcSpecClick(event=None):
    global w, materials, spheres
    wls = get_wavelengths()
    # create SPR object
    spr = SPR(wls)
    spr.environment_material = get_matrix_material()
    spr.set_spheres(spheres)
    # calculate!
    spr.simulate(w.model_fn)
    # tkMessageBox.showinfo('MSTM studio', 'Calculation finished')
    btPlotSpecClick(event)

def btSaveSpecClick(event=None):
    global w, root
    x, y = load_spec(w.model_fn)
    global root
    ftypes = [('Text files', '*.txt'), ('Dat files', '*.dat'), ('All files', '*')]
    fn = tkFileDialog.asksaveasfilename(filetypes=ftypes)
    if fn:
        try:
            f = open(fn, 'w')
            f.write('#Lambda(nm)\tExctinction\r\n')
            for i in xrange(len(x)):
                f.write(' %.3f\t%.6f\r\n' % (x[i], y[i]))
            print('Saved to %s' % fn)
        finally:
            f.close()

def btPlotSpecClick(event=None):
    global w
    x, y = load_spec(w.model_fn)
    axs = w.TPanedwindow3_p1.axs
    axs.clear()
    axs.plot(x, y, 'b-', label='Model')
    axs.set_ylabel('Intensity')
    axs.set_xlabel('Wavelength, nm')
    axs.legend()
    w.TPanedwindow3_p1.canvas.draw()

def load_spec(filename):
    global spheres
    try:
        model = np.genfromtxt(filename)
    except Exception as err:
        tkMessageBox.showerror('Error', 'Error loading file "%s"\n%s'
                               % (filename, str(err)))
        return
    wls = get_wavelengths()
    if (spheres is None) or (len(spheres)==0):
        print('* No spheres case *')
        x = wls
        y = np.zeros_like(x)
    else:
        x = model[:, 0]
        y = model[:, 1]
        if (len(wls) != len(x)) or (np.abs(wls[0]-x[0])>1E-3) or (np.abs(wls[1]-x[1])>1E-3):
            result = tkMessageBox.askquestion('Grid mismatch',
              'Requested and stored wavelengths are different.\nProbably you will need to Calcualte first.\nProceed and interpolate to requeste scale?')
            if result == 'yes':
                f = interpolate.interp1d(x, y)
                x = wls
                y = f(wls)
            else:
                return
        scale = get_scale()
        y = scale * y
    cbBkgMethodSelect()
    y += background.get_bkg(get_bkg_params())
    return x, y


background = None

def cbBkgMethodSelect(event=None):
    global w, background
    s = w.cbBkgMethod.get()  # ['Constant', 'Linear', 'Lorentz']
    print(s)
    wls = get_wavelengths()
    if s == 'Linear':
        background = LinearBackground(wls)
    elif s == 'Lorentz':
        background = LorentzBackground(wls)
    else:  # 'Constant':
        background = Background(wls)

def btPlotBkgClick(event=None):
    global w, background
    #~ if background is None:
    cbBkgMethodSelect(event)
    w.TPanedwindow3_p1.axs.clear()
    params = get_bkg_params()
    background.plot(params, fig=w.TPanedwindow3_p1.fig, axs=w.TPanedwindow3_p1.axs)
    w.TPanedwindow3_p1.canvas.draw()

def get_bkg_params():
    global w, background
    result = []
    n = background.number_of_params()
    try:
        if n > 0:
            result.append(float(w.edBkg1.get()))
        if n > 1:
            result.append(float(w.edBkg2.get()))
        if n > 2:
            result.append(float(w.edBkg3.get()))
    except ValueError as err:
        tkMessageBox.showerror('Error', 'Bad value.\n%s' % err)
    return result


spheres = None

def btAddSphClick(master=None):
    global w, root
    global spheres
    dial = SphereDialog(root) # 10, 0.0, 0.0, 0.0, 'm0')
    if dial.result is None:
        return
    a, x, y, z, key = dial.result
    try:
        sphere = SingleSphere(a=a, x=x, y=y, z=z, mat_filename=materials[key][0])
    except Exception as err:
        tkMessageBox.showerror('Error', err)
        return
    if spheres is None:
        spheres = sphere
    else:
        spheres.append(sphere)
    print('len(spheres) = %i' % len(spheres))
    update_spheres_tree()
    btPlotSphClick(master)

def btGenerateSpheresClick(event=None):
    global materials, spheres
    dial = GenerateSpheresDialog(root, 8, 10.0, 5.0, materials.keys())
    if dial.result is None:
        return
    N, a, d, key = dial.result
    try:
        sphere = LogNormalSpheres(N=N, mu=a, sigma=1E-3, d=d, mat_filename=materials[key][0])
    except Exception as err:
        tkMessageBox.showerror('Error', err)
        return
    if spheres is None:
        spheres = sphere
    else:
        spheres.extend(sphere)
    print('len(spheres) = %i' % len(spheres))
    update_spheres_tree()
    btPlotSphClick(event)

def btImportSpheres(master=None):
    global root, spheres, materials
    # open file dialog
    ftypes = [('Text files', '*.txt'), ('Data files', '*.dat'),
             ('Input files', '*.inp'), ('All files', '*')]
    #~ dlg = tkFileDialog.Open(root, filetypes=ftypes)
    #~ fn = dlg.show()
    fn = tkFileDialog.askopenfilename(filetypes=ftypes)
    if fn:
        try:
            data = np.genfromtxt(fn)
        except Exception as err:
            tkMessageBox.showerror('Error', 'Can not read %s\n%s' % (fn, err))
        print(data.shape[1])
        if data.shape[1] == 4:  # a,x,y,z,n,k
            if len(materials) == 0:
                #materials['m0'] = Material('4+2j')
                add_material('m0', Material('4+2j'))
                update_materials_tree()
            mat_key = next(iter(materials))  # set same material for all spheres
            for row in data:
                print(row)
                sphere = SingleSphere(a=row[0], x=row[1], y=row[2], z=row[3],
                                      mat_filename=materials[mat_key][0])
                if spheres is None:
                    spheres = sphere
                else:
                    spheres.append(sphere)
            update_spheres_tree()
            update_spheres_canvas()
        elif data.shape[1] == 6:  # a,x,y,z,n,k
            for row in data:
                print(row)
                mat = Material('%.3f%+.3fj' % (row[4], row[5]))
                mat_key = find_mat_key(mat)
                if mat_key is None:
                    mat_key = gen_mat_key()
                    #materials[mat_key] = mat
                    add_material(mat_key, mat)
                    print('Created material: %s' % materials[mat_key][0])
                sphere = SingleSphere(a=row[0], x=row[1], y=row[2], z=row[3],
                                      mat_filename=materials[mat_key][0])
                if spheres is None:
                    spheres = sphere
                else:
                    spheres.append(sphere)
            update_materials_tree()
            update_spheres_tree()
            update_spheres_canvas()
        else:
            tkMessageBox.showerror('Imort failed', 'Expected 4- or 6- columns in file.\n%s' % fn)

def btExportSpheres(event=None):
    global root, spheres, materials
    # open file dialog
    ftypes = [('Text files', '*.txt'), ('Data files', '*.dat'),
             ('Input files', '*.inp'), ('All files', '*')]
    #~ dlg = tkFileDialog.Open(root, filetypes=ftypes)
    #~ fn = dlg.show()
    fn = tkFileDialog.asksaveasfilename(filetypes=ftypes)
    if fn:
        try:
            f = open(fn, 'w')
            f.write('#radius\tx\ty\tz\tn\tk\r\n')
            for i in xrange(len(spheres)):
                wl = float(w.edLambdaMin.get())
                a = spheres.a[i]
                x = spheres.x[i]
                y = spheres.y[i]
                z = spheres.z[i]
                n = spheres.materials[i].get_n(wl)
                k = spheres.materials[i].get_n(wl)
                f.write('%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\r\n' % (a, x, y, z, n, k))
        finally:
            f.close()


def btEditSphClick(master=None):
    global w, top_level, root
    global spheres

    tree = w.stvSpheres
    sel = tree.selection()
    if len(sel)>0:
        key = tree.item(sel[0], 'text')
        print(key)
        i = int(key[1:])  # remove leading 's' symbol from key to get index
        dial = SphereDialog(root, spheres.a[i], spheres.x[i], spheres.y[i], spheres.z[i],
                            find_mat_key(spheres.materials[i]))
        if dial.result is not None:
            a, x, y, z, mat = dial.result
            spheres.a[i] = a
            spheres.x[i] = x
            spheres.y[i] = y
            spheres.z[i] = z
            spheres.materials[i] = materials[mat][0]
            update_spheres_tree()
            update_spheres_canvas()

def btDelSphClick(master=None):
    global w
    tree = w.stvSpheres
    sel = tree.selection()
    if len(sel)>0:
        key = tree.item(sel[0], 'text')
        i = int(key[1:])
        spheres.delete(i)
        update_spheres_tree()
        btPlotSphClick(master)

def btPlotSphClick(master=None):
    global w, spheres
    w.canvas.camera.viewpoint = spheres.get_center('mass')
    w.canvas.camera.axes = np.identity(3)
    update_spheres_canvas()

def sync_spheres_materials():
    global w
    if (spheres is None) or (materials is None):
        return
    tree = w.stvSpheres
    for child in tree.get_children():
        sphid = tree.item(child, 'text')
        matkey = tree.item(child, 'values')[4]
        i = int(sphid[1:])
        spheres.materials[i] = materials[matkey][0]

def update_spheres_tree():
    global w
    if spheres is None:
        return
    tree = w.stvSpheres
    tree.delete(*tree.get_children())
    for i in xrange(len(spheres)):
        matkey = find_mat_key(spheres.materials[i])
        tree.insert('' , 0, text='s%i'%i, values=(spheres.a[i], spheres.x[i],
                    spheres.y[i], spheres.z[i], matkey))

def update_spheres_canvas():
    global w
    if spheres is None:
        return
    cv = w.canvas
    W = cv.winfo_width()
    H = cv.winfo_height()
    cv.delete('all')
    positions = np.stack((spheres.x, spheres.y, spheres.z), axis=-1)
    projected = cv.camera.project(positions)
    indices = projected[:, 2].argsort()  # sorted by Z-buffer
    for i in indices:
        a = spheres.a[i] * cv.camera.scale
        x = W/2 + projected[i, 0]
        y = H/2 - projected[i, 1]
        key = find_mat_key(spheres.materials[i])
        col = materials[key][1]  # was '#5544FF'
        cv.create_oval(x-a, y-a, x+a, y+a, outline='#004500', width=3, fill=col, stipple='gray50')

def mouse_wheel(event):
    global w
    if event.num == 4 or event.delta < 0:    # Lin or Win
        w.canvas.camera.zoom_in()
    elif event.num == 5 or event.delta > 0:
        w.canvas.camera.zoom_out()
    w.lbZoom['text'] = 'x%.2f' % w.canvas.camera.scale
    update_spheres_canvas()

def mouse_down(event):
    """ the idea and some code borrowed from ASE <https://wiki.fysik.dtu.dk/ase> """
    global w
    cam = w.canvas.camera
    cam.xy = (event.x, event.y)
    #~ cam.t0 = event.time
    cam.axes0 = cam.axes

def mouse_move(event):
    """ the idea and some code borrowed from ASE <https://wiki.fysik.dtu.dk/ase> """
    global w
    cam = w.canvas.camera
    x = event.x
    y = event.y
    x0, y0 = cam.xy
    a = x - x0
    b = y0 - y
    t = np.sqrt(a * a + b * b)
    if t > 0:
        a /= t
        b /= t
    else:
        a = 1.0
        b = 0.0
    c =  np.cos(0.02 * t)  # was 0.01
    s = -np.sin(0.02 * t)  # was 0.01
    rotation = np.array([(c * a * a + b * b, (c - 1) * b * a,    s * a),
                         ((c - 1) * a * b,    c * b * b + a * a, s * b),
                         (-s * a,            -s * b,             c)])
    cam.axes = np.dot(cam.axes0, rotation)
    update_spheres_canvas()

def mouse_up(event):
    #~ print('mouse up?')
    pass


materials = {}

def btDelMatClick(master=None):
    global w
    tree = w.stvMaterial
    sel = tree.selection()
    if sel:
        key = tree.item(sel[0], 'text')
        materials.pop(key)
        update_materials_tree()

def btAddMatClick(master=None):
    global w, material
    mat_str = tkSimpleDialog.askstring('Material data', 'Enter material name or refraction index')
    if mat_str:
        try:  # try create material
            mat = Material(mat_str)
        except Exception as err:
            tkMessageBox.showerror('Error', err)
            return
        key = gen_mat_key(True)
        add_material(key, mat)
        update_materials_tree()

def btLoadMatClick(master=None):
    ftypes = [('Text files', '*.txt'), ('All files', '*')]
    fn = tkFileDialog.askopenfilename(filetypes=ftypes)
    if fn:
        try:
            mat = Material(file_name=fn)  # encode() ?
        except Exception as err:
            tkMessageBox.showerror('Error', str(err))
            return
        key = gen_mat_key(True)
        add_material(key, mat)
        update_materials_tree()

def btPlotMatClick(master=None):
    global w, top_level, root
    tree = w.stvMaterial
    sel = tree.selection()
    if sel:
        key = tree.item(sel[0], 'text')
        mat = materials[key][0]
        w.TPanedwindow3_p1.axs.clear()
        mat.plot(wls=get_wavelengths(), fig=w.TPanedwindow3_p1.fig, axs=w.TPanedwindow3_p1.axs)
        w.TPanedwindow3_p1.canvas.draw()
    else:
        tkMessageBox.showwarning('Warning', 'Material not selected')

def btChangeMatColClick(master=None):
    global w, root
    tree = w.stvMaterial
    sel = tree.selection()
    if sel:
        key = tree.item(sel[0], 'text')
        _, res = askcolor(color=materials[key][1], parent=root, title='Color for material %s'%key)  #, alpha=True)
        if res:
            materials[key][1] = res
            update_spheres_canvas()

def add_material(key, material):
    global materials
    if key in materials:
       materials[key][0] = material
       sync_spheres_materials()
    else:
        color = w.color_pool.next()
        materials[key] = [material, color]

def find_mat_key(material):
    mat_name = str(material)
    try:
        key = materials.keys()[[str(m[0]) for m in materials.values()].index(mat_name)]
    except:
        # tkMessageBox.showerror('Error', 'Material key error for "%s".' % mat_name)
        return None
    return key

def gen_mat_key(ask_replace=False):
    """ generate new key for material dict """
    global w, materials
    if not ask_replace:
        i = len(materials)
        while "m%i"%i in materials:
            i += 1
        return "m%i"%i
    else:
        tree = w.stvMaterial
        sel = tree.selection()
        if sel:
            key = tree.item(sel[0], 'text')
            result = tkMessageBox.askquestion('Replace', 'Replace material %s?' % key)
            if result == 'yes':
                return key
        return gen_mat_key(False)

def update_materials_tree():
    global w, materials
    tree = w.stvMaterial
    tree.delete(*tree.get_children())
    for key in sorted(materials.iterkeys()):
        tree.insert('' , 0, text=key, values=(materials[key][0]))

    w.cbEnvMat.configure(values=materials.keys())
    if w.cbEnvMat.get() not in materials.keys():
        w.cbEnvMat.current(0)

def get_matrix_material():
    return materials[w.cbEnvMat.get()][0]

def get_scale():
    global w
    try:
        scale = float(w.edSpecScale.get())
    except ValueError as err:
        tkMessageBox.showerror('Error', 'Bad scale value. \n %s' % err)
    return scale

def get_wavelengths():
    global w
    try:
        xmin = float(w.edLambdaMin.get())
        xmax = float(w.edLambdaMax.get())
        count = int(w.edLambdaCount.get())
    except ValueError as err:
        tkMessageBox.showerror('Error', 'Bad value. \n %s' % err)
    assert count > 0
    return np.linspace(xmin, xmax, count)

def fitter_callback(fitter, values):
    global w, spheres

    w.edSpecScale.delete(0, 'end')
    w.edSpecScale.insert(0, fitter.params['scale'].value)

    n = fitter.background.number_of_params()
    if n > 0:
        w.edBkg1.delete(0, 'end')
        w.edBkg1.insert(0, fitter.params['bkg0'].value)
    if n > 1:
        w.edBkg2.delete(0, 'end')
        w.edBkg2.insert(0, fitter.params['bkg1'].value)
    if n > 2:
        w.edBkg3.delete(0, 'end')
        w.edBkg3.insert(0, fitter.params['bkg2'].value)

    if spheres is not None:
        spheres = copy.deepcopy(fitter.spheres)  # copy is not enough :(
        update_spheres_tree()
        update_spheres_canvas()

    btPlotExpClick()

    w.lbChiSq['text'] = 'ChiSq: %.6f' % fitter.chisq

def create_fitter(wls, fn, bkg):
    fitter = Fitter(fn, wl_min=wls.min(), wl_max=wls.max(),
                    wl_npoints=len(wls), bkg_method=bkg,
                    plot_progress=False)
    fitter.set_callback(fitter_callback)
    return fitter


def initialize_plot(widget):
    global fig, axs, canvas
    #~ f = Figure(figsize=(5, 4), dpi=100)
    widget.fig = Figure(dpi=75)
    widget.axs = widget.fig.add_subplot(111)
    widget.canvas = FigureCanvasTkAgg(widget.fig, master=widget)
    widget.canvas.show()
    widget.toolbar_frame = Frame(widget)
    widget.toolbar_frame.pack(side='top', fill='x')
    widget.toolbar_frame.toolbar = NavigationToolbar2TkAgg(widget.canvas, widget.toolbar_frame)
    widget.toolbar_frame.toolbar.update()
    widget.canvas.get_tk_widget().pack(side='top', fill='both', expand=False)
    widget.canvas.draw()

def TODO():
    print('test_paned_support.TODO')
    sys.stdout.flush()

def init(top, gui, *args, **kwargs):
    global w, top_level, root
    w = gui
    top_level = top
    root = top
    initialize_plot(w.TPanedwindow3_p1)
    w.canvas.camera = Camera()
    w.color_pool = cycle(['aqua', 'lime', 'blue', 'fuchsia', 'green', 'maroon',
                        'orange', 'pink', 'purple', 'red', 'yellow',
                        'violet', 'indigo', 'chartreuse', '#f55c4b'])
    cbBkgMethodSelect()
    w.model_fn = 'extinction.txt'
    w.constr_win = Toplevel(root)
    w.constr_win_app = ConstraintsWindow(w.constr_win)
    w.constr_win.withdraw()

def destroy_window():
    # Function which closes the window.
    global top_level
    if fitter is not None:
        fitter.stop()
        fitter.join()
    top_level.destroy()
    top_level = None


class SphereDialog(tkSimpleDialog.Dialog):
    def __init__(self, master, data_a=10, data_x=0, data_y=0, data_z=0, data_mat='m0'):
        self.data_x = data_x
        self.data_y = data_y
        self.data_z = data_z
        self.data_a = data_a
        self.data_mat = data_mat
        tkSimpleDialog.Dialog.__init__(self, master)

    def body(self, master):
        Label(master, text="X:").grid(row=1)
        Label(master, text="Y:").grid(row=2)
        Label(master, text="Z:").grid(row=3)
        Label(master, text="Radius:").grid(row=0)
        Label(master, text="Material ID:").grid(row=4)

        self.eX = Entry(master)
        self.eY = Entry(master)
        self.eZ = Entry(master)
        self.eR = Entry(master)
        self.emat = Entry(master)

        self.eX.insert(0, self.data_x)
        self.eY.insert(0, self.data_y)
        self.eZ.insert(0, self.data_z)
        self.eR.insert(0, self.data_a)
        self.emat.insert(0, self.data_mat)

        self.eX.grid(row=1, column=1)
        self.eY.grid(row=2, column=1)
        self.eZ.grid(row=3, column=1)
        self.eR.grid(row=0, column=1)
        self.emat.grid(row=4, column=1)
        return self.eR  # initial focus

    def validate(self):
        try:
            X = float(self.eX.get())
            Y = float(self.eY.get())
            Z = float(self.eZ.get())
            R = float(self.eR.get())
            mat_key = self.emat.get()
            assert mat_key in materials
            self.result = R, X, Y, Z, mat_key
            return True
        except ValueError as err:
            tkMessageBox.showerror("Error", str(err))
            return False

    def apply(self):
        pass


class GenerateSpheresDialog(tkSimpleDialog.Dialog):
    def __init__(self, master, data_N=0, data_a=10, data_d=5, materials=['m0']):
        self.data_N = data_N
        self.data_a = data_a
        self.data_d = data_d
        self.materials = materials
        tkSimpleDialog.Dialog.__init__(self, master)

    def body(self, master):
        Label(master, text='Number of spheres:').grid(row=1)
        Label(master, text='Spheres radius:').grid(row=2)
        Label(master, text='Gap between spheres:').grid(row=3)
        Label(master, text='Material ID:').grid(row=4)

        self.eN = Entry(master)
        self.eN.insert(0, self.data_N)
        self.ea = Entry(master)
        self.ea.insert(0, self.data_a)
        self.ed = Entry(master)
        self.ed.insert(0, self.data_d)
        self.cbmat = ttk.Combobox(master, values=self.materials)
        self.cbmat.insert(0, self.materials[-1])

        self.eN.grid(row=1, column=1)
        self.ea.grid(row=2, column=1)
        self.ed.grid(row=3, column=1)
        self.cbmat.grid(row=4, column=1)
        return self.eN  # initial focus

    def validate(self):
        try:
            N = int(self.eN.get())
            a = float(self.ea.get())
            d = float(self.ed.get())
            mat_key = self.cbmat.get()
            assert mat_key in materials
            self.result = N, a, d, mat_key
            return True
        except ValueError as err:
            tkMessageBox.showerror("Error", str(err))
            return False

    def apply(self):
        pass


class ConstraintsWindow:

    constr_types = ['Fix', 'Equality', 'Concentric']

    def __init__(self, master=None):
        self.master = master
        self.nspheres = 0
        master.title('Constraints')
        #~ master.geometry('600x450')
        #~ self.padWE = dict(sticky=('w', 'e'), padx='0.5mm', pady='0.5mm')
        self.padWE = dict(padx='0.5mm', pady='0.5mm')
        self.frame = ttk.Frame(self.master)
        self.create_widgets()
        self.configure_widgets()
        self.frame.grid(row=0, column=0)
        # binds
        self.count.trace('w', self.change_count)
        self.master.bind('<Configure>', self.configure_widgets)
        self.master.protocol('WM_DELETE_WINDOW', self.hide_window)
        self.master.bind('<Destroy>', self.hide_window)

    def create_widgets(self):
        self.count = StringVar()
        self.count.set('0')
        self.lbCount = ttk.Label(self.frame, text='Number of constraints:')
        self.sbCount = Spinbox(self.frame, from_=0, to=1000, textvariable=self.count, width=10)
        self.lbType = ttk.Label(self.frame, text='Type')
        self.lbPrm1 = ttk.Label(self.frame, text='Parameter#1')
        self.lbPrm2 = ttk.Label(self.frame, text='Parameter#2')
        self.cbTypes = []
        self.cbPrm1s = []
        self.cbPrm2s = []
        self.btOk = ttk.Button(self.frame, text='Ok', command=self.hide_window)
        self.btHelp = ttk.Button(self.frame, text='Help', command=self.show_help)

    def configure_widgets(self, event=None):
        self.lbCount.grid(row=0, column=0, columnspan=2, **self.padWE)
        self.sbCount.grid(row=0, column=2, columnspan=2, **self.padWE)
        self.lbType.grid(row=1, column=0, **self.padWE)
        self.lbPrm1.grid(row=1, column=1, **self.padWE)
        self.lbPrm2.grid(row=1, column=2, **self.padWE)
        n = int(self.count.get())
        for i in xrange(n):
            self.cbTypes[i].grid(row=2+i, column=0, **self.padWE)
            self.cbPrm1s[i].grid(row=2+i, column=1, **self.padWE)
            self.cbPrm2s[i].grid(row=2+i, column=2, **self.padWE)
        self.btOk.pack_forget()
        self.btOk.grid(row=3+n, column=0, **self.padWE)
        self.btHelp.grid(row=3+n, column=2, **self.padWE)

    def hide_window(self, event=None):
        if (event is not None) and (event.widget != self.master):
            return  # skip events from destruction of widgets
        self.master.withdraw()

    def show_window(self, nspheres):
        assert nspheres >= 0
        self.nspheres = nspheres
        print('  Number of spheres passed to Constrains window: %i' % nspheres)
        self.master.deiconify()

    def change_count(self, var, blank, mode):
        n = int(self.count.get())
        while len(self.cbTypes) > n:
            self.cbTypes[-1].destroy()
            self.cbTypes.pop()
            self.cbPrm1s[-1].destroy()
            self.cbPrm1s.pop()
            self.cbPrm2s[-1].destroy()
            self.cbPrm2s.pop()
        while len(self.cbTypes) < n:
            i = len(self.cbTypes)
            self.cbTypes.append(ttk.Combobox(self.frame, width=10))
            self.cbTypes[i].configure(values=self.constr_types)
            self.cbTypes[i].bind('<<ComboboxSelected>>', lambda event: self.select_type(event, i))
            self.cbPrm1s.append(ttk.Combobox(self.frame, width=10))
            self.cbPrm2s.append(ttk.Combobox(self.frame, width=10))
        print(self.cbTypes)
        self.configure_widgets()

    def select_type(self, event, irow):
        stype = self.cbTypes[irow].get()
        prms = ['scale', 'bkg0', 'bkg1', 'bkg2']
        for i in xrange(self.nspheres):
            prms.append('a%i'%i)
            prms.append('x%i'%i)
            prms.append('y%i'%i)
            prms.append('z%i'%i)
        if stype == 'Fix':
            self.cbPrm1s[irow].configure(values=prms)
            self.cbPrm2s[irow].configure(values=[])
        elif stype == 'Equality':
            self.cbPrm1s[irow].configure(values=prms)
            self.cbPrm2s[irow].configure(values=prms)
        elif stype == 'Concentric':
            prms = ['s%i' % i for i in xrange(self.nspheres)]
            self.cbPrm1s[irow].configure(values=prms)
            self.cbPrm2s[irow].configure(values=prms)
        else:
            raise Exception('Unknonw Constraint: "%s"' % stype)
        self.cbPrm1s[irow].delete(0, 'end')
        self.cbPrm2s[irow].delete(0, 'end')

    def get_constraints_list(self):
        n = int(self.count.get())
        result = []
        for i in xrange(n):
            stype = self.cbTypes[i].get()
            if stype == 'Fix':
                p1 = self.cbPrm1s[i].get()
                print('  %sConstraint(%s)' % (stype, p1))
                result.append(FixConstraint(p1))
            elif stype == 'Equality':
                p1 = self.cbPrm1s[i].get()
                p2 = self.cbPrm2s[i].get()
                print('  %sConstraint(%s, %s)' % (stype, p1, p2))
                result.append(EqualityConstraint(p1, p2))
            elif stype == 'Concentric':
                i1 = int(self.cbPrm1s[i].get()[1:])
                i2 = int(self.cbPrm2s[i].get()[1:])
                print('  %sConstraint(%i, %i)' % (stype, i1, i2))
                result.append(ConcentricConstraint(i1, i2))
            else:
                raise Exception('Unknonw Constraint: "%s"' % stype)
        return result

    def show_help(self):
        tkMessageBox.showinfo('Constraints Help',
        '''Types of constraints
        *Fix* -- don't vary the parameter during fitting. Can be applied to any parameter.
        *Equality* -- keep the values of to parameters equal. Can be applied to any parameters pair.
        *Concentric* -- maintain the centers of two spheres at the same position. This position is still varied. Can be applied to spheres pair only.
        ''')


class Camera(object):
    def __init__(self, scale=1, viewpoint=(0,0,0), alpha=0, beta=0):
        self.scale = scale
        self.viewpoint = np.array(viewpoint)
        self.alpha = alpha
        self.beta = beta
        self.axes = np.identity(3)

    def project_point(self, x,y,z):
        """
        Apply shift and project 3D positions to screen
        """
        X = (np.array([x, y, z]) - self.viewpoint) * self.scale
        X = np.dot(X, self.axes)
        #self.indices = X[:, 2].argsort()
        return X

    def project(self, positions):
        """
        Apply shift and project 3D positions to screen
        """
        X = (positions - self.viewpoint) * self.scale
        X = np.dot(X, self.axes)
        return X

    def zoom_in(self):
        self.scale *= 1.25

    def zoom_out(self):
        self.scale *= 1/1.25


class SplashWindow(Toplevel):
    def __init__(self, master):
        Toplevel.__init__(self, master)
        self.title('Splash')

        ## required to make window show before the program gets to the mainloop
        self.update()

if __name__ == '__main__':
    import mstm_studio
    mstm_studio.vp_start_gui()


























