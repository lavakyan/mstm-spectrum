#! /usr/bin/env python
#
# Support module generated by PAGE version 4.9
# In conjunction with Tcl version 8.6
#    Dec 03, 2017 12:51:25 AM
#    Dec 03, 2017 09:18:10 PM
#    Dec 03, 2017 11:35:42 PM


import sys
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from itertools import cycle
import numpy as np

from mstm_spectrum import Material, Spheres, SingleSphere #, FilmBackground, LorentzBackground
from fit_spheres_optic import Fitter


try:
    from Tkinter import *
    from tkColorChooser import askcolor
except ImportError:
    from tkinter import *
    from tkinter.colorchooser import askcolor
import tkFileDialog, tkSimpleDialog, tkMessageBox
try:
    import ttk
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    py3 = 1








spheres = None
def btAddSphClick(master=None):
    global w, root
    global spheres
    dial = SphereDialog(root) # 10, 0.0, 0.0, 0.0, 'm0')
    if dial.result is None:
        return
    a, x, y, z, key = dial.result
    print(materials)
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

def btImportSpheres(master=None):
    global root, spheres, materials
    # open file dialog
    ftypes = [('Text files', '*.txt'), ('Data files', '*.dat'),
             ('Input files', '*.inp'), ('All files', '*')]
    dlg = tkFileDialog.Open(root, filetypes=ftypes)
    fn = dlg.show()
    if fn != '':
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
    #~ for i in xrange(len(spheres)):
    for i in indices:
        a = spheres.a[i] * cv.camera.scale
        x = W/2 + projected[i, 0]
        y = H/2 - projected[i, 1]
        key = find_mat_key(spheres.materials[i])
        col = materials[key][1]  # was '#5544FF'
        cv.create_oval(x-a, y-a, x+a, y+a, outline='#004500', width=3, fill=col, stipple='gray75')

def mouse_wheel(event):
    global w
    if event.num == 4 or event.delta == -120:    # Lin or Win
        w.canvas.camera.scale *= 1.25
    elif event.num == 5 or event.delta == 120:
        w.canvas.camera.scale *= 1/1.25
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
    global root
    ftypes = [('Text files', '*.txt'), ('All files', '*')]
    dlg = tkFileDialog.Open(root, filetypes=ftypes)
    fn = dlg.show()
    if fn:
        try:
            mat = Material(file_name=fn)
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
        axs.clear()
        mat.plot(fig=fig, axs=axs)
        canvas.draw()
    else:
        tkMessageBox.showwarning('Warning', 'Material not selected')

def btChangeMatColClick(master=None):
    global w, root
    tree = w.stvMaterial
    sel = tree.selection()
    if sel:
        key = tree.item(sel[0], 'text')
        res = askcolor(color=materials[key][1], parent=root, title='Color for material %s'%key)  #, alpha=True)
        if res:
            print(res[1])
            materials[key][1] = res[1]
            update_spheres_canvas()

def add_material(key, material):
    global materials
    if key in materials:
       materials[key][0] = material
    else:
        color = w.color_pool.next()
        print(color)
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
    global w
    tree = w.stvMaterial
    tree.delete(*tree.get_children())
    for key in sorted(materials.iterkeys()):
        tree.insert("" , 0, text=key, values=(materials[key][0]))

def initialize_plot(widget):
    global fig, axs, canvas
    #~ f = Figure(figsize=(5, 4), dpi=100)
    fig = Figure(dpi=75)
    axs = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=widget)
    canvas.show()
    widget.toolbar_frame = Frame(widget)
    widget.toolbar_frame.pack(side='top', fill='x')
    toolbar = NavigationToolbar2TkAgg(canvas, widget.toolbar_frame)
    toolbar.update()
    #~ canvas._tkcanvas.pack(side='top', fill='both', expand=False)
    canvas.get_tk_widget().pack(side='top', fill='both', expand=False)
    canvas.draw()

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
    w.color_pool = cycle(['aqua', 'blue', 'fuchsia', 'green', 'maroon',
                        'orange', 'pink', 'purple', 'red','yellow',
                        'violet', 'indigo', 'chartreuse', 'lime', '#f55c4b'])

def destroy_window():
    # Function which closes the window.
    global top_level
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

if __name__ == '__main__':
    import mstm_studio
    mstm_studio.vp_start_gui()


























