#! /usr/bin/env python
#
# GUI module was first generated by PAGE version 4.9
# In conjunction with Tcl version 8.6
#    Dec 04, 2017 01:09:15 AM
# manually edited afterwards
import sys
import os
try:
    from Tkinter import Tk, Toplevel, Canvas, Menu, Pack, Grid, Place
except ImportError:
    from tkinter import Tk, Toplevel, Canvas, Menu, Pack, Grid, Place
from PIL import Image, ImageTk
try:
    import ttk
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    py3 = 1

import mstm_studio_support as sup

def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = Tk()
    top = MSTM_studio (root)
    sup.init(root, top)
    root.mainloop()

w = None
def create_MSTM_studio(root, *args, **kwargs):
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = Toplevel (root)
    top = MSTM_studio (w)
    sup.init(w, top, *args, **kwargs)
    return (w, top)

def destroy_MSTM_studio():
    global w
    w.destroy()
    w = None


class MSTM_studio:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        top.withdraw()
        self.splash = sup.SplashWindow(top)
        self.style = ttk.Style()
        if sys.platform == 'win32':
            self.style.theme_use('winnative')
        self.style.configure('.',font='TkDefaultFont')

        #~ top.geometry('838x455+364+117')
        top.geometry('850x435')
        top.title('MSTM studio')
        #~ top.configure(highlightcolor='black')
        self.load_images()
        self.TPanedwindow1 = ttk.Panedwindow(top, orient='horizontal')
        self.TPanedwindow1.place(relx=0.0, rely=0.0, relheight=1.0, relwidth=1.0)

        self.TPanedwindow1.configure(width=200)
        self.TPanedwindow1_p1 = ttk.Frame(width=220.0)
        self.TPanedwindow1.add(self.TPanedwindow1_p1)
        self.TPanedwindow1_p2 = ttk.Labelframe(width=350, text='View')
        self.TPanedwindow1.add(self.TPanedwindow1_p2)
        self.TPanedwindow1_p3 = ttk.Frame()
        self.TPanedwindow1.add(self.TPanedwindow1_p3)
        self.__funcid0 = self.TPanedwindow1.bind('<Map>', self.__adjust_sash0)

        self.TPanedwindow2 = ttk.Panedwindow(self.TPanedwindow1_p1, orient='vertical')
        self.TPanedwindow2.place(relx=0.0, rely=0.0, relheight=1.0, relwidth=1.0)

        self.TPanedwindow2.configure(width=200)
        self.TPanedwindow2_p1 = ttk.Labelframe(height=145, text='Materials')
        self.TPanedwindow2.add(self.TPanedwindow2_p1)
        self.TPanedwindow2_p2 = ttk.Labelframe(text='Spheres')
        self.TPanedwindow2.add(self.TPanedwindow2_p2)
        self.__funcid1 = self.TPanedwindow2.bind('<Map>', self.__adjust_sash1)

        self.style.configure('Treeview.Heading',  font='TkDefaultFont')
        self.stvMaterial = ScrolledTreeView(self.TPanedwindow2_p1)
        self.stvMaterial.place(relx=0.0, y=30, relheight=0.9, relwidth=1.0)
        self.stvMaterial.configure(columns='Col1')
        self.stvMaterial.heading('#0',text='MatID')
        self.stvMaterial.heading('#0',anchor='center')
        self.stvMaterial.column('#0',width='46')
        self.stvMaterial.column('#0',minwidth='20')
        self.stvMaterial.column('#0',stretch='1')
        self.stvMaterial.column('#0',anchor='w')
        self.stvMaterial.heading('Col1',text='Name')
        self.stvMaterial.heading('Col1',anchor='center')
        self.stvMaterial.column('Col1',width='150')
        self.stvMaterial.column('Col1',minwidth='20')
        self.stvMaterial.column('Col1',stretch='1')
        self.stvMaterial.column('Col1',anchor='w')
        self.stvMaterial.bind('<Double-1>', sup.btChangeMatColClick)

        self.btAddMat = ttk.Button(self.TPanedwindow2_p1, command=sup.btAddMatClick,
                                   text='A', image=self.imAdd)
        self.btAddMat.place(x=5, y=0, height=25, width=25)

        self.btLoadMat = ttk.Button(self.TPanedwindow2_p1, command=sup.btLoadMatClick,
                                    text='L', image=self.imLoad)
        self.btLoadMat.place(x=30, y=0, height=25, width=25)

        self.btPlotMat = ttk.Button(self.TPanedwindow2_p1, command=sup.btPlotMatClick,
                                    text='P', image=self.imPlot)
        self.btPlotMat.place(x=55, y=0, height=25, width=25)

        self.btDelMat = ttk.Button(self.TPanedwindow2_p1, command=sup.btDelMatClick,
                                   text='D', image=self.imDelete)
        self.btDelMat.place(relx=1, x=-30, rely=0, height=25, width=25)

        self.stvSpheres = ScrolledTreeView(self.TPanedwindow2_p2)
        self.stvSpheres.place(relx=0.0, y=30, relheight=0.9, relwidth=1.0)
        self.stvSpheres.configure(columns='Col1 Col2 Col3 Col4 Col5')
        self.stvSpheres.heading('#0', text='ID')
        self.stvSpheres.heading('#0', anchor='center')
        self.stvSpheres.column('#0', width='34')
        self.stvSpheres.column('#0', minwidth='20')
        self.stvSpheres.column('#0', stretch='1')
        self.stvSpheres.column('#0', anchor='w')
        self.stvSpheres.heading('Col1', text='R')
        self.stvSpheres.heading('Col1', anchor='center')
        self.stvSpheres.column('Col1', width='37')
        self.stvSpheres.column('Col1', minwidth='20')
        self.stvSpheres.column('Col1', stretch='1')
        self.stvSpheres.column('Col1', anchor='w')
        self.stvSpheres.heading('Col2', text='X')
        self.stvSpheres.heading('Col2', anchor='center')
        self.stvSpheres.column('Col2', width='31')
        self.stvSpheres.column('Col2', minwidth='20')
        self.stvSpheres.column('Col2', stretch='1')
        self.stvSpheres.column('Col2', anchor='w')
        self.stvSpheres.heading('Col3', text='Y')
        self.stvSpheres.heading('Col3', anchor='center')
        self.stvSpheres.column('Col3', width='34')
        self.stvSpheres.column('Col3', minwidth='20')
        self.stvSpheres.column('Col3', stretch='1')
        self.stvSpheres.column('Col3', anchor='w')
        self.stvSpheres.heading('Col4', text='Z')
        self.stvSpheres.heading('Col4', anchor='center')
        self.stvSpheres.column('Col4', width='34')
        self.stvSpheres.column('Col4', minwidth='20')
        self.stvSpheres.column('Col4', stretch='1')
        self.stvSpheres.column('Col4', anchor='w')
        self.stvSpheres.heading('Col5', text='mID')
        self.stvSpheres.heading('Col5', anchor='center')
        self.stvSpheres.column('Col5', width='32')
        self.stvSpheres.column('Col5', minwidth='20')
        self.stvSpheres.column('Col5', stretch='1')
        self.stvSpheres.column('Col5', anchor='w')
        self.stvSpheres.bind('<Double-1>', sup.btEditSphClick)

        self.btAddSph = ttk.Button(self.TPanedwindow2_p2, command=sup.btAddSphClick,
                                   text='A', image=self.imAdd)
        self.btAddSph.place(x=5, y=0, height=25, width=25)

        self.btEditSph = ttk.Button(self.TPanedwindow2_p2, command=sup.btEditSphClick,
                                    text='E', image=self.imEdit)
        self.btEditSph.place(x=30, y=0, height=25, width=25)

        self.btPlotSph = ttk.Button(self.TPanedwindow2_p2, command=sup.btPlotSphClick,
                                    text='R', image=self.imRefresh)
        self.btPlotSph.place(x=55, y=0, height=25, width=25)

        self.btDelSph = ttk.Button(self.TPanedwindow2_p2, command=sup.btDelSphClick,
                                   text='D', image=self.imDelete)
        self.btDelSph.place(relx=1.0, y=0, x=-30, height=25, width=25)

        self.canvas = Canvas(self.TPanedwindow1_p2)
        self.canvas.place(relx=0.0, rely=0, relheight=0.92, relwidth=1.0)
        self.canvas.configure(background='white')
        self.canvas.configure(borderwidth='2')
        self.canvas.configure(relief='ridge')
        self.canvas.configure(selectbackground='#c4c4c4')
        self.canvas.bind('<Button-4>',   sup.mouse_wheel)  # for Linux
        self.canvas.bind('<Button-5>',   sup.mouse_wheel)  # for Linux
        self.canvas.bind('<MouseWheel>', sup.mouse_wheel)  # for Windowz
        self.canvas.bind('<Button-3>',   sup.mouse_down)
        self.canvas.bind('<B3-Motion>',  sup.mouse_move)
        self.canvas.bind('<ButtonRelease-3>',  sup.mouse_up)

        self.lbZoom = ttk.Label(self.TPanedwindow1_p2, text='x1.00') #font=('courier', 18, 'bold'), width=10)
        self.lbZoom.place(relx=1.0, x=-50, rely=1.0, y=-25)

        self.TPanedwindow3 = ttk.Panedwindow(self.TPanedwindow1_p3, orient='vertical')
        self.TPanedwindow3.place(relx=0.0, rely=0.0, relheight=1.0, relwidth=1.0)

        self.TPanedwindow3.configure(width=200)
        self.TPanedwindow3_p1 = ttk.Labelframe(height=200, text='Plot')
        self.TPanedwindow3.add(self.TPanedwindow3_p1)
        self.TPanedwindow3_p2 = ttk.Labelframe(height=90, text='Spectrum')
        self.TPanedwindow3.add(self.TPanedwindow3_p2)
        self.TPanedwindow3_p3 = ttk.Labelframe(height=50, text='Background')
        self.TPanedwindow3.add(self.TPanedwindow3_p3)
        self.TPanedwindow3_p4 = ttk.Labelframe(text='Fitting')
        self.TPanedwindow3.add(self.TPanedwindow3_p4)
        self.__funcid2 = self.TPanedwindow3.bind('<Map>', self.__adjust_sash2)

        # Spectrum pane
        self.lbLambdaMin = ttk.Label(self.TPanedwindow3_p2, text='min')
        self.lbLambdaMin.place(x=5, y=0)
        self.edLambdaMin = ttk.Entry(self.TPanedwindow3_p2)
        self.edLambdaMin.place(x=5, y=15, width=35)
        self.edLambdaMin.insert(0, '300')

        self.lbLambdaMin = ttk.Label(self.TPanedwindow3_p2, text='max')
        self.lbLambdaMin.place(x=45, y=0)
        self.edLambdaMax = ttk.Entry(self.TPanedwindow3_p2)
        self.edLambdaMax.place(x=45, y=15, width=35)
        self.edLambdaMax.insert(0, '800')

        self.lbLambdaCount = ttk.Label(self.TPanedwindow3_p2, text='count')
        self.lbLambdaCount.place(x=85, y=0)
        self.edLambdaCount = ttk.Entry(self.TPanedwindow3_p2)
        self.edLambdaCount.place(x=85, y=15, width=35)
        self.edLambdaCount.insert(0, '51')

        self.lbSpecScale = ttk.Label(self.TPanedwindow3_p2, text='Scale')
        self.lbSpecScale.place(relx=1, x=-115, y=0)
        self.edSpecScale = ttk.Entry(self.TPanedwindow3_p2)
        self.edSpecScale.place(relx=1, x=-115, y=15, width=50)
        self.edSpecScale.insert(0, '1')

        self.lbEnvMat = ttk.Label(self.TPanedwindow3_p2, text='Matrix')
        self.lbEnvMat.place(relx=1, x=-60, y=0)
        self.cbEnvMat = ttk.Combobox(self.TPanedwindow3_p2)
        self.cbEnvMat.place(relx=1, x=-60, y=15, width=55)

        self.btCalcSpec = ttk.Button(self.TPanedwindow3_p2, command=sup.btCalcSpecClick,
                                     text='Calculate', image=self.imCalc, compound='left')
        self.btCalcSpec.place(x=5, y=40, width=90, height=25)

        self.btSaveSpec = ttk.Button(self.TPanedwindow3_p2, command=sup.btSaveSpecClick,
                                     text='S', image=self.imSave)
        self.btSaveSpec.place(relx=1, x=-55, y=40, width=25, height=25)

        self.btPlotSpec = ttk.Button(self.TPanedwindow3_p2, command=sup.btPlotSpecClick,
                                     text='P', image=self.imPlot)
        self.btPlotSpec.place(relx=1, x=-30, y=40, width=25, height=25)

        # Background pane
        self.cbBkgMethod = ttk.Combobox(self.TPanedwindow3_p3)
        self.BkgMethod_list = ['Constant', 'Linear', 'Lorentz']
        self.cbBkgMethod.configure(values=self.BkgMethod_list)
        self.cbBkgMethod.current(0)
        self.cbBkgMethod.place(x=5, y=5, width=80)
        self.cbBkgMethod.bind('<<ComboboxSelected>>',  sup.cbBkgMethodSelect)

        self.edBkg1 = ttk.Entry(self.TPanedwindow3_p3)
        self.edBkg1.place(x=85, y=5, width=45)
        self.edBkg1.insert(0, '0')

        self.edBkg2 = ttk.Entry(self.TPanedwindow3_p3)
        self.edBkg2.place(x=85+45, y=5, width=45)
        self.edBkg2.insert(0, '0')

        self.edBkg3 = ttk.Entry(self.TPanedwindow3_p3)
        self.edBkg3.place(x=85+45+45, y=5, width=45)
        self.edBkg3.insert(0, '0')

        self.btPlotBkg = ttk.Button(self.TPanedwindow3_p3, command=sup.btPlotBkgClick,
                                    text='P', image=self.imPlot)
        self.btPlotBkg.place(relx=1.0, x=-30, y=0, height=25, width=25)

        # Fitting pane
        self.edExpFileName = ttk.Entry(self.TPanedwindow3_p4, text='Exp. file name')
        self.edExpFileName.place(x=5, y=0, height=25, relwidth=0.8)

        self.btLoadExp = ttk.Button(self.TPanedwindow3_p4, command=sup.btLoadExpClick,
                                    text='L', image=self.imLoad)
        self.btLoadExp.place(relx=1.0, x=-55, y=0, height=25, width=25)

        self.btPlotExp = ttk.Button(self.TPanedwindow3_p4, command=sup.btPlotExpClick,
                                    text='P', image=self.imPlot)
        self.btPlotExp.place(relx=1.0, x=-30, y=0, height=25, width=25)

        self.btStartFit = ttk.Button(self.TPanedwindow3_p4, command=sup.btStartFitClick,
                                     text='>', image=self.imPlay)
        self.btStartFit.place(x=5, y=30, height=25, width=25)

        self.btStopFit = ttk.Button(self.TPanedwindow3_p4, command=sup.btStopFitClick,
                                    text='|', image=self.imStop)
        self.btStopFit.place(x=30, y=30, height=25, width=25)

        self.lbChiSq = ttk.Label(self.TPanedwindow3_p4, text='ChiSq:')
        self.lbChiSq.place(x=60, y=35)

        self.btConstraints = ttk.Button(self.TPanedwindow3_p4, command=sup.btConstraintsClick, text='Constraints...')
        self.btConstraints.place(relx=1, x=-90, y=30, height=25, width=85)

        self._create_menu(top)

        top.deiconify()
        self.splash.withdraw()

    def load_images(self):
        def tryload(fn):
            try:
                #im = PhotoImage(file=os.path.normpath(os.path.join('images', fn)))
                im = ImageTk.PhotoImage(Image.open(os.path.join('images', fn)))
            except Exception as err:
                print('Can not load %s\n%s' % (fn, err))
                return None
            return im
        self.imLoad    = tryload('folder_open_icon&16.png')
        self.imDelete  = tryload('delete_icon&16.png')
        self.imPlot    = tryload('chart_bar_icon&16.png')
        self.imAdd     = tryload('sq_plus_icon&16.png')
        self.imSave    = tryload('save_icon&16.png')
        self.imRefresh = tryload('refresh_icon&16.png')
        self.imExport  = tryload('export_icon&16.png')
        self.imImport  = tryload('import_icon&16.png')
        self.imPlay    = tryload('playback_play_icon&16.png')
        self.imStop    = tryload('playback_stop_icon&16.png')
        self.imCalc    = tryload('cogs_icon&16.png')
        self.imEdit    = tryload('doc_edit_icon&16.png')
        self.imExit    = tryload('on-off_icon&16.png')
        self.imBrush   = tryload('brush_icon&16.png')
        self.imZoomIn  = tryload('round_plus_icon&16.png')
        self.imZoomOut = tryload('round_minus_icon&16.png')

    def __adjust_sash0(self, event):
        paned = event.widget
        pos = [220, 575, ]
        i = 0
        for sash in pos:
            paned.sashpos(i, sash)
            i += 1
        paned.unbind('<map>', self.__funcid0)
        del self.__funcid0

    def __adjust_sash1(self, event):
        paned = event.widget
        pos = [145, ]
        i = 0
        for sash in pos:
            paned.sashpos(i, sash)
            i += 1
        paned.unbind('<map>', self.__funcid1)
        del self.__funcid1

    def __adjust_sash2(self, event):
        paned = event.widget
        pos = [200, ]
        i = 0
        for sash in pos:
            paned.sashpos(i, sash)
            i += 1
        paned.unbind('<map>', self.__funcid2)
        del self.__funcid2

    def _create_menu(self, top):
        self.menubar = Menu(top)

        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label='Import spheres...', command=sup.btImportSpheres,
                                  image=self.imImport, compound='left')
        self.filemenu.add_command(label='Export spheres...', command=sup.btExportSpheres,
                                  image=self.imExport, compound='left')
        self.filemenu.add_separator()
        self.filemenu.add_command(label='Exit', command=sup.destroy_window,
                                  image=self.imExit, compound='left')
        self.menubar.add_cascade(label='File', menu=self.filemenu)

        self.matmenu = Menu(self.menubar, tearoff=0)
        self.matmenu.add_command(label='Add constant...', command=sup.btAddMatClick,
                                 image=self.imAdd, compound='left')
        self.matmenu.add_command(label='Load function...', command=sup.btLoadMatClick,
                                 image=self.imLoad, compound='left')
        self.matmenu.add_separator()
        self.matmenu.add_command(label='Delete selected', command=sup.btDelMatClick,
                                 image=self.imDelete, compound='left')
        self.matmenu.add_separator()
        self.matmenu.add_command(label='Plot selected', command=sup.btPlotMatClick,
                                 image=self.imPlot, compound='left')
        self.matmenu.add_separator()
        self.matmenu.add_command(label='Change view color...', command=sup.btChangeMatColClick,
                                 image=self.imBrush, compound='left')
        self.menubar.add_cascade(label='Materials', menu=self.matmenu)

        self.sphmenu = Menu(self.menubar, tearoff=0)
        self.sphmenu.add_command(label='Add...', command=sup.btAddSphClick,
                                 image=self.imAdd, compound='left')
        self.sphmenu.add_separator()
        self.sphmenu.add_command(label='Edit selected...', command=sup.btEditSphClick,
                                 image=self.imEdit, compound='left')
        self.sphmenu.add_command(label='Delete selected', command=sup.btDelSphClick,
                                 image=self.imDelete, compound='left')
        self.sphmenu.add_separator()
        self.sphmenu.add_command(label='Generate on mesh...', command=sup.btGenerateSpheresClick)
        self.menubar.add_cascade(label='Spheres', menu=self.sphmenu)

        self.viewmenu = Menu(self.menubar, tearoff=0)
        self.viewmenu.add_command(label='Zoom in',
            command=lambda : sup.mouse_wheel(type('', (), {'num':4, 'delta':0})()),
            image=self.imZoomIn, compound='left')
        self.viewmenu.add_command(label='Zoom out',
            command=lambda : sup.mouse_wheel(type('', (), {'num':5, 'delta':0})()),
            image=self.imZoomOut, compound='left')
        self.viewmenu.add_command(label='Reset view', command=sup.btPlotSphClick,
                                  image=self.imRefresh, compound='left')
        self.menubar.add_cascade(label='View', menu=self.viewmenu)

        self.opticsmenu = Menu(self.menubar, tearoff=0)
        #~ self.opticsmenu.add_command(label='Setup...', command=sup.TODO)
        #~ self.opticsmenu.add_separator()
        self.opticsmenu.add_command(label='Calculate', command=sup.btCalcSpecClick,
                                    image=self.imCalc, compound='left')
        self.menubar.add_cascade(label='Spectrum', menu=self.opticsmenu)

        self.fittingmenu = Menu(self.menubar, tearoff=0)
        self.fittingmenu.add_command(label='Load experiment...', command=sup.btLoadExpClick,
                                     image=self.imLoad, compound='left')
        self.fittingmenu.add_separator()
        self.fittingmenu.add_command(label='Constraints...', command=sup.TODO)
        self.fittingmenu.add_separator()
        self.fittingmenu.add_command(label='Start fit', command=sup.btStartFitClick,
                                     image=self.imPlay, compound='left')
        self.fittingmenu.add_command(label='Stop fit', command=sup.btStopFitClick,
                                     image=self.imStop, compound='left')
        self.menubar.add_cascade(label='Fitting', menu=self.fittingmenu)

        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label='About', command=sup.TODO)
        self.menubar.add_cascade(label='Help', menu=self.helpmenu)
        # display the menu
        top.config(menu=self.menubar)




# The following code is added to facilitate the Scrolled widgets you specified.
class AutoScroll(object):
    '''Configure the scrollbars for a widget.'''

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)

        #self.configure(yscrollcommand=_autoscroll(vsb),
        #    xscrollcommand=_autoscroll(hsb))
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))

        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')

        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)

        # Copy geometry methods of master  (taken from ScrolledText.py)
        if py3:
            methods = Pack.__dict__.keys() | Grid.__dict__.keys() \
                  | Place.__dict__.keys()
        else:
            methods = Pack.__dict__.keys() + Grid.__dict__.keys() \
                  + Place.__dict__.keys()

        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        '''Hide and show scrollbar as needed.'''
        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)
        return wrapped

    def __str__(self):
        return str(self.master)

def _create_container(func):
    '''Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget.'''
    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        return func(cls, container, **kw)
    return wrapped

class ScrolledTreeView(AutoScroll, ttk.Treeview):
    '''A standard ttk Treeview widget with scrollbars that will
    automatically show/hide as needed.'''
    @_create_container
    def __init__(self, master, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)


if __name__ == '__main__':
    vp_start_gui()



