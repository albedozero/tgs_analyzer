#!/usr/bin/env python

VERSION_STR = 'TGSA v0.99 18-Apr-2013'

import os
from math import sin,cos,atan,pi,log
import numpy as np
from scipy.ndimage.interpolation import rotate
from scipy.signal import medfilt
import pyfits
import ConfigParser

import wx
import wx.html as html
import wx.grid as grid
import wx.aui as aui

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
import matplotlib.pyplot as plt
import matplotlib.font_manager
import matplotlib.figure
import matplotlib.gridspec as gridspec

# MKS Conversions
cm = 1.e-2
mm = 1.e-3
micron = 1.e-6
nm = 1.e-9
deg = pi/180.
arcmin = deg/60.
arcsec = deg/3600.

# User-configurable program settings
DEFAULT_APP_WIDTH=900	# default application window width
DEFAULT_APP_HEIGHT=650	# default application window height
PIXEL_START = 400		# start of spectrum from zeroth order
PIXEL_END = 1200		# end of spectrum from zeroth order
PIXEL_SCALE = 1.5		# nm to pixel scale
DEFAULT_WIDTH = 30		# default spectrum width in pixels
DEFAULT_WMIN = 350		# lower wavelength limit to plot
DEFAULT_WMAX = 750		# upper wavelength limit to plot
DEFAULT_WSPLIT = '483.3, 616.7'	# default wavelengths to stitch
STITCH_WIDTH = 3		# pixel window to use for scaling when stitching
REDSHIFT_Z = 0.0		# default redshift to apply
MEDAVG_WIDTH = 1		# default median averaging width
ZOOM_MAX = 4.0			# max zoom multiplier (power of two)
ZOOM_MIN = 0.03125		# min zoom multiplier (power of two)
TILT_INC = 0.25*deg		# spectrum tilt increment, in degrees
PLOT_BALMER = True
PLOT_HELIUM = False
PLOT_METALLIC = False
PLOT_TELLURIC = False
 
# Telescope and grating parameters
f_ratio = 14		# Telescope focal length
Diam = 37*cm        # Telescope diameter 
L = 38.8*mm   		# Distance from grating to CCD sensor  
lpmm = 600/mm		# Grating lines per mm
npixel = 2048		# Number of pixels along dispersion direction
pixel = 18*micron	# Pixel size

# Derived quantities
d_g = 1/lpmm
FL = Diam*f_ratio
x_ctr = npixel/2

myEVT_BROADCAST = wx.NewEventType()
EVT_BROADCAST = wx.PyEventBinder(myEVT_BROADCAST, 1)

class MsgEvent(wx.PyCommandEvent):
    def __init__(self, evtType, id):
        wx.PyCommandEvent.__init__(self, evtType, id)
        msg = None

    def setMsg(self, msg):
        self.msg = msg

    def getMsg(self):
        return self.msg
		
class FigPanel(wx.Panel):
	def __init__(self, *args, **kwargs):
		super(FigPanel, self).__init__(*args, **kwargs)
		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.fig = plt.Figure()
		self.canvas = FigureCanvas(self, -1, self.fig)
		self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
		self.SetSizer(self.sizer)
		self.toolbar = NavigationToolbar2Wx(self.canvas)
		self.toolbar.Hide()

		self.wmin = DEFAULT_WMIN
		self.wmax = DEFAULT_WMAX
		
	def plot(self, wave=[], intensity_wt=[]):
		self.fig.clear()
		if len(wave)>0:
			self.wave = wave
			self.intensity_wt = intensity_wt

		# run Median filter if requested
		if MEDAVG_WIDTH != 1:
			ampl = medfilt(self.intensity_wt,MEDAVG_WIDTH)
		else:
			ampl = self.intensity_wt
		
		if hasattr(self, 'spec'):
			gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
		else:
			gs = gridspec.GridSpec(1,1)
		ax1 = self.fig.add_subplot(gs[0])
		ax1.plot(self.wave,ampl, 'r-')
		ax1.set_xlabel('Wavelength [nm]')
		ax1.set_ylabel('Normalized intensity')

		# Set plot limits
		ax1.set_xlim(self.wmin, self.wmax)
		ymin = np.min(self.intensity_wt); ymax = 1.2*np.max(self.intensity_wt)
		ax1.set_ylim(ymin,ymax)

		# Plot reference lines if requested
		if PLOT_BALMER: plt_Balmer(ax1, REDSHIFT_Z)
		if PLOT_HELIUM: plt_Helium(ax1, REDSHIFT_Z)
		if PLOT_METALLIC: plt_Metals(ax1, REDSHIFT_Z)
		if PLOT_TELLURIC: plt_Telluric(ax1, REDSHIFT_Z)
		if PLOT_BALMER or PLOT_HELIUM or PLOT_METALLIC or PLOT_TELLURIC:
			ax1.legend(loc='upper right',prop={'size':10})

		# Add title
		if not hasattr(self, 'title'):
			title = 'Rigel TGS Spectrum Plot'
			ax1.set_title(title, fontsize=12)						
		else:
			ax1.set_title(self.title, fontsize=12)			
		ax1.grid(True)

		# Strip Spectrum
		if hasattr(self, 'spec'):
			ax2 = self.fig.add_subplot(gs[1])
			ax2.imshow(self.spec)
			ax2.set_xlabel("Pixel offset from zeroth order (+%d)" % (PIXEL_START))
			ax2.get_yaxis().set_visible(False)

		# Add version info in lower right corner
		ax1.text(0.85,-0.5,VERSION_STR,ha ='left',fontsize=8, transform = ax1.transAxes)
		
		self.canvas.draw()
		self.is_saved = False
		
class ImgViewer(wx.SplitterWindow):
	def __init__(self, path, data, hdr, *args, **kwargs):
		super(ImgViewer, self).__init__(*args, **kwargs)
		self.scroll = wx.ScrolledWindow(self)
		self.info = wx.grid.Grid(self)
		self.info.CreateGrid(1,3)
		self.info.SetColLabelSize(0)
		self.info.SetRowLabelSize(0)
		self.SplitVertically(self.info, self.scroll, 200)
		self.Unsplit(self.info)
	
		# Load the FITS image
		self.fname = os.path.basename(path)
		self.data = data
		self.hdr = hdr
		(w, h) = self.data.shape
		self.imin = np.min(self.data)
		self.imax = self.imin + 16.0*np.std(self.data)
		b = np.clip( 255.99*(self.data-self.imin)/(self.imax-self.imin),0,255.99 )
		self.img = wx.EmptyImage(w, h)
		self.img.SetData( np.dstack((b,b,b)).astype('uint8').tostring() )
		self.bmp = wx.BitmapFromImage(self.img)
		i=0
		for k in hdr.keys():
			if hdr.get(k,False):
				while i>=self.info.GetNumberRows():
					self.info.AppendRows()
				self.info.SetCellValue(i,0,"%s" % k)
				self.info.SetCellValue(i,1,"%s" % (hdr.get(k)) )
				if hdr.comments[k] != '':
					self.info.SetCellValue(i,2,"%s" % (hdr.comments[k]) )
				for j in range(3): self.info.SetReadOnly(i,j,True)
				i+=1
		self.info.AutoSize()

		# Object data
		self.zoomLevel = 1.0
		self.was_zoomed=False
		self.alignDegrees = 0.0
		xymax = np.unravel_index(np.argmax(self.data),self.data.shape)
		self.specZero = wx.Point(xymax[1],xymax[0])
		self.specWidth = DEFAULT_WIDTH
		
		# Event bindings
		self.scroll.Bind(wx.EVT_LEFT_DOWN, self.onLeftDown)
		self.scroll.Bind(wx.EVT_PAINT, self.onPaint)
		self.scroll.Bind(wx.EVT_MOTION, self.onMotion)
		self.scroll.Bind(wx.EVT_LEAVE_WINDOW, self.onLeave)
		self.scroll.Bind(wx.EVT_IDLE, self.zoom_redraw)

		self.scroll.SetVirtualSize((w,h))
		self.scroll.SetScrollRate(20,20)
		x,y=self.scroll.CalcScrolledPosition(xymax[1]-20,xymax[0]-20)
		dx,dy=self.scroll.GetScrollPixelsPerUnit()
		self.scroll.Scroll(x/dx,y/dy)
				
	def onMotion(self, event):
		pos = event.GetPosition()
		(x, y) = self.scroll.CalcUnscrolledPosition(pos.x, pos.y)
		info = "(%d, %d) %d%%" % (x/self.zoomLevel, y/self.zoomLevel, int(self.zoomLevel*100))
		event = MsgEvent(myEVT_BROADCAST, self.GetId())
		event.setMsg(info)
		self.GetEventHandler().ProcessEvent(event)

	def onLeave(self, event):
		event = MsgEvent(myEVT_BROADCAST, self.GetId())
		event.setMsg('')
		self.GetEventHandler().ProcessEvent(event)	
	
	def onLeftDown(self, event):
		pos = event.GetPosition()
		(x, y) = self.scroll.CalcUnscrolledPosition(pos.x, pos.y)
		self.specZero = wx.Point(x/self.zoomLevel, y/self.zoomLevel)
		self.scroll.Refresh()

	def onPaint(self, event):
		dc = wx.PaintDC(self.scroll)
		self.scroll.DoPrepareDC(dc)
		self.draw(dc)
		
	def adjustImage(self, min=-1, max=-1):
		if min<0: min=self.imin
		if max<0: max=self.imax
		if max==min:max+=1.0
		b = np.clip( 256.0*(self.data-min)/(max-min),0,255.99 )
		self.img.SetData( np.dstack((b,b,b)).astype('uint8').tostring() )
		self.zoom(mult=1)

	def draw(self, dc):
		dc.DrawBitmap(self.bmp,0,0,False)
		dc.SetBrush(wx.Brush('#000000', wx.TRANSPARENT))
		dc.SetPen(wx.Pen('RED', 1, wx.SOLID))
		x = self.specZero.x*self.zoomLevel
		y = self.specZero.y*self.zoomLevel
		cosa = cos(-self.alignDegrees)
		sina = sin(-self.alignDegrees)
		dw = 0.5*self.specWidth*self.zoomLevel
		p0 = PIXEL_START*self.zoomLevel
		pf = PIXEL_END*self.zoomLevel
		dc.DrawLine(x-10,y,x+10,y)
		dc.DrawLine(x,y-10,x,y+10)
		sel =  [ wx.Point( x + p0*cosa + dw*sina, y + p0*sina - dw*cosa ) ]
		sel += [ wx.Point( x + p0*cosa - dw*sina, y + p0*sina + dw*cosa ) ]
		sel += [ wx.Point( x + pf*cosa - dw*sina, y + pf*sina + dw*cosa ) ]
		sel += [ wx.Point( x + pf*cosa + dw*sina, y + pf*sina - dw*cosa ) ]
		sel += [ wx.Point( x + p0*cosa + dw*sina, y + p0*sina - dw*cosa ) ]
		dc.DrawLines(sel)

	def zoom(self, mult=1.0, zoom=1.0):
		if mult != 0:
			zoom = 2**(int(log(self.zoomLevel*mult,2)))
		if zoom>ZOOM_MAX:
			zoom=ZOOM_MAX
		if zoom<ZOOM_MIN:
			zoom=ZOOM_MIN
		w,h = self.img.GetWidth(), self.img.GetHeight()
		w *= zoom
		h *= zoom
		self.bmp = wx.BitmapFromImage(self.img.Scale(w, h))
		self.scroll.SetVirtualSize((w,h))
		self.scroll.SetScrollRate(20,20)
		self.scroll.Refresh(True)

		info = "Zoom level %d%%" % (int(self.zoomLevel*100))
		event = MsgEvent(myEVT_BROADCAST, self.GetId())
		event.setMsg(info)
		self.GetEventHandler().ProcessEvent(event)
		self.zoomLevel=zoom
		self.was_zoomed=True

	def zoom_redraw(self, event):
		# This is a hack for Mac OS X - selection box shows up in
		# old spot - maybe scrolledwindow virtual size isn't updating
		# until idle?
		if self.was_zoomed:
			self.Refresh(True)
			self.was_zoomed = False

	def zoom_fit(self):
		vw,vh = self.scroll.GetClientSize()
		w,h = self.img.GetWidth(), self.img.GetHeight()
		z = 1.*vw/w
		if 1.*vh/h<z: z=1.*vh/h
		if z>1.0: z=1.0
		self.zoom(0, z)

	def sel_width(self, inc):
		self.specWidth += inc
		if self.specWidth<1: self.specWidth=1
		self.Refresh()

	def sel_tilt(self, inc):
		self.alignDegrees += inc
		if self.alignDegrees<-60: self.alignDegrees=-60
		if self.alignDegrees> 60: self.alignDegrees= 60
		self.Refresh()
		
	def sel_nudge(self, dir):
		(dx, dy) = dir
		self.specZero.x += dx
		self.specZero.y += dy
		w,h = self.img.GetWidth(), self.img.GetHeight()
		if self.specZero.x > w: self.specZero.x = w
		if self.specZero.x < 0: self.specZero.x = 0
		if self.specZero.y > h: self.specZero.y = h
		if self.specZero.y < 0: self.specZero.y = 0
		self.Refresh()

	def extract(self):
		# Define subimage containing dispersed spectrum
		xmin = PIXEL_START + self.specZero.x
		xmax = PIXEL_END + self.specZero.x
		ymin  = self.specZero.y - self.specWidth/2.
		ymax = self.specZero.y + self.specWidth/2.
		w, h = self.data.shape
		if xmin>=w: return
		if xmax>=w: xmax=w-1
		if ymin<0: ymin=0
		if ymax>=h: ymax=h-1
		
		if self.alignDegrees != 0.0: im = rotate(self.data, -self.alignDegrees/deg)
		else: im = np.copy(self.data)
		spec = im[ymin:ymax,xmin:xmax]
		# Get 'off source' spectrum to find background
		yoff = self.specWidth
		spec_off = im[ymin+yoff:ymax+yoff,xmin:xmax]

		# Subtract median background from each pixel
		xmed =  np.median(spec_off,axis=0)
		spec -= xmed
		spec_sum = np.sum(spec,axis=0)

		# Fill arrays for plotting, using efficiency curves to correct for sensitivity vs wavelength
		npts = xmax - xmin
		wave = np.zeros(npts); intensity = np.zeros(npts)
		intensity_wt = np.zeros(npts)
		for n in range(npts):
			wave[n] = f_lambda(xmin-self.specZero.x+n,self.specZero.x)
			wt = ccd_sens(wave[n])*grating_sens(wave[n])
			intensity[n] = float(spec_sum[n])
			intensity_wt[n] = intensity[n]/wt
		intensity_wt /= max(intensity_wt)

		return (wave, intensity_wt, spec)
		
class MainWindow(wx.Frame):
	def __init__(self, filename='noname.txt'):
		super(MainWindow, self).__init__(None)
		super(MainWindow, self).SetTitle('TGS Analyzer')
		self.SetSize(size=wx.Size(DEFAULT_APP_WIDTH,DEFAULT_APP_HEIGHT))
#		_icon = wx.Icon('tgsa.ico', wx.BITMAP_TYPE_ICO)
#		self.SetIcon(_icon)

		# Make interior window components
		self.tabs = aui.AuiNotebook(self)
		self.CreateStatusBar()

		# Event Bindings
		self.Bind(EVT_BROADCAST, self.onBroadcast)
		self.Bind(wx.EVT_CLOSE, self.onExit)
		self.Bind(wx.EVT_CHAR_HOOK, self.onKeyPress)
		self.tabs.Bind(aui.EVT_AUINOTEBOOK_BUTTON, self.onClose)

		# Create Menus		
		fileMenu = wx.Menu()
		tmp = fileMenu.Append(wx.ID_ANY,'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = fileMenu.Append(wx.ID_OPEN, '&Open\tCtrl+O', 'Open FITS spectrum file(s)')
		self.Bind(wx.EVT_MENU, self.onOpen, item)
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN))
		fileMenu.Remove(tmp.GetId()) # deals with wxPython bug where first menu item with a bitmap doesn't show up
		item = fileMenu.Append(wx.ID_SAVEAS, '&Save Plot\tCtrl+S', 'Save spectrum plot')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE))
		self.Bind(wx.EVT_MENU, self.onSave, item)
		item = fileMenu.Append(wx.ID_ANY, '&Close\tCtrl+W', 'Close current window')
		self.Bind(wx.EVT_MENU, self.onClose, item)
		item = fileMenu.Append(wx.ID_ANY, '&Import Data', 'Load spectrum data from a text file')
		self.Bind(wx.EVT_MENU, self.onImport, item)
		item = fileMenu.Append(wx.ID_ANY, 'Export &Data', 'Export spectrum data as text')
		self.Bind(wx.EVT_MENU, self.onExport, item)
		fileMenu.AppendSeparator()
		item = fileMenu.Append(wx.ID_EXIT, 'E&xit\tCtrl+Q', 'Terminate the program')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_QUIT))
		self.Bind(wx.EVT_MENU, self.onExit, item)

		selMenu = wx.Menu()
		tmp = selMenu.Append(wx.ID_ANY,'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = selMenu.Append(wx.ID_ANY, '&Decrease width\tCtrl+Left', 'Decrease selection width')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_GO_BACK))
		selMenu.Remove(tmp.GetId())
		self.Bind(wx.EVT_MENU, lambda x: self.onWidth(x, -1), item)
		item = selMenu.Append(wx.ID_ANY, '&Increase width\tCtrl+Right', 'Increase selection width')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD))
		self.Bind(wx.EVT_MENU, lambda x: self.onWidth(x, 1), item)
		item = selMenu.Append(wx.ID_ANY, 'Tilt &up\tCtrl+Up', 'Tilt selection up')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_GO_UP))
		self.Bind(wx.EVT_MENU, lambda x: self.onTilt(x, TILT_INC), item)
		item = selMenu.Append(wx.ID_ANY, 'Tilt dow&n\tCtrl+Down', 'Tilt selection down')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN))
		self.Bind(wx.EVT_MENU, lambda x: self.onTilt(x, -TILT_INC), item)
		item = selMenu.Append(wx.ID_ANY, '&Set Position', 'Set zeroth order diffraction point')
		self.Bind(wx.EVT_MENU, self.onSetPosition, item)

		imgMenu = wx.Menu()
		tmp = imgMenu.Append(wx.ID_ANY, 'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = imgMenu.Append(wx.ID_ZOOM_IN, 'Zoom &in\tCtrl++', 'Zoom in')
		self.Bind(wx.EVT_MENU, lambda x: self.onZoom(x, 2.0), item)
		item = imgMenu.Append(wx.ID_ZOOM_OUT, 'Zoom &out\tCtrl+-', 'Zoom out')
		self.Bind(wx.EVT_MENU, lambda x: self.onZoom(x, 0.5), item)
		item = imgMenu.Append(wx.ID_ZOOM_100, 'Zoo&m 100%\tCtrl+0', 'Zoom to original size')
		self.Bind(wx.EVT_MENU, lambda x: self.onZoom(x, 0), item)
		item = imgMenu.Append(wx.ID_ANY, 'Zoom &Fit', 'Zoom to fit in window')
		self.Bind(wx.EVT_MENU, self.onZoomFit, item)
		imgMenu.AppendSeparator()
		imgMenu.AppendSubMenu(selMenu, '&Selection', 'Modify selection')
		item = imgMenu.Append(wx.ID_ANY, '&Adjust Image Levels...', 'Adjust brightness ranges in image')
		self.Bind(wx.EVT_MENU, self.onAdjImage, item)
		item = imgMenu.Append(wx.ID_ANY, 'Show/Hide FITS &Header\tCtrl+F', 'Toggle visibility of the FITS header data')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_REPORT_VIEW))
		imgMenu.Remove(tmp.GetId())
		self.Bind(wx.EVT_MENU, self.onShowFITS, item)

		linesMenu = wx.Menu()
		self.Balmer = linesMenu.Append(wx.ID_ANY, 'Balmer Series', 'Show hydrogen Balmer lines on plot', kind=wx.ITEM_CHECK)
		self.Helium = linesMenu.Append(wx.ID_ANY, 'Helium Lines', 'Show helium lines on plot', kind=wx.ITEM_CHECK)
		self.Metallic = linesMenu.Append(wx.ID_ANY, 'Metallic', 'Show metal lines on plot', kind=wx.ITEM_CHECK)
		self.Telluric = linesMenu.Append(wx.ID_ANY, 'Telluric', 'Show atmospheric absorption lines on plot', kind=wx.ITEM_CHECK)
		if PLOT_BALMER: self.Balmer.Check(True)
		if PLOT_HELIUM: self.Helium.Check(True)
		if PLOT_METALLIC: self.Metallic.Check(True)
		if PLOT_TELLURIC: self.Telluric.Check(True)		
		
		modMenu = wx.Menu()
		tmp = modMenu.Append(wx.ID_ANY, 'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = modMenu.Append(wx.ID_ANY, 'Change &Title...', 'Change title of current plot')
		self.Bind(wx.EVT_MENU, self.onSetTitle, item)
		item = modMenu.Append(wx.ID_ANY, 'Change &Range...', 'Change wavelength range of current plot')
		self.Bind(wx.EVT_MENU, self.onSetRange, item)
		item = modMenu.Append(wx.ID_ANY, '&Adjust Plot', 'Adjust size/margins of current plot')
		self.Bind(wx.EVT_MENU, self.onAdjust, item)
		modMenu.Remove(tmp.GetId())
		
		plotMenu = wx.Menu()
		tmp = plotMenu.Append(wx.ID_ANY, 'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = plotMenu.Append(wx.ID_ANY, 'Create/Redo &Plot', 'Plot spectrum from image or refresh current plot')
		self.Bind(wx.EVT_MENU, self.onPlot, item)
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_NEW))
		plotMenu.Remove(tmp.GetId())
		item = plotMenu.Append(wx.ID_ANY, '&Stitch Spectra...', 'Stitch multiple spectra together')
		self.Bind(wx.EVT_MENU, self.onStitch, item)
		plotMenu.AppendSubMenu(linesMenu, '&Lines', 'Select spectral lines to show')
		item = plotMenu.Append(wx.ID_ANY, 'Set Redshift(&z)...', 'Set redshift of lines')
		self.Bind(wx.EVT_MENU, self.onSetRedshift, item)
		item = plotMenu.Append(wx.ID_ANY, 'Set A&veraging...', 'Set median average width')
		self.Bind(wx.EVT_MENU, self.onSetAverage, item)
		plotMenu.AppendSubMenu(modMenu, '&Modify Plot', 'Modify the current plot')
		
		helpMenu = wx.Menu()
		tmp = helpMenu.Append(wx.ID_ANY,'X')
		tmp.SetBitmap(wx.EmptyBitmap(1,1))
		item = helpMenu.Append(wx.ID_HELP, '&Help\tF1', 'Get help with this program')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_QUESTION, size=(16,16)))
		helpMenu.Remove(tmp.GetId())
		self.Bind(wx.EVT_MENU, self.onHelp, item)
		item = helpMenu.Append(wx.ID_ABOUT, '&About', 'Information about this program')
		item.SetBitmap(wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, size=(16,16)))
		self.Bind(wx.EVT_MENU, self.onAbout, item)

		menuBar = wx.MenuBar()
		menuBar.Append(fileMenu, '&File')
		menuBar.Append(imgMenu, '&Image')
		menuBar.Append(plotMenu, '&Plot')
		menuBar.Append(helpMenu, '&Help')
		self.SetMenuBar(menuBar)

		# Make Toolbars
		fileTool = self.CreateToolBar()
		tool = fileTool.AddLabelTool(wx.ID_OPEN, 'Open', wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN), shortHelp='Open', longHelp='Open a file')
		self.Bind(wx.EVT_TOOL, self.onOpen, tool)
		tool = fileTool.AddLabelTool(wx.ID_ANY, 'Save', wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE), shortHelp='Save', longHelp='Save plot')
		self.Bind(wx.EVT_TOOL, self.onSave, tool)
		fileTool.AddSeparator()
		b = wx.ArtProvider.GetBitmap(wx.ART_GO_BACK)
		tool = fileTool.AddLabelTool(wx.ID_ANY, 'Decrease', b, shortHelp='Decrease width', longHelp='Decrease selection width')
		self.Bind(wx.EVT_TOOL, lambda x: self.onWidth(x, -1), tool)
		b = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD)
		tool = fileTool.AddLabelTool(wx.ID_ANY, 'Increase', b, shortHelp='Increase width', longHelp='Increase selection width')
		self.Bind(wx.EVT_TOOL, lambda x: self.onWidth(x, 1), tool)
		b = wx.ArtProvider.GetBitmap(wx.ART_GO_UP)
		tool = fileTool.AddLabelTool(wx.ID_ANY, 'Tilt Up', b, shortHelp='Tilt Up', longHelp='Tilt selection up')
		self.Bind(wx.EVT_TOOL, lambda x: self.onTilt(x, TILT_INC), tool)
		b = wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN)
		tool = fileTool.AddLabelTool(wx.ID_ANY, 'Tilt Down', b, shortHelp='Tilt Down', longHelp='Tilt selection down')
		self.Bind(wx.EVT_TOOL, lambda x: self.onTilt(x, -TILT_INC), tool)
		b = wx.ArtProvider.GetBitmap(wx.ART_NEW)
		tool = fileTool.AddLabelTool(wx.ID_APPLY, 'Plot', b, shortHelp='Plot', longHelp='Plot spectrum from image or refresh current plot')
		self.Bind(wx.EVT_TOOL, self.onPlot, tool)
		fileTool.Realize()
				
		# Create dialog boxes
		self.stitch = StitchDialog(self)
		self.stitch.Show(False)
		self.stitch.Bind(wx.EVT_BUTTON, self.do_stitch, id=wx.ID_OK)
		
		self.help = HelpDialog(self)
		self.help.Show(False)
				
	def readCheckItems(self):
		global PLOT_BALMER, PLOT_HELIUM, PLOT_METALLIC, PLOT_TELLURIC
		if self.Balmer.IsChecked(): PLOT_BALMER=True
		else: PLOT_BALMER=False
		if self.Helium.IsChecked(): PLOT_HELIUM=True
		else: PLOT_HELIUM=False
		if self.Metallic.IsChecked(): PLOT_METALLIC=True
		else: PLOT_METALLIC=False
		if self.Telluric.IsChecked(): PLOT_TELLURIC=True
		else: PLOT_TELLURIC=False		

	def onKeyPress(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer) and event.GetModifiers()<=0:
			if event.GetKeyCode() == wx.WXK_UP:
				cur.sel_nudge(( 0,-1))
				return
			if event.GetKeyCode() == wx.WXK_DOWN:
				cur.sel_nudge(( 0, 1))
				return
			if event.GetKeyCode() == wx.WXK_LEFT:
				cur.sel_nudge((-1, 0))
				return
			if event.GetKeyCode() == wx.WXK_RIGHT:
				cur.sel_nudge(( 1, 0))
				return
		event.Skip()
		
	def onOpen(self, event):
		wildcard = "FITS image files (*.fts,*.fits,*.fit)|*.fts;*.fits;*.fit"
		dialog = wx.FileDialog(None, "Choose a file", wildcard=wildcard, style=wx.FD_OPEN|wx.FD_MULTIPLE)
		if dialog.ShowModal() == wx.ID_OK:
			for path in dialog.GetPaths():
				fname = os.path.basename(path)
				try:
					data,hdr = pyfits.getdata(path,0,header=True)
				except:
					msg = 'Error opening %s' % (fname)
					errmsg = wx.MessageDialog(self, msg ,'File Error', style=wx.OK|wx.ICON_ERROR)
					errmsg.ShowModal()
					continue
				newim = ImgViewer(parent=self.tabs, path=path, data=data, hdr=hdr)
				self.tabs.AddPage(newim, fname, select=True)
		dialog.Destroy()

	def onImport(self, event):
		wildcard = "Text file (*.csv)|*.csv"
		dialog = wx.FileDialog(None, "Choose a file", wildcard=wildcard, style=wx.FD_OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			fname = os.path.basename(path)
			newplot = FigPanel(self.tabs)
			newplot.fname = fname
			self.tabs.AddPage(newplot, fname, select=True)
			fn = open(path,'r')
			lines = fn.readlines()
			wave = []; ampl = []
			for line in lines:
				if line[0] == '#': continue
				s = [float(t) for t in line.split()]
				wave.append(s[0]); ampl.append(s[1])
			fn.close()
			newplot.title = 'Rigel TGS Spectrum from %s' % (fname)
			newplot.plot(wave, ampl)
		dialog.Destroy()

	def onSave(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, FigPanel): return
		if not hasattr(cur, 'fname'): fname = '.pdf'
		else: fname = os.path.splitext(cur.fname)[0] + '-spec.pdf'
		wildcard = "PDF File (*.pdf)|*.pdf"
		dialog = wx.FileDialog(None, "Save Plot", defaultFile=fname, wildcard=wildcard, style=wx.FD_SAVE)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			cur.canvas.print_figure(path)
			cur.is_saved=True
			dialog.Destroy()
			return wx.ID_OK
		else:
			dialog.Destroy()
			return wx.ID_CANCEL

	def onExport(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, FigPanel): return
		if not hasattr(cur, 'fname'):
			fname = '.pdf'
			fromfile = 'unknown'
		else:
			csvname = os.path.splitext(cur.fname)[0] + '-data.csv'
			fromfile = cur.fname
		wildcard = "Text file (*.csv)|*.csv"
		dialog = wx.FileDialog(None, "Export Spectral Data", defaultFile=csvname, wildcard=wildcard, style=wx.FD_SAVE)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			fn = open(path,'w')
			fn.write('# File %s\n' % fromfile)
			fn.write('# Wavelength [nm]    Normalized Amplitude\n')
			for n in range(len(cur.wave)):
				fn.write('   %5.2f         %6.4f\n' % (cur.wave[n], cur.intensity_wt[n]))
			fn.close()
		dialog.Destroy()

	def onClose(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		p = self.tabs.GetPage(n)
		if isinstance(p, FigPanel) and not p.is_saved:
			name = self.tabs.GetPageText(n)
			dialog = wx.MessageDialog(self, 'Save %s before closing?' % (name), 'Close', wx.YES_NO|wx.CANCEL)
			response = dialog.ShowModal()
			dialog.Destroy()
			if response == wx.ID_CANCEL:
				return wx.ID_CANCEL
			if response == wx.ID_YES:
				saveok = self.onSave(None)
				if saveok == wx.ID_CANCEL:
					return wx.ID_CANCEL
			self.tabs.DeletePage(n)
		else:
			self.tabs.DeletePage(n)
			return wx.ID_OK

	def onStitch(self, event):
		if self.stitch.IsShown():
			self.stitch.Raise()
			return
		plotnames = {}
		for k in range(self.tabs.GetPageCount()):
			if isinstance(self.tabs.GetPage(k), FigPanel):
				id = self.tabs.GetPage(k).GetId()
				plotnames[id] = self.tabs.GetPageText(k)
		if len(plotnames.keys())<2: return
		self.stitch.repopulate(plotnames)
		self.stitch.Show(True)
		
	def do_stitch(self, event):
		ids, wsplits = self.stitch.get_stitched()
		ids = ids[:len(wsplits)+1] # in case user didn't specify enough splits
		for id in ids:
			plot = wx.FindWindowById(id)
			if plot == None: return
		waves = [wx.FindWindowById(id).wave for id in ids]
		vals = [wx.FindWindowById(id).intensity_wt for id in ids]
		wb = [waves[0][0]] + wsplits[0:] + [waves[-1][-2]]
		nf = [j for j in range(len(waves[0])) if waves[0][j]>wb[1]][0]
		wave_stitch = waves[0][0:nf]
		intensity_stitch = vals[0][0:nf]
		for i in range(len(wb)-2):
			# get indices of array segments to stitch
			nf = [j for j in range(len(waves[i])) if waves[i][j]>wb[i+1]][0]
			mi = [j for j in range(len(waves[i+1])) if waves[i+1][j]>wb[i+1]][0]
			mf = [j for j in range(len(waves[i+1])) if waves[i+1][j]>wb[i+2]][0]
			# get median value a few pixels on each side of stitch for scaling
			scale = np.median(vals[i][nf-STITCH_WIDTH:nf])/np.median(vals[i+1][mi:mi+STITCH_WIDTH])
			vals[i+1] = [x * scale for x in vals[i+1]]
			x = waves[i+1][mi:mf]
			wave_stitch = np.concatenate((wave_stitch, waves[i+1][mi:mf]))
			intensity_stitch = np.concatenate((intensity_stitch, vals[i+1][mi:mf]))
		# done! now make the plot
		newplot = FigPanel(self.tabs)
		newplot.fname = 'stitched'
		newplot.title = 'Rigel TGS Stitched Spectrum'
		newplot.plot(wave_stitch, intensity_stitch)
		curnames = [self.tabs.GetPageText(i) for i in range(self.tabs.GetPageCount())]
		n = 1
		while "Stitched - %d" % (n) in curnames:
			n+=1
		self.tabs.AddPage(newplot, "Stitched - %d" % (n), select=True)
		self.stitch.Show(False)

	def onExit(self, event):
		while self.tabs.GetSelection()>=0:
			if self.onClose(None) == wx.ID_CANCEL: return
		self.Destroy()
		
	def onZoom(self, event, mult):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer): cur.zoom(mult)
		
	def onZoomFit(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer): cur.zoom_fit()

	def onSetPosition(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, ImgViewer): return
		dialog = wx.TextEntryDialog(None, "Enter x,y coordinates (pixels)", "Set Position", "%d,%d" % (cur.specZero.x, cur.specZero.y))
		if dialog.ShowModal() == wx.ID_OK:
			pt=dialog.GetValue().split(',')
			cur.specZero=wx.Point(int(pt[0]),int(pt[1]))
			cur.Refresh()
		dialog.Destroy()
	
	def onWidth(self, event, inc):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer):
			cur.sel_width(inc)
				
	def onTilt(self, event, inc):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer):
			cur.sel_tilt(inc)
		
	def onShowFITS(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer):
			if cur.IsSplit():
				cur.Unsplit(cur.info)
			else:
				cur.SplitVertically(cur.info, cur.scroll, 200)

	def onPlot(self, event):
		self.readCheckItems()
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if isinstance(cur, ImgViewer):
			name=self.tabs.GetPageText(self.tabs.GetSelection())
			(wave, intensity_wt, spec) = cur.extract()
			newplot = FigPanel(self.tabs)
			curnames = [self.tabs.GetPageText(i) for i in range(self.tabs.GetPageCount())]
			n = 1
			while "%s - Plot %d" % (name, n) in curnames:
				n+=1
			self.tabs.AddPage(newplot, "%s - Plot %d" % (name, n), select=True)
			newplot.spec = spec
			newplot.fname = cur.fname
			newplot.title = 'Rigel TGS Spectrum from %s' % (cur.fname)
			newplot.plot(wave, intensity_wt)
		elif isinstance(cur, FigPanel):
			cur.plot()
			
	def onSetTitle(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, FigPanel): return
		if not hasattr(cur, 'title'): cur.title=''
		dialog = wx.TextEntryDialog(None, "Enter plot title", "Set title", cur.title)
		if dialog.ShowModal() == wx.ID_OK:
			cur.title=dialog.GetValue()
			self.readCheckItems()
			cur.plot()
		dialog.Destroy()

	def onAdjImage(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, ImgViewer): return
		dialog = AdjustImageDialog(parent=None, imgv=cur)
		if dialog.ShowModal() == wx.ID_OK:
			cur.imin, cur.imax = dialog.get_minmax()
			cur.adjustImage()
		else:
			cur.adjustImage()
		dialog.Destroy()

	def onSetRange(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, FigPanel): return
		dialog = SetRangeDialog(parent=None, wmin=cur.wmin, wmax=cur.wmax)
		if dialog.ShowModal() == wx.ID_OK:
			cur.wmin, cur.wmax = dialog.get_range()
			cur.plot()
		dialog.Destroy()
		
	def onAdjust(self, event):
		n = self.tabs.GetSelection()
		if n<0: return
		cur = self.tabs.GetPage(n)
		if not isinstance(cur, FigPanel):
			return
		cur.toolbar.configure_subplots(None)
			
	def onSetAverage(self, event):
		global MEDAVG_WIDTH
		val = "%s" % MEDAVG_WIDTH
		dialog = wx.TextEntryDialog(None, "Enter median average width (in pixels)", "Set Averaging", val)
		if dialog.ShowModal() == wx.ID_OK:
			a = int(dialog.GetValue())
			if np.mod(a,2) == 0: a+=1
			if a<1: a=1
			MEDAVG_WIDTH=a
		dialog.Destroy()
	
	def onSetRedshift(self, event):
		global REDSHIFT_Z
		val = "%s" % REDSHIFT_Z
		dialog = wx.TextEntryDialog(None, "Enter redshift value (z)", "Set redshift", val)
		if dialog.ShowModal() == wx.ID_OK:
			REDSHIFT_Z=float(dialog.GetValue())
		dialog.Destroy()

	def onHelp(self, event):
		self.help.Show(True)
		self.help.Raise()

	def onAbout(self, event):
		description = """Extracts spectra from FITS images,
applies site-specific corrections, and
produces plots and calibrated spectral data.

Provides options to median-smooth, apply
redshift, and overlay common spectral
lines. Code can be easily modified to
apply corrections for different sites.
"""
		info = wx.AboutDialogInfo()
		
		info.SetName('TGS Analyzer')
		info.SetVersion(VERSION_STR)
		info.SetDescription(description)
		info.SetCopyright('(C) 2013 University of Iowa Physics and Astronomy')
		info.SetWebSite('http://astro.physics.uiowa.edu/rigel')
		info.AddDeveloper('Robert Mutel (robert-mutel@uiowa.edu)')
		info.AddDeveloper('Bill Peterson (bill.m.peterson@gmail.com)')
		
		wx.AboutBox(info)
		return
	
	def onBroadcast(self, event):
		msg = event.getMsg()
		self.SetStatusText(msg)

class StitchDialog(wx.Dialog):
	def __init__(self, *args, **kw):
		super(StitchDialog, self).__init__(style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE, *args, **kw) 			
		self.SetSize((550, 450))
		self.SetTitle("TGS Stitch")

		ctrl = wx.FlexGridSizer(2, 3, vgap=8, hgap=8)
		ctrl2 = wx.BoxSizer(wx.VERTICAL)
		hbox = wx.BoxSizer(wx.VERTICAL)
		
		self.plots = wx.ListBox(self, style=wx.LB_EXTENDED)
		self.parts = wx.ListBox(self, style=wx.LB_EXTENDED)
		addButton = wx.Button(self, label='Add >>')
		remButton = wx.Button(self, label='<< Remove')
		rb1 = self.useSplits = wx.RadioButton(self, label='Split spectra equally by wavelength', style=wx.RB_GROUP)
		rb2 = self.pickSplits = wx.RadioButton(self, label='Specify wavelengths at which to split (nm):')
		self.splits = wx.TextCtrl(self, size=wx.Size(350,22))
		self.splits.Enable(False)
		
		ctrl.Add(wx.StaticText(self, label='Available Plots:'), flag=wx.ALL|wx.ALIGN_LEFT, border=5)
		ctrl.Add(wx.StaticText(self, label=''))
		ctrl.Add(wx.StaticText(self, label='Plots to Stitch:'), flag=wx.ALL|wx.ALIGN_LEFT, border=5)
		ctrl.Add(self.plots, border=5, flag=wx.ALL|wx.EXPAND)
		ctrl.Add(ctrl2, flag=wx.ALIGN_CENTER)
		ctrl.Add(self.parts, border=5, flag=wx.ALL|wx.EXPAND)
		ctrl.AddGrowableCol(0,1)
		ctrl.AddGrowableCol(2,1)
		ctrl.AddGrowableRow(1,1)

		ctrl2.Add(addButton, border=5, flag=wx.ALL|wx.ALIGN_CENTER)
		ctrl2.Add(remButton, border=5, flag=wx.ALL|wx.ALIGN_CENTER)
		
		hbox.Add(self.useSplits, border=5, flag=wx.ALIGN_LEFT|wx.ALL)
		hbox.Add(self.pickSplits, border=5, flag=wx.ALIGN_LEFT|wx.ALL)
				
		vbox = wx.BoxSizer(wx.VERTICAL)
		vbox.Add(ctrl, flag=wx.ALL|wx.EXPAND, border=5, proportion=1)
		vbox.Add(hbox, flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		vbox.Add(self.splits, flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		vbox.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL), flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		self.SetSizer(vbox)
		
		addButton.Bind(wx.EVT_BUTTON, self.add_sel, addButton)
		remButton.Bind(wx.EVT_BUTTON, self.rem_sel, remButton)
		self.useSplits.Bind(wx.EVT_RADIOBUTTON, self.toggle_entry, rb1)
		self.pickSplits.Bind(wx.EVT_RADIOBUTTON, self.toggle_entry, rb2)
		self.Centre()
		
	def repopulate(self, name_dict):
		# takes figpanel object ids and tab names from main window
		self.name_dict = name_dict
		self.plots_dict = {}
		self.parts_dict = {}
		self.plots.Clear()
		self.parts.Clear()
		for id in self.name_dict.keys():
			pos = self.plots.GetCount()
			self.plots.Insert(name_dict[id], pos)
			self.plots_dict[pos] = id
		
	def add_sel(self, event):
		for k in self.plots.GetSelections():
			id = self.plots_dict[k]
			pos = self.parts.GetCount()
			self.parts.Insert(self.name_dict[id], pos)
			self.parts_dict[pos] = id
			del self.plots_dict[k]
		plots = self.plots_dict.values()
		self.plots.Clear()
		self.plots_dict = {}
		for id in plots:
			pos = self.plots.GetCount()
			self.plots.Insert(self.name_dict[id], pos)
			self.plots_dict[pos] = id
	
	def rem_sel(self, event):
		for k in self.parts.GetSelections():
			id = self.parts_dict[k]
			pos = self.plots.GetCount()
			self.plots.Insert(self.name_dict[id], pos)
			self.plots_dict[pos] = id
			del self.parts_dict[k]
		parts = self.parts_dict.values()
		self.parts.Clear()
		self.parts_dict = {}
		for id in parts:
			pos = self.parts.GetCount()
			self.parts.Insert(self.name_dict[id], pos)
			self.parts_dict[pos] = id
	
	def toggle_entry(self, event):
		if self.pickSplits.GetValue():
			self.splits.Enable(True)
		else:
			self.splits.Enable(False)
		
	def get_stitched(self):
		ids = self.parts_dict.values()
		if self.pickSplits.GetValue():
			wsplits = [float(x) for x in self.splits.GetValue().split(',')]
		else:
			wsplits = [DEFAULT_WMIN + (DEFAULT_WMAX - DEFAULT_WMIN)*(x+1.)/len(ids) for x in range(len(ids)-1)]
		return ids, wsplits
		
class SetRangeDialog(wx.Dialog):
	def __init__(self, wmin, wmax, *args, **kw):
		super(SetRangeDialog, self).__init__(*args, **kw) 	
		self.SetTitle("Set Plot Range")
		self.SetSize((350,150))
		mintext = "%s" % wmin
		maxtext = "%s" % wmax

		self.min = wx.TextCtrl(self, style=wx.TE_RIGHT,size=wx.Size(150,22), value=mintext)
		self.max = wx.TextCtrl(self, style=wx.TE_RIGHT,size=wx.Size(150,22), value=maxtext)		
		ctrl = wx.FlexGridSizer(2, 2, vgap=8, hgap=8)
		vbox = wx.BoxSizer(wx.VERTICAL)

		szf=wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL
		ctrl.Add(wx.StaticText(self, label='Minimum (nm):'), flag=szf)
		ctrl.Add(self.min, flag=szf)
		ctrl.Add(wx.StaticText(self, label='Maximum (nm):'), flag=szf)
		ctrl.Add(self.max, flag=szf)
		ctrl.AddGrowableCol(0,0)

		vbox.Add(ctrl, border=5, flag=wx.ALL|wx.EXPAND, proportion=1)
		vbox.Add(self.CreateButtonSizer(wx.OK|wx.CANCEL), flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		self.SetSizer(vbox)
		self.Centre()
		
	def get_range(self):
		wmin = float(self.min.GetValue())
		wmax = float(self.max.GetValue())
		return wmin, wmax
		
class AdjustImageDialog(wx.Dialog):
	def __init__(self, imgv, *args, **kw):
		super(AdjustImageDialog, self).__init__(style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE, *args, **kw)
		self.SetTitle("Set Image Levels")
		self.SetSize((450,250))
		amax=np.max(imgv.data)
		amin=np.min(imgv.data)
		self.imgv = imgv
		
		self.min = wx.Slider(self, style=wx.SL_LABELS, minValue=amin, maxValue=amax, value=imgv.imin)
		self.max = wx.Slider(self, style=wx.SL_LABELS, minValue=amin, maxValue=amax, value=imgv.imax)
		ctrl = wx.FlexGridSizer(2, 2, vgap=8, hgap=8)
		vbox = wx.BoxSizer(wx.VERTICAL)

		ctrl.Add(wx.StaticText(self, label='Minimum Intensity:'), flag=wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
		ctrl.Add(self.min, flag=wx.EXPAND)
		ctrl.Add(wx.StaticText(self, label='Maximum Intensity:'), flag=wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
		ctrl.Add(self.max, flag=wx.EXPAND)
		ctrl.AddGrowableCol(1,1)

		vbox.Add(ctrl, border=5, flag=wx.ALL|wx.EXPAND, proportion=1)
		bs=self.CreateButtonSizer(wx.OK|wx.CANCEL)
		apply=wx.Button(self,label='Preview')
		bs.Add(apply)
		vbox.Add(bs, flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		self.SetSizer(vbox)
		self.Bind(wx.EVT_BUTTON, self.onApply, apply)
		self.Centre()
		
	def onApply(self, event):
		min,max=self.get_minmax()
		self.imgv.adjustImage(min,max)
		
	def get_minmax(self):
		return (self.min.GetValue(), self.max.GetValue())
		
class HelpDialog(wx.Dialog):
	def __init__(self, *args, **kw):
		super(HelpDialog, self).__init__(*args, **kw) 	
		self.SetTitle("Help")
		self.SetSize((400,500))
		
		help = html.HtmlWindow(self)
		help.LoadPage('tgsa_help.html')
		
		vbox = wx.BoxSizer(wx.VERTICAL)
		vbox.Add(help, border=5, flag=wx.ALL|wx.EXPAND, proportion=1)
		vbox.Add(self.CreateButtonSizer(wx.OK), flag=wx.ALL|wx.ALIGN_CENTER, border=5)
		self.SetSizer(vbox)
		self.Centre()

### Nuts and Bolts

def f_lambda(x,x_star):
    '''
    Input: displacement from the zeroth order image (x, pixels) and star at x_star
    Returns: wavelength in nanometers
    '''
    global L, d_g, FL, x_ctr
    dx = (x_ctr - x_star) * pixel
    x *= pixel
    a0 = 0.004; a1 = 2.2; a2 = -40 # Determined by best-fit to actual spectra
    x2 = x + a1*(x**2) + a2 * (x**3)
    lam = (d_g/nm) * ( sin( atan(x2/L) ) + a0 )
    return lam

def grating_sens(lam):
	'''
	Models eff. (0.0-1.0) of Edmunds 600 lpmm grating using published efficiency curve
	Input wavelength, nm
	'''
	sigma = 120; a= 110; lam0 = 250
	t = np.abs(float(lam) -lam0) / sigma
	eff =  ( (lam - lam0)/a)**2 * np.exp(-t)
	eff /= 100.
	return eff

def ccd_sens(lam):
	'''
	Models FLI CCD QE (0.0 - 1.0) vs wavelength (nm) - approximate fit to QE curve
	'''
	sigma = 130; a= 125; lam0 = 260
	t = np.abs(float(lam) -lam0) / sigma
	eff =  ( (lam - lam0)/a)**2 * np.exp(-t)
	eff /= 100.
	return eff

def plt_Balmer(axis,z):
	z1 = 1 + z
	(y1,y2) = axis.get_ylim()
	axis.vlines(656.3*z1,y1,y2,linestyle='solid', linewidth='1',color ='r', label=r'H${\alpha}$ 656.3')
	axis.vlines(486.1*z1,y1,y2,linestyle='solid', linewidth='1',color = 'g',label=r'H${\beta}$ 486.1')
	axis.vlines(434.3*z1,y1,y2,linestyle='solid', linewidth='1',color = 'b',label=r'H${\gamma}$ 434.3')
	axis.vlines(410.2*z1,y1,y2,linestyle='solid', linewidth='1',color = 'm',label=r'H${\delta}$ 410.2')
	axis.vlines(397.0*z1,y1,y2,linestyle='solid', linewidth='1',color = 'm',label=r'H${\epsilon}$ 397.0')

def plt_Helium(axis,z):
	z1 = 1 + z
	(y1,y2) = axis.get_ylim()
	axis.vlines(501.5*z1,y1,y2,linestyle='dashdot', linewidth='1',color = 'r',label=r'HeI 501.5')
	axis.vlines(587.6*z1,y1,y2,linestyle='dashdot', linewidth='1',color = 'g',label=r'HeI 587.6')
	axis.vlines(667.8*z1,y1,y2,linestyle='dashdot', linewidth='1',color = 'b',label=r'HeI 667.8')
	axis.vlines(706.5*z1,y1,y2,linestyle='dashdot', linewidth='1',color = 'm',label=r'HeI 706.5')

def plt_Metals(axis,z):
	z1 = 1 + z
	(y1,y2) = axis.get_ylim()
	axis.vlines(715.0*z1,y1,y2,linestyle='dotted', linewidth='2',color ='k',label='TiO 715.0')
	axis.vlines(410.0*z1,y1,y2,linestyle='dotted', linewidth='2',color ='b',label='HeII 410.0')
	axis.vlines(464.0*z1,y1,y2,linestyle='dotted', linewidth='2',color ='g',label='NIII 464.0 ')
	axis.vlines(465.0*z1,y1,y2,linestyle='dotted', linewidth='2',color ='r',label='CIV 465.0')
	axis.vlines(468.6*z1,y1,y2,linestyle='dotted', linewidth='2',color ='b',label='HeII 468.6')
	axis.vlines(541.1*z1,y1,y2,linestyle='dotted', linewidth='2',color ='b',label='HeII 541.1')
	axis.vlines(569.6*z1,y1,y2,linestyle='dotted', linewidth='2',color ='r',label='CIII 569.6')
	axis.vlines(580.5*z1,y1,y2,linestyle='dotted', linewidth='2',color ='r',label='CIV 580.5')
    
def plt_Telluric(axis,z):
	z1 = 1 + z
	(y1,y2) = axis.get_ylim()
	axis.vlines(759.0,y1,y2,linestyle='solid', linewidth='1',color ='m',label='O$_{2}$ (telluric)')

def mkconfigfile():
	newcfg = ConfigParser.ConfigParser()
	cfgfile = open('tgsa_cfg.ini', 'w')
	newcfg.add_section('Defaults')
	newcfg.set('Defaults','WindowWidth',DEFAULT_APP_WIDTH)
	newcfg.set('Defaults','WindowHeight',DEFAULT_APP_HEIGHT)
	newcfg.set('Defaults','SelWidth',DEFAULT_WIDTH)
	newcfg.set('Defaults','WavelengthMin',DEFAULT_WMIN)
	newcfg.set('Defaults','WavelengthMax',DEFAULT_WMAX)
	newcfg.set('Defaults','RedshiftZ',REDSHIFT_Z)
	newcfg.set('Defaults','PlotBalmer',PLOT_BALMER)
	newcfg.set('Defaults','PlotHelium',PLOT_HELIUM)
	newcfg.set('Defaults','PlotMetallic',PLOT_METALLIC)
	newcfg.set('Defaults','PlotTelluric',PLOT_TELLURIC)
	
	newcfg.add_section('Telescope')
	newcfg.set('Telescope','FocalRatio',f_ratio)
	newcfg.set('Telescope','Diameter_cm',Diam/cm)
	newcfg.set('Telescope','SensorDistance_mm',L/mm)
	newcfg.set('Telescope','GratingLines_mm',lpmm*mm)
	newcfg.set('Telescope','NPixels',npixel)
	newcfg.set('Telescope','PixelSize_um',pixel/micron)
	newcfg.set('Telescope','PixelStart',PIXEL_START)
	newcfg.set('Telescope','PixelEnd',PIXEL_END)
	newcfg.set('Telescope','PixelScale',PIXEL_SCALE)

	newcfg.add_section('Advanced')
	newcfg.set('Advanced','StitchWidth',STITCH_WIDTH)
	newcfg.set('Advanced','MedAvgWidth',MEDAVG_WIDTH)
	newcfg.set('Advanced','ZoomMax',ZOOM_MAX)
	newcfg.set('Advanced','ZoomMin',ZOOM_MIN)
	newcfg.set('Advanced','TiltInc_deg',TILT_INC/deg)
	newcfg.write(cfgfile)
	cfgfile.close()

def envpath(envvar, *paths):
	a = os.getenv(envvar)
	if a == None: return ''
	return os.path.join(a,*paths)
	
### MAIN

# hunt for configuration files in various locations
cfgfiles = [ os.path.join( os.getcwd(), 'tgsa_cfg.ini'),
envpath('USERPROFILE', 'AppData', 'Local', 'TGS Analyzer', 'tgsa_cfg.ini'),
envpath('HOME', 'Library', 'tgsa_cfg.ini'),
envpath('HOME', '.tgsarc')
]
cfg = ConfigParser.ConfigParser()
if cfg.read(cfgfiles) == []:
	mkconfigfile()
else:
	try:
		DEFAULT_APP_WIDTH=cfg.getint('Defaults','WindowWidth')
		DEFAULT_APP_HEIGHT=cfg.getint('Defaults','WindowHeight')
		DEFAULT_WIDTH=cfg.getint('Defaults','SelWidth')
		DEFAULT_WMIN=cfg.getfloat('Defaults','WavelengthMin')
		DEFAULT_WMAX=cfg.getfloat('Defaults','WavelengthMax')
		REDSHIFT_Z=cfg.getfloat('Defaults','RedshiftZ')
		PLOT_BALMER=cfg.getboolean('Defaults','PlotBalmer')
		PLOT_HELIUM=cfg.getboolean('Defaults','PlotHelium')
		PLOT_METALLIC=cfg.getboolean('Defaults','PlotMetallic')
		PLOT_TELLURIC=cfg.getboolean('Defaults','PlotTelluric')
		
		f_ratio=cfg.getfloat('Telescope','FocalRatio')
		Diam=cfg.getfloat('Telescope','Diameter_cm')*cm
		L=cfg.getfloat('Telescope','SensorDistance_mm')*mm
		lpmm=cfg.getfloat('Telescope','GratingLines_mm')/mm
		npixel=cfg.getint('Telescope','NPixels')
		pixel=cfg.getfloat('Telescope','PixelSize_um')*micron
		PIXEL_START=cfg.getint('Telescope','PixelStart')
		PIXEL_END=cfg.getint('Telescope','PixelEnd')
		PIXEL_SCALE=cfg.getfloat('Telescope','PixelScale')

		STITCH_WIDTH=cfg.getint('Advanced','StitchWidth')
		MEDAVG_WIDTH=cfg.getint('Advanced','MedAvgWidth')
		ZOOM_MAX=cfg.getfloat('Advanced','ZoomMax')
		ZOOM_MIN=cfg.getfloat('Advanced','ZoomMin')
		TILT_INC=cfg.getfloat('Advanced','TiltInc_deg')*deg
	except ConfigParser.NoOptionError:
		# make a new config file if one doesn't exist
		# or if a token is missing
		mkconfigfile()

app = wx.App(redirect=False)
frame = MainWindow()

frame.Show()
app.MainLoop()
