# Use PyGMT to make map of S40RTS globally
import pygmt


gproj = "Ks20c"
fig = pygmt.Figure()
fig.basemap(region='g', projection=gproj, frame="afg")
fig.grdimage(grid='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS_2800km.grd',
              cmap='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS.cpt')
fig.colorbar(frame=['a0.5', 'x+l"dVs (%)"' ], cmap='/Users/ja17375/DiscrePy/Data/S40RTS/S40RTS.cpt')
fig.coast(shorelines=True)
fig.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter01/Figs/Global_S40RTS.eps',
            crop=True, show=True)