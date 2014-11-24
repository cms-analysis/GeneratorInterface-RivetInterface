# BEGIN PLOT /CMS_SMP_13_004/d01-x01-y01
Title= $p_T$ Z distribution (GeV)
XLabel= $p_{T}^Z$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}p_{T}^Z$
LogX=0
LogY=0
#YMin= 0.000001 when log to see the data
YMin= -0.001
RatioPlotYMin= -7.0
RatioPlotYMax= 3.0
LegendXPos=0.37
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d02-x01-y01
Title= $p_T$ leading b jet distribution (GeV)
XLabel= $p_{T}^{b-lead}$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}p_{T}^{b-lead}$
LogX=0
LogY=1
RatioPlotYMax= 1.7
LegendXPos=0.37
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d03-x01-y01
Title= $p_T$ subleading b jet distribution (GeV)
XLabel= $p_{T}^{b-sublead}$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}p_{T}^{b-sublead}$
LogX=0
LogY=1
YMin= 0.0000001
RatioPlotYMax= 17.0
LegendXPos=0.37
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d04-x01-y01
Title= $p_T$ di-bjets distribution (GeV)
XLabel= $p_{T}^{bb}$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}p_{T}^{bb}$
LogX=0
LogY=1
RatioPlotYMin= 0.2
LegendXPos=0.35
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d05-x01-y01
Title= Mass di-bjets distribution (GeV)
XLabel= $M_{bb}$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}M_{bb}$
LogX=0
LogY=1
RatioPlotYMin= 0.0
LegendXPos=0.35
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d06-x01-y01
Title= $\Delta R_{b,b}$  distribution 
XLabel= $\Delta R_{b,b}$
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}\Delta R_{b,b}$
LogX=0
LogY=0
RatioPlotYMin= 0.35
RatioPlotYMax= 1.6
LegendXPos=0.01
#LegendYPos=0.01
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d07-x01-y01
Title= $\Delta\phi_{Z,bb}$  distribution 
XLabel= $\Delta\phi_{Z,bb}$
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}\Delta\phi_{Z,bb}$
LogX=0
LogY=0
YMin= -0.2
RatioPlotYMin= -8.2
RatioPlotYMax= 3.2
LegendXPos=0.05
# END PLOT

# BEGIN PLOT /CMS_SMP_13_004/d08-x01-y01
Title= Mass Zbb system distribution (GeV)
XLabel= $M_{Zbb}$ (GeV)
YLabel= $1/\sigma \; \text{d}\sigma/\text{d}M_{Zbb}$
LogX=0
LogY=1
RatioPlotYMin= 0.0
RatioPlotYMax= 1.8
LegendXPos=0.40
# END PLOT

# + any additional plot settings you might like, see make-plots documentation
# ... add more histograms as you need them ...
