# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1125, 801]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [120.0, 25.0, 170.0]
renderView1.StereoType = 0
renderView1.CameraPosition = [-200.54799428779967, -608.4265779325958, 565.9157604526779]
renderView1.CameraFocalPoint = [150.751294739542, 25.308187189361615, 200.0617652356396]
renderView1.CameraViewUp = [-0.9007503588501264, 0.3541685364912382, -0.25142282869903215]
renderView1.CameraParallelScale = 210.08711746916933
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.Background2 = [1.0, 1.0, 1.0]
renderView1.UseGradientBackground = 1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = 'Z (m)'
renderView1.AxesGrid.YTitle = 'Y (m)'
renderView1.AxesGrid.ZTitle = 'X (m)'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleFontSize = 18
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleFontSize = 18
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleFontSize = 18
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelFontSize = 18
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelFontSize = 18
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelFontSize = 18

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Image Data Reader'
imageEz = XMLImageDataReader(FileName=['/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.11.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.21.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.31.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.41.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.51.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.61.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.71.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.81.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.91.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.101.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.111.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.121.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.131.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.141.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.151.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.161.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.171.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.181.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.191.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.201.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.211.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.221.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.231.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.241.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.251.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.261.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.271.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.281.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.291.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.301.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.311.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.321.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.331.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.341.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.351.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.361.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.371.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.381.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.391.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.401.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.411.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.421.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.431.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.441.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.451.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.461.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.471.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.481.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.491.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.501.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.511.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.521.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.531.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.541.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.551.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.561.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.571.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.581.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.591.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.601.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.611.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.621.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.631.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.641.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.651.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.661.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.671.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.681.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.691.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.701.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.711.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.721.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.731.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.741.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.751.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.761.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.771.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.781.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.791.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.801.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.811.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.821.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.831.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.841.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.851.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.861.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.871.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.881.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.891.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.901.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.911.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.921.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.931.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.941.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.951.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.961.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.971.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.981.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.991.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1001.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1011.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1021.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1031.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1041.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1051.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1061.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1071.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1081.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1091.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1101.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1111.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1121.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1131.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1141.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1151.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1161.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1171.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1181.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1191.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1201.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1211.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1221.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1231.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1241.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1251.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1261.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1271.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1281.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1291.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1301.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1311.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1321.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1331.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1341.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1351.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1361.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1371.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1381.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1391.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1401.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1411.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1421.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1431.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1441.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1451.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1461.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1471.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1481.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1491.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1501.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1511.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1521.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1531.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1541.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1551.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1561.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1571.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1581.vti', '/home/bernsen/SeidarT/EXAMPLES/dipping_bed/imageEz.1591.vti'])
imageEz.CellArrayStatus = ['Displacement']

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'Displacement'
displacementLUT = GetColorTransferFunction('Displacement')
displacementLUT.RGBPoints = [-8.071549382293597e-05, 0.231373, 0.298039, 0.752941, -9.04645276023075e-07, 0.865003, 0.865003, 0.865003, 7.890620327088982e-05, 0.705882, 0.0156863, 0.14902]
displacementLUT.ColorSpace = 'RGB'
displacementLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Displacement'
displacementPWF = GetOpacityTransferFunction('Displacement')
displacementPWF.Points = [-8.071549382293597e-05, 1.0, 0.5, 0.0, 9.926350230671233e-08, 0.0, 0.5, 0.0, 7.890620327088982e-05, 1.0, 0.5, 0.0]
displacementPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from imageEz
imageEzDisplay = Show(imageEz, renderView1)
# trace defaults for the display properties.
imageEzDisplay.Representation = 'Volume'
imageEzDisplay.ColorArrayName = ['CELLS', 'Displacement']
imageEzDisplay.LookupTable = displacementLUT
imageEzDisplay.OSPRayScaleArray = 'Displacement'
imageEzDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
imageEzDisplay.SelectOrientationVectors = 'None'
imageEzDisplay.ScaleFactor = 34.0
imageEzDisplay.SelectScaleArray = 'None'
imageEzDisplay.GlyphType = 'Arrow'
imageEzDisplay.GlyphTableIndexArray = 'None'
imageEzDisplay.DataAxesGrid = 'GridAxesRepresentation'
imageEzDisplay.PolarAxes = 'PolarAxesRepresentation'
imageEzDisplay.ScalarOpacityUnitDistance = 2.6232065800759536
imageEzDisplay.ScalarOpacityFunction = displacementPWF
imageEzDisplay.Slice = 170

# show color legend
imageEzDisplay.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for displacementLUT in view renderView1
displacementLUTColorBar = GetScalarBar(displacementLUT, renderView1)
displacementLUTColorBar.AutoOrient = 0
displacementLUTColorBar.Orientation = 'Horizontal'
displacementLUTColorBar.WindowLocation = 'AnyLocation'
displacementLUTColorBar.Position = [0.3233999159840363, 0.09248439450686707]
displacementLUTColorBar.Title = 'Displacement'
displacementLUTColorBar.ComponentTitle = ''
displacementLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
displacementLUTColorBar.TitleFontSize = 18
displacementLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
displacementLUTColorBar.LabelFontSize = 18
displacementLUTColorBar.ScalarBarThickness = 20
displacementLUTColorBar.ScalarBarLength = 0.33000000000000074

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(imageEz)
# ----------------------------------------------------------------