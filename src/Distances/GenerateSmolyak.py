import numpy as np
import muq.Approximation as ma

# Define the 1d family of quadrature rules
orthoPoly = ma.Legendre()
scalarRule = ma.GaussQuadrature(orthoPoly)

# Define the 2d smolyak quadrature rule
smolyQuad = ma.SmolyakQuadrature(2, scalarRule)

for order in range(0,5):
  smolyQuad.Compute(order)
  pts = smolyQuad.Points()
  wts = smolyQuad.Weights()

  out = 'Eigen::Matrix2Xd pts_{:02d}(2,{});\n'.format(order,pts.shape[1])
  out += 'Eigen::VectorXd wts_{:02d}({});\n'.format(order,pts.shape[1])
  out += 'pts_{:02d} << '.format(order)

  for i in range(0, pts.shape[1]):
      out += '{:.15f}, '.format(0.5*(pts[0,i]+1))

  out += '\n          {:.15f}'.format(0.5*(pts[1,0]+1))
  for i in range(1, pts.shape[1]):
      out += ', {:.15f}'.format(0.5*(pts[1,i]+1))

  out += ';\nwts_{:02d} << {:.15f}'.format(order, 0.25*wts[0])
  for i in range(1,pts.shape[1]):
      out += ', {:0.15f}'.format(0.25*wts[i])
  out += ';'
  print(out)
  print('\n\n')
#  print(pts,wts)
