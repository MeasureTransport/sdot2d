import numpy as np
import muq.Approximation as ma

# Define the 1d family of quadrature rules
orthoPoly = ma.Legendre()
scalarRule = ma.GaussQuadrature(orthoPoly)

for order in range(0,7):
  scalarRule.Compute(order)
  pts = scalarRule.Points()
  wts = scalarRule.Weights()

  out = 'Eigen::VectorXd pts_{:02d}({});\n'.format(order,pts.shape[1])
  out += 'Eigen::VectorXd wts_{:02d}({});\n'.format(order,pts.shape[1])
  out += 'pts_{:02d} << '.format(order)

  out += '{:.15f}'.format(0.5*(pts[0,0]+1))
  for i in range(1, pts.shape[1]):
      out += ', {:.15f}'.format(0.5*(pts[0,i]+1))

  out += ';\nwts_{:02d} << {:.15f}'.format(order, 0.5*wts[0])
  for i in range(1,pts.shape[1]):
      out += ', {:0.15f}'.format(0.5*wts[i])
  out += ';\nreturn std::make_pair(pts_{:02d},wts_{:02d});'.format(order,order)
  print(out)
  print('\n\n')
#  print(pts,wts)
