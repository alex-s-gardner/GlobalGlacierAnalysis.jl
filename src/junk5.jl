err = 700
VX = 1.0
VY = 1.0
V = VX.^2
sqrt((err * VX / V)^2 + (err * VY / V)^2)