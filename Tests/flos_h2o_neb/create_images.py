import sisl

init = sisl.Geometry.read('h2o_neb.fdf')
fixed = init.sub([0,1,2])
sub = init.sub([3,4,5])

# Create images
o = sub.xyz[0, :]
for i, ang in enumerate([0, 30, 60, 90, 120, 150, 180]):
    new = fixed.add(sub.translate(-o).rotatec(ang, 'xyz').translate(o))
    new.write('image_{}.xyz'.format(i))
