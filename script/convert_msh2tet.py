import meshio
import sys
import os

input_file = sys.argv[1]
output_file = os.path.splitext(input_file)[0]+'.tet'

mesh = meshio.read(input_file)
p_size = len(mesh.points)
t_size = len(mesh.cells_dict['tetra'])

f = open(output_file, "w")

# write size
f.write("%d %d\n" % (p_size, t_size))

# write vertices
for p in mesh.points:
    # print(p)
    f.write("%f %f %f\n" % (p[0], p[1], p[2]))

# write tets
for t in mesh.cells_dict['tetra']:
    f.write("4 %d %d %d %d\n" % (t[0], t[1], t[2], t[3]))
    # break

f.close()

print('save .tet file to %s' % (output_file))
