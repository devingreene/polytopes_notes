load('triangulate.sage')

P = [[0, 0, 0, 1, 1, 1, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 1, 1, 1],
     [0, 1, 0, 0, 1, 0, 0, 1, 0],
     [0, 0, 1, 0, 0, 1, 0, 0, 1]]
    
P = matrix(P)

N = 100 if not sys.argv[1:] else int(sys.argv[1])
import math
progress_interval = 10**int(math.log10(N//10))
triangulations = []
degeneracies = 0
for i in range(1,N+1):
    if i % progress_interval == 0:
        print(f"...{i}",end="",file=sys.stderr,flush=True)
    f = random_matrix(ZZ,1,P.ncols(),x = 5, y = 100)
    try:
        tri = find_triangulation(P,f)
    except Degeneracy as e:
        degeneracies += 1
        continue

    if tri in triangulations:
        continue
    triangulations += [tri]
else:
    print()

print(f"Number of triagulations found: {len(triangulations)}")
print(f"Number of degeneracies: {degeneracies}")
