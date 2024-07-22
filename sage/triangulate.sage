HEIGHT_MAX = 1000000

class Degeneracy(ValueError):
    pass

def all_nonneg(W):
    return all(x >= 0 for x in W.list())

def all_neg(W):
    return all(x < 0 for x in W.list())

def all_ltez(W):
    return all(x <= 0 for x in W.list())

def make_W_matrix(P,f):
    f = f.list()
    assert len(f) == P.ncols()
    P = matrix([1]*P.ncols()).stack(P)
    P = matrix(P.nrows(),1,[0]*P.nrows()).change_ring(QQ).augment(P)
    P = matrix([-1] + f).stack(P)
    W = P.right_kernel_matrix()
    if W[:,0] != matrix([1] + [0]*(W.nrows() - 1)).T:
        raise Degeneracy(f"Degeneracy detected: kernel of\n"
                         f"{P}\nwith\n{f}\nis\n{W}")
    return W

def iterator_n_k(n,k):
    assert 0 <= k <= n
    if n == 0:
        yield []
        return

    b_s = [1]*k + [0]*(n - k)

    yield b_s.copy()
    while True:
        i = n - 1
        while i >= 0 and b_s[i] == 1:
            i -= 1
        if i < 0:
            return
        while i >= 0 and b_s[i] == 0:
            i -= 1
        if i < 0:
            return
        b_s[i],b_s[i+1] = b_s[i+1],b_s[i]
        b_s[i+2:] = b_s[n-1:i+1:-1]
        yield b_s.copy()

def reduce_column(V, r, c, height_max = HEIGHT_MAX):
    assert V.base_ring() == QQ
    assert 0 <= r <= V.nrows() - 2
    assert 0 <= c <= V.ncols() - 1
    assert V[r,c] != 0
    W = V[:-1,:]
    W[r,:] /= W[r,c]
    for i in range(W.nrows()):
        if i == r: continue
        W.add_multiple_of_row(i,r,-W[i,c])
    if W.height() > height_max:
        raise ValueError(f"maxheight exceeded:\n{V}")
    V[:-1,:] = W

def haswall(V):
    def alert():
        raise Degeneracy(f"Degeneracy detected\n{V}")
        
    assert V.nrows() >= 1
    pivots = [0]
    while True:
        if len(pivots) == V.nrows() - 1:
            break
        for j in range(1,V.ncols()):
            col = vector(V[:,j])
            if col[0] < 0 and col[-1] == 0:
                if col[1:-1] == 0:
                    alert()
                break
        else:
            return False
        for i in range(1,len(col) - 1):
            if i in pivots: continue
            if col[i] != 0:
                pivots += [i]
                reduce_column(V,i,j,HEIGHT_MAX)
                break
        else:
            alert()

    assert len(pivots) == V.nrows() - 1

    for j in range(1,V.ncols()):
        col = vector(V[:,j])
        if all_neg(col[:-1]) and col[-1] == 0:
            return True
    return False

def check_all_triangles(P,W,k):
    envelope = []
    P = P.stack(matrix([1]*P.ncols()))
    for cntpt in iterator_n_k(W.ncols() - 1,k + 1):
        # Pass over faces
        support = [i for i,bt in enumerate(cntpt) if bt]
        if rank(P[:,support]) < k + 1:
            continue
        V = W.stack(matrix([0] + cntpt))
        if haswall(V):
            envelope += [cntpt]
    return envelope

def find_triangulation(P,f):
    W = make_W_matrix(P,f)
    k = P.nrows()
    return check_all_triangles(P,W,k)
