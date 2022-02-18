from diagrams import* 

class ClosedMesh:
    def __init__(self, vertices = []):
        self.vertices = vertices

    def intersectLine(self, o, d, i):
        # TODO: Parallele Richtungsvektoren (singulare Matrix, determinante = 0)
        A, B, _ = self.getSurfaceAndNormal(i)

        beta, t = solve2x2(A - B, d, A - o)
        if beta and t:
            if beta >= 0 and beta <= 1 and t > 0:
                return t, i
        return None

    def intersect(self, ct, o, d):
        values = minIntersect([ self.intersectLine(o, d, i) for i in range(len(self.vertices)) ])
        if(type(values) == tuple):
            t, i = values
            _, _, n = self.getSurfaceAndNormal(i)
            return t, n
        return None

    def getSurfaceAndNormal(self, i):
        A, B = self.getSurface(i)
        return A, B, normalize(getOrthoVec(np.subtract(B, A)))

    def getSurface(self, i):
        # A = vertex(i), B = vertex((i+1)%N)
        A = self.vertices[i]
        B = self.vertices[( i + 1) % len(self.vertices) ]
        return np.array(A), np.array(B)
    
    def draw(self, ct, pos, drawNormals = True, drawExtends=False):
        # für alle subaschnitte linien zeichnen
        # für alle subabschnitte normale zeichnen.
        polygon(ct, pos, self.vertices, (Hue.VIOLET,.2,1.))
        if drawNormals:
            for i in range(len(self.vertices)):
                A, B, n = self.getSurfaceAndNormal(i)
                mid = A + (B - A)/2
                vector(ct, mid, n, 0.2)
        
        if drawExtends:
            a, b, c, d = self.getExtends()
            ct.rectangle(a, b, c - a, d - b)
            ct.stroke()

    def getExtends(self):
        xmin = min(x for (x, y) in self.vertices)
        ymin = min(y for (x, y) in self.vertices)
        xmax = max(x for (x, y) in self.vertices)
        ymax = max(y for (x, y) in self.vertices)
        return (xmin, ymin, xmax, ymax)