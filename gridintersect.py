from closedmesh import ClosedMesh
from diagrams import *
import math
from colorsys import hsv_to_rgb
import sys

FLOATMAX = sys.float_info.max

# z-ordering for cairo
class DrawQueue:
    def __init__(self):
        self.drawCalls = [] 
    def addCall(self, z_order, call):
        if call == None:
            raise ValueError('call is null')
        self.drawCalls.append( (z_order, call) )
    def draw(self):
        for z, call in sorted(self.drawCalls, key=lambda x:x[0]):
            print(z)
            call()
    def reset(self):
        self.drawCalls = []
    
drawQueue = DrawQueue()

class Grid:
    def __init__(self, ct, gridres = (5, 5), gridpos = (1, 1), gridsize=(5, 5)):
        self.ct = ct
        self.gridres = gridres
        self.gridpos = gridpos
        self.gridsize = gridsize

        self.cells = [ [ [] for _ in range(gridres[1]) ] for _ in range(gridres[0]) ]

    @classmethod
    def setupGrid(GridCls, ct , meshes, mink = 5):
        def getSceneExtends(meshes, epsilon = 0.001):
            allExtends = [mesh.getExtends() for mesh in meshes ]
            xmin = min([ xmin for xmin, _, _, _ in allExtends ])
            ymin = min([ ymin for _, ymin, _, _ in allExtends ])
            xmax = max([ xmax for _, _, xmax, _ in allExtends ])
            ymax = max([ ymax for _, _, _, ymax in allExtends ])
            return xmin - epsilon, ymin - epsilon, xmax + epsilon, ymax + epsilon
        
        if len(meshes) < 1:
            return
        k = max(math.floor(np.sqrt(len(meshes))), mink)

        xmin, ymin, xmax, ymax = getSceneExtends(meshes)
        grid = Grid(ct, gridres=(k, k), gridpos=(xmin, ymin),gridsize=(xmax-xmin, ymax-ymin))
        for mesh in meshes:
            grid.placeIntoGrid(mesh)

        return grid

    def getAllDims(self):
        left, bottom, xres, yres, w, h = (*self.gridpos,  *self.gridres, *self.gridsize)
        right, top = left+w, bottom+h
        return left, bottom, right, top, xres, yres, w, h

    def getGridCellIndicesAt(self, x, y):
        left, bottom, _, _, xres, yres, w, h = self.getAllDims()
        ix = min(max(math.floor((x - left) * xres / w), 0), xres-1)
        iy = min(max(math.floor((y - bottom) * yres / h), 0), yres-1)
        return ix, iy

    def placeIntoGrid(self, mesh):
        mesh_xmin, mesh_ymin, mesh_xmax, mesh_ymax = mesh.getExtends()
        ix_start, iy_start = self.getGridCellIndicesAt(mesh_xmin, mesh_ymin)
        ix_end, iy_end = self.getGridCellIndicesAt(mesh_xmax, mesh_ymax)

        for ix in range(ix_start, ix_end+1):
            for iy in range(iy_start, iy_end+1):
                self.cells[ix][iy].append(mesh)

    def traverseGrid(self,o, d):
        left, bottom, right, top, xres, yres, w, h = self.getAllDims()

        # ==================  does collide with box
        tx_left = (left - o[0]) / d[0]
        tx_right = (right - o[0]) / d[0]
        ty_bottom = (bottom - o[1]) / d[1] # bottom = min
        ty_top = (top - o[1]) / d[1]       # top = max
        
        if d[0] < 0:
            tx_min = tx_right
            tx_max = tx_left
        else:
            tx_min = tx_left
            tx_max = tx_right

        if d[1] < 0:
            ty_min = ty_top
            ty_max = ty_bottom
        else:
            ty_min = ty_bottom
            ty_max = ty_top

        t0 = max(tx_min, ty_min)
        t1 = min(tx_max, ty_max)

        if not (t0 < t1 and t1 > 0):
            return

        # =========== isHIt
        if self.isInside(o):
            p = o
        else:
            p = o + t0 * d
        
        ix, iy = self.getGridCellIndicesAt(*p)

        # =============== Grid Intersects
        dtx = (tx_max - tx_min) / xres;
        dty = (ty_max - ty_min) / yres;

        if d[0] > 0:
            gridtx = ix + 1
            ix_step = 1
            ix_stop = xres
            tx_next = tx_min + gridtx * dtx
        elif d[0] < 0:
            gridtx = xres - ix
            ix_step = -1
            ix_stop = -1
            tx_next = tx_min + gridtx * dtx
        else:
            ix_step = -1
            ix_stop = -1
            tx_next = FLOATMAX;

        if d[1] > 0:
            gridty = iy + 1
            iy_step = 1
            iy_stop = yres
            ty_next = ty_min + gridty * dty
        elif d[1] < 0:
            gridty = yres - iy
            iy_step = -1
            iy_stop = -1
            ty_next = ty_min + gridty * dty
        else:
            ty_next = FLOATMAX;
            iy_step = -1
            iy_stop = -1

        # traversial as generator (cleaner for loops for multiple passes)
        def gridTravesial(ix, iy, tx_next, ty_next):
            while ix != ix_stop and iy != iy_stop:
                yield (ix, iy, tx_next, ty_next)
                if tx_next < ty_next:
                    ix += ix_step
                    tx_next += dtx
                else:
                    iy += iy_step
                    ty_next += dty

        # go through all geometries in cell and intersect

        cellStyles = [] # (x, y, syle), cashe cell coordinates and draw style /traverse information
        overall_hits = [] # cashe all hits to get the closest hit, normally you could break at the first one because of spacial coherance, but we keep on traversing for drawing
        for index_x, index_y, my_tx_next, my_ty_next in gridTravesial(ix, iy, tx_next, ty_next):
            geometries = self.cells[index_x][index_y]

            mailbox_skip = 0
            if len(geometries) > 0: 
                hitsInCell = []
                for geometry in geometries:
                    if not hasattr(geometry, 'isMailboxed'): # only non mailboxed geometries
                        hit_info = geometry.intersect(self.ct, o, d)
                        if type(hit_info) == tuple: # we have valid hit?
                            hit_t, _ = hit_info
                            if hit_t < my_tx_next and hit_t < my_ty_next: # only count hits located in cur cell
                                hitsInCell.append(hit_t)
                                geometry.isMailboxed = True # mailbox future intersections, normally not needed, because you could break
                            else:
                                pass
                        else:
                            geometry.isMailboxed = True      # missed geometry, no need to reintersect in future
                    else:
                        mailbox_skip += 1 # skipped because of mailbox

                # some geoms where hit, get min for cell
                if len(hitsInCell) > 0:
                    t_min = min(hitsInCell)
                    cellStyles.append( (index_x, index_y, 'hit') )
                    overall_hits += hitsInCell
                    #break;
                elif mailbox_skip > 0:
                    cellStyles.append( (index_x, index_y, 'mailskip') )
                else:
                    cellStyles.append( (index_x, index_y, 'tested') )
            else:
                # cell is empty
                cellStyles.append( (index_x, index_y, 'traversed') )

        def markCells():
            for x_index, y_index, style in cellStyles:
                if style == 'hit':
                    self.markCell(x_index, y_index,(Hue.GREEN, 0.3, .9) )
                elif style == 'traversed':
                    self.markCell(x_index, y_index, (Hue.AZURE, 0.1, 1))
                elif style == 'tested':
                    self.markCell(x_index, y_index, (Hue.AZURE, 0.3, 1))
                elif style == 'mailskip':
                    self.markCell(x_index, y_index, (Hue.RED, 0.3, 1))
        
        drawQueue.addCall(2, markCells)

        if len(overall_hits):
            def drawRayObjectIntersection():
                self.ct.save()
                self.ct.set_source_rgb(0, 0, 0)
                hitpoint = np.add(o, np.multiply(min(overall_hits), d))
                point(self.ct, hitpoint)
                self.ct.restore()

            drawQueue.addCall(5, drawRayObjectIntersection)

        # intersections with grid cells
        def drawGridIntersections():
            for _, _, my_tx_next, my_ty_next in gridTravesial(ix, iy, tx_next, ty_next):
                ct = self.ct
                ct.save()
                ct.set_source_rgb(0, 0, 1)
                point(self.ct, np.add(o, np.multiply(my_tx_next, d)))
                ct.set_source_rgb(1, 0, 0)
                point(self.ct, np.add(o, np.multiply(my_ty_next, d)))
                ct.restore()
        drawQueue.addCall(2, drawGridIntersections)

    def drawFilledCells(self):
        for x_index in range(self.gridres[0]):
            for y_index in range(self.gridres[1]):
                if len(self.cells[x_index][y_index]) > 0:
                    self.markCell(x_index, y_index, (Hue.YELLOW, 0.1, .93))

    def markCell(self, ix, iy, hsvcolor=None):
        left, bottom, _, _, xres, yres, w, h = self.getAllDims()

        cellsize_x = (w / xres)
        cellsize_y = (h / yres)
        ct = self.ct
        ct.save()

        if type(hsvcolor) == tuple:
            ct.set_source_rgb(*hsv_to_rgb(*hsvcolor))
        else:
            ct.set_source_rgb(*hsv_to_rgb(Hue.AZURE, 0.15, 1))

        ct.rectangle(ix*cellsize_x + left, iy*cellsize_y + bottom, cellsize_x, cellsize_y)
        ct.fill()
        ct.restore()

    def isInside(self, o):
        left, bottom, right, top, _, _, _, _ = self.getAllDims()
        # jetzt nur noch die Formel und dann gucken, wo landet!
        if o[0] > left and o[0] < right and o[1] > bottom and o[1] < top:
            # box anmalen!
            # getCell
            return True
        return False

    def drawGrid(self) :
        #arrow(ct, (0,0), gridpos)
        x, y, _, _, xres, yres, w, h = self.getAllDims()
        dx = w / xres
        dy = h / yres
        
        ct = self.ct
        ct.save()
        ct.set_line_width(ct.get_line_width() * 0.5)
        ct.set_dash([.1])
        for i in range(xres+1):
            ct.move_to(x + i*dx , y)
            ct.rel_line_to(0, h)
            ct.stroke()
        for j in range(yres+1):
            ct.move_to(x, y + j*dy)
            ct.rel_line_to(w, 0)
            ct.stroke()
        ct.restore()


if __name__=='__main__':

    def myClick(x, y, b):
        global rayOrigin, rayDir
        if b == 1:
            rayOrigin = (x, y)
        elif b == 3:
            if not (rayOrigin[0] == x and rayOrigin[1] == y):
                rayDir = normalize(np.subtract((x, y), rayOrigin))

    def myscene(ct):
        drawQueue.reset()
        # geometries
        meshes = [
            ClosedMesh(getCircleVertices(1, res = 3)), 
            ClosedMesh([ (p[0] + 5, p[1] + 5) for p in getCircleVertices(1, res = 4)]),
            ClosedMesh([ (p[0] - 5, p[1] - 5) for p in getCircleVertices(1, res = 5)]),
            #
            ClosedMesh([ (p[0] - 2, p[1] - 1) for p in getCircleVertices(1, res = 6)]),
            ClosedMesh([ (p[0] + 3, p[1] - 3) for p in getCircleVertices(1, res = 7)]),
            #   
            ClosedMesh([ (p[0] - 7, p[1] - 1) for p in getCircleVertices(1, res = 3)]),
            ClosedMesh([ (p[0] + 2, p[1] - 5) for p in getCircleVertices(1, res = 3)]),
            #
            ClosedMesh([ (p[0] - 5, p[1] + 5) for p in getCircleVertices(2, res = 8)]),
            ClosedMesh([ (p[0], p[1]-2) for p in kasten(3) ]),
        ]
        # Grid
        #grid = Grid(ct, gridres = (6, 6), gridpos = (-6, -6), gridsize=(12, 12))
        grid = Grid.setupGrid(ct, meshes)
        grid.traverseGrid(rayOrigin, rayDir)

        # draw grid and coordinate frame
        drawQueue.addCall(3, lambda : grid.drawGrid())
        drawQueue.addCall(3, lambda : drawCoordinateFrame(ct))

        # vector origin and direction
        drawQueue.addCall(3, lambda : point(ct, rayOrigin))
        drawQueue.addCall(3, lambda : vector(ct, rayOrigin, rayDir) )

        drawQueue.addCall(0, lambda : grid.drawFilledCells() )
        def shit():
            for mesh in meshes:
                mesh.draw(ct, (0,0))
        drawQueue.addCall(4, shit)

        # dashed vector line 
        def drawDashedLine():
            ct.save()
            ct.set_line_width(ct.get_line_width() / 2)
            ct.set_dash([.1])
            lineSection(ct, np.add(rayOrigin, rayDir), np.add(rayOrigin, np.multiply(rayDir, 20)))
            ct.restore()
        drawQueue.addCall(5, drawDashedLine)

        drawQueue.draw()

    global rayOrigin, rayDir
    rayOrigin = (0, 3.5)
    rayDir = normalize((1, 0))  # drin, negativ, von rechts
    Window((500,500), myscene, myClick)
