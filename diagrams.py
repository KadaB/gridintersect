#!/usr/bin/env python3
import cairo
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

import numpy as np
import numpyadd as npa

from numpy import linalg as la
from colorsys import hsv_to_rgb

def lerp(a, b, alpha):
    return (1-alpha) * a + alpha * b

class Hue:
    RED = 0
    ORANGE = 30/360
    YELLOW = 60/360
    OLIVE = 90/360
    GREEN = 120/360
    EMERALD = 150/360
    CYAN = 180/360
    AZURE = 210/360
    BLUE = 240/360
    VIOLET=270/360
    MAGENTA = 300/360

def solve2x2(a, b, c): # x*A + y*B = C
    M = np.array( [a, b] ).transpose()
    if la.det(M) != 0:
        return la.solve(M, np.array(c) )
    else: 
        return None, None

def getCircleVertices(r, res = 20):
    pis = reversed([np.pi*2/res*i for i in range(res)])
    return [ (np.cos(rad) * r, np.sin(rad) *r ) for rad in pis]

def kasten(size, wall_thickness = .3):
    a = -size
    b = a - wall_thickness
    return [ (a, a), (-a, a), (-a, -a), (a, -a), (a, -b), (-b, -b), (-b,b), (a, b)]

def minIntersect(intersections):
    x = [ v for v in intersections if type(v) == tuple ]
    return min(x) if len(x) > 0 else None
def getOrthoVec(v):
    return np.array([-v[1], v[0]])
            
def point(ct, pos):
    linewidth = ct.get_line_width()
    ct.save()
    ct.translate(*pos)
    ct.scale(linewidth, linewidth)

    ct.arc(0, 0 , 3, 0, 2*np.pi)

    ct.set_line_width(1)
    ct.fill()
    ct.restore() 

def drawCoordinateFrame(ct):
    linewidth = ct.get_line_width()
    ct.save()
    arrow(ct, (0, 0), (0,9))
    arrow(ct, (0, 0), (0,-9))
    for i in range(-9, 10):
        ct.save()
        ct.translate(0, i)
        ct.scale(linewidth, linewidth)
        ct.set_line_width(1)
        ct.move_to(-2, 0)
        ct.line_to(2, 0)
        ct.stroke()
        ct.restore()

    arrow(ct, (0, 0), (9, 0))
    arrow(ct, (0, 0), (-9, 0))
    for i in range(-9, 10):
        ct.save()
        ct.translate(i, 0)
        ct.scale(linewidth, linewidth)
        ct.set_line_width(1)
        ct.move_to(0, -2)
        ct.line_to(0, 2)
        ct.stroke()
        ct.restore()
    ct.restore()

def lineSection(ct, start, end):
    ct.move_to(*start)
    ct.line_to(*end)
    ct.stroke()

def arrow_head(ct, pos, dir_vec, size=1):#, setback, angle = np.pi/4, dir_vec):
    linewidth = ct.get_line_width()
    dir_vec_ortho = (-dir_vec[1], dir_vec[0])
    a = npa.normalize(dir_vec)
    b = npa.normalize(dir_vec_ortho)
    
    ct.save()
    ct.transform(cairo.Matrix(*a, *b, *pos)) # orient arrow head
    ct.scale(linewidth*size, linewidth*size)           # set coordinates in rel to lw
    ct.set_line_width(1)                     # keep world space lw

    #draw geometry
    ct.move_to(0,2)
    ct.line_to(10,0)
    ct.line_to(0,-2)
    ct.close_path()
    ct.fill()

    ct.restore()

def arrow(ct, start, end, head_size=1):
    lineSection(ct, start, end)
    vec = np.subtract(end, start)
    arrow_head(ct, end, vec, head_size)

# o stützvector, d richtungsvektor, l länge 
def vector(ct, o, d, l = 1, head_size=1):
    arrow(ct, o, np.add(o, np.multiply(d, l)), head_size)

def polygon(ct, pos, vertices, hsv_color=None):
    ct.save()
    ct.translate(*pos)
    ct.move_to(*(vertices[0]))
    for p in vertices[1:]:
        ct.line_to(*p)
    ct.close_path()
    ct.stroke_preserve()
    if type(hsv_color) == tuple:
        ct.set_source_rgb(*hsv_to_rgb(*hsv_color))
        ct.fill()
    ct.restore()

class Window(Gtk.Window):
    def __init__(self, res, drawFunction=None, onClick=None):
        super().__init__()
        self.connect('destroy', lambda w: Gtk.main_quit())

        def onkey(w, e):
            if e.keyval == Gdk.KEY_q:
                Gtk.main_quit()
            if e.keyval == Gdk.KEY_Return:
                with cairo.SVGSurface("diagram.svg", *(w.get_size())) as surface:
                    context = cairo.Context(surface)
                    self.draw(None, context)
                    print("file written")
        
        self.connect("key-press-event", onkey)


        self.set_default_size(*res)

        drawingarea = Gtk.DrawingArea()
        drawingarea.connect('draw', self.draw)
        drawingarea.connect('button-press-event', self.onButtonPress)
        drawingarea.set_events(Gdk.EventMask.BUTTON_PRESS_MASK)
        self.add(drawingarea)
        self.show_all()
        self.drawFunction = drawFunction
        self.clickFunction = onClick

        self.coordSize = 10 # from -10 to 10
        Gtk.main()

    def getRes(self):
        coordSize = self.coordSize
        w, h = self.get_size()
        aspect = w / h
        res = ((w/coordSize)/2/aspect,-(h/coordSize)/2) if aspect >= 1 else ((w/coordSize)/2,-(h/coordSize)/2*aspect)
        center = (w / 2, h / 2)
        return res, center

    def onButtonPress(self, w, e):
        if type(self.clickFunction) != type(None):
            res, center = self.getRes()
            M = np.array([ [res[0], 0, center[0]], 
                            [0, res[1], center[1]],
                            [0, 0, 1] ])
            p = npa.inv(M).dot((e.x, e.y, 1))
            self.clickFunction(p[0], p[1], e.button)
            w.queue_draw()

    def draw(self, da, ct):
        w, h = self.get_size()

        linewidth = self.coordSize*2/min(w, h) # keep line width same for every zoom level
        ct.set_line_width(linewidth) # typical value: 0.036

        res, center = self.getRes()
        ct.save()
        # todo: richtige matrix für top,left,bottom,right
        ct.transform(cairo.Matrix(res[0], 0, 0, res[1], center[0], center[1]))
        if type(self.drawFunction) != type(None):
            self.drawFunction(ct)

        ct.restore()

def myscene(ct):
    # geometry
    drawCoordinateFrame(ct)

def main():
    win = Window((650, 550), myscene)

if __name__ == '__main__':
    main()
