from Bio import Nexus
from Bio import Phylo
#from PIL import Image, ImageDraw, ImageChops, ImageFont
#import numpy
import math
from Point import Point
from Line import Line
from Cache import Cache
import copy
import re

try:
    import cairocffi as cairo
except:
    print "cairocffi not present, trying loading cairo"
    import cairo


class RadialTree:
    """
    Render an unrooted tree
    Adapted in python/Pillow from figtree/treeviewer/treelayouts/RadialTreeLayout.java
    https://github.com/rambaut/figtree/
    
    Copyright (C) 2006-2014 Andrew Rambaut
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    """
    spread = 0.0;
    #fnt_large = ImageFont.truetype('/Library/Fonts/Microsoft/Calibril.ttf', 35)
    #fnt_small = ImageFont.truetype('/Library/Fonts/Microsoft/Calibril.ttf', 25)
    fnt_large = 14 #ImageFont.truetype('Pillow/Tests/fonts/DejaVuSans.ttf', 14)
    fnt_small = 10 #ImageFont.truetype('Pillow/Tests/fonts/DejaVuSans.ttf', 10)

    tip_border= 20
    tip_distance=10
    image_border=0
    point_radius=4

    def render_png(self,tree, width, height, out_file):
        #image = Image.new('RGBA', (width*2,height*2), (255,255,255,255))
        image = cairo.ImageSurface (cairo.FORMAT_ARGB32, width,height)
        d = cairo.Context (image)
        #d.scale (width, height) # Normalizing the canvas
        #d = ImageDraw.Draw(image)
        cache = Cache()
        root = tree.node(tree.root)
        hasBranchLengths = self.total_branchlength(tree)
        print str(hasBranchLengths)+" "+str(self.total_branchlength(tree))
        self.constructNode(tree, root, 0.0, math.pi * 2, 0.0, 0.0, 0.0,hasBranchLengths, cache)
        self.drawTree(cache,d,width,height)
        #image = image.resize((width,height), Image.ANTIALIAS)
        #image.save(out_file, "PNG")
        d.stroke ()
        image.write_to_png (out_file) # Output to PNG

    def render_pdf(self,tree, width, height, out_file):
        image = cairo.PDFSurface(out_file,width,height)
        d = cairo.Context (image)
        cache = Cache()
        root = tree.node(tree.root)
        hasBranchLengths = self.total_branchlength(tree)
        self.constructNode(tree, root, 0.0, math.pi * 2, 0.0, 0.0, 0.0,hasBranchLengths, cache)
        self.drawTree(cache,d,width,height)
        d.stroke ()

    def taxNumber(self,tree,node):
        num = 0
        if(len(node.succ)==0):
            num = 1
        else:
            for n in node.succ:
                num += self.taxNumber(tree,tree.node(n))
        return num

    def constructNode(self,
                      tree, node,
                      angleStart, angleFinish,
                      xPosition, yPosition, length, hasBranchLengths, cache):
        branchAngle = (angleStart + angleFinish) / 2.0

        directionX = math.cos(branchAngle)
        directionY = math.sin(branchAngle)
        nodePoint = Point(xPosition + (length * directionX), yPosition + (length * directionY))

        if(len(node.succ)!=0) :
            children = node.succ
            leafCounts = []
            sumLeafCount = 0
            
            for n in children :
                child=tree.node(n)
                numT = self.taxNumber(tree,child)
                leafCounts.append(numT)
                sumLeafCount += numT

            span = (angleFinish - angleStart)

            if (node.get_id() != tree.root):
                span *= 1.0 + (self.spread / 10.0)
                angleStart = branchAngle - (span / 2.0)
                angleFinish = branchAngle + (span / 2.0)

            a2 = angleStart

            rotate = False;
            for i in range(0,len(children)):
                index = i
                if (rotate):
                    index = len(children.size()) - i - 1

                child = tree.node(children[index])

                childLength = child.data.branchlength
                if not hasBranchLengths:
                    childLength= 1.0

                a1 = a2;
                a2 = a1 + (span * leafCounts[index] / sumLeafCount)

                childPoint = self.constructNode(tree, child, a1, a2,
                                                nodePoint.x, nodePoint.y, childLength,hasBranchLengths, cache)

                branchLine = Line(childPoint,nodePoint)

                # add the branchLine to the map of branch paths
                cache.branchPaths[child] = branchLine
                cache.branchLabelPaths[child] =copy.copy(branchLine)

            nodeLabelPoint = Point(xPosition + ((length + 1.0) * directionX),
                                   yPosition + ((length + 1.0) * directionY))

            nodeLabelPath = Line(nodePoint, nodeLabelPoint)
            cache.nodeLabelPaths[node] = nodeLabelPath
        else:
            taxonPoint = Point(xPosition + ((length + 1.0) * directionX),
                               yPosition + ((length + 1.0) * directionY))

            taxonLabelPath = Line(nodePoint, taxonPoint);
            cache.tipLabelPaths[node] =taxonLabelPath

        nodeShapePoint = Point(xPosition + ((length - 1.0) * directionX),
                               yPosition + ((length - 1.0) * directionY))
        nodeShapePath = Line(nodePoint, nodeShapePoint)
        cache.nodeShapePaths[node] = nodeShapePath
        cache.nodePoints[node] = nodePoint
        print nodeShapePath

        return nodePoint

    def drawTree(self,cache,image_draw,width, height):
        xscale=self.xscale(cache,width)
        yscale=self.yscale(cache,height)
        scale = min(xscale,yscale)
        xoffset=self.xoffset(cache,width,scale)
        yoffset=self.yoffset(cache,height,scale)

        for line in cache.branchPaths.values():
            image_draw.set_source_rgb (0, 0, 0) # Solid color
            image_draw.set_line_width (2)
            image_draw.move_to(line.x1()*scale+xoffset,line.y1()*scale+yoffset);
            image_draw.line_to(line.x2()*scale+xoffset,line.y2()*scale+yoffset);

            #image_draw.line([(line.x1()*xscale+xoffset,line.y1()*yscale+yoffset),(line.x2()*xscale+xoffset,line.y2()*yscale+yoffset)],(0,0,0,255),2)
            
        for node,line in cache.tipLabelPaths.iteritems():
            # Knowing the line from the taxon tip : we find the point in the same direction, at distance 15
            point3 = self.label_coordinates(line.point1,line.point2,self.tip_distance,scale,xoffset,yoffset)
            image_draw.set_source_rgb (0.0, 0.0, 0.0)
            image_draw.select_font_face ("sans-serif")#,cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            image_draw.set_font_size(self.fnt_large)
            x_bearing, y_bearing, w, h = image_draw.text_extents(node.data.taxon)[:4]
            
            image_draw.move_to (point3.x - w / 2 - x_bearing,
                                point3.y - h / 2 - y_bearing)
            image_draw.show_text (node.data.taxon)
            #tw,th = image_draw.textsize(node.data.taxon, font=self.fnt_small)
            #image_draw.text([line.x1()*xscale+xoffset-tw/2,line.y1()*yscale+yoffset-th/2], node.data.taxon, (0,0,0,255), font=self.fnt_large)

        # for node,line in cache.nodeLabelPaths.iteritems():
        #     image_draw.ellipse([line.x1()*xscale+xoffset-self.point_radius,
        #                         line.y1()*yscale+yoffset-self.point_radius,
        #                         line.x1()*xscale+xoffset+self.point_radius,
        #                         line.y1()*yscale+yoffset+self.point_radius],
        #                        fill=(0,0,0,255))
            

        # for line in cache.nodeLabelPaths.values() :
        #     image_draw.line([((line.x1()+xoffset)*xscale,(line.y1()+yoffset)*yscale),((line.x2()+xoffset)*xscale,(line.y2()+yoffset)*yscale)],(0,0,0,255),2)
        # for line in cache.tipLabelPaths.values():
        #     image_draw.line([((line.x1()+xoffset)*xscale,(line.y1()+yoffset)*yscale),((line.x2()+xoffset)*xscale,(line.y2()+yoffset)*yscale)],(0,0,0,255),2)
        # for line in cache.nodeShapePaths.values():
        #     image_draw.line([((line.x1()+xoffset)*xscale,(line.y1()+yoffset)*yscale),((line.x2()+xoffset)*xscale,(line.y2()+yoffset)*yscale)],(0,0,0,255),2)

    def xscale(self,cache,width):
        xmin=100000
        xmax=0
        for line in cache.branchPaths.values():
            xmin = min(xmin,line.x1(),line.x2())
            xmax = max(xmax,line.x1(),line.x2())
        # for line in  cache.nodeLabelPaths.values() :
        #     xmin = min(xmin,line.x1(),line.x2())
        #     xmax = max(xmax,line.x1(),line.x2())
        # for line in cache.tipLabelPaths.values():
        #     xmin = min(xmin,line.x1(),line.x2())
        #     xmax = max(xmax,line.x1(),line.x2())
        # for line in cache.nodeShapePaths.values():
        #     xmin = min(xmin,line.x1(),line.x2())
        #     xmax = max(xmax,line.x1(),line.x2())
        return (width-self.tip_border*2-self.image_border*2)*1.0/(xmax-xmin)

    def yscale(self,cache,height):
        ymin=100000
        ymax=0
        for line in cache.branchPaths.values():
            ymin = min(ymin,line.y1(),line.y2())
            ymax = max(ymax,line.y1(),line.y2())
        # for line in cache.nodeLabelPaths.values() :
        #     ymin = min(ymin,line.y1(),line.y2())
        #     ymax = max(ymax,line.y1(),line.y2())
        # for line in cache.tipLabelPaths.values():
        #     ymin = min(ymin,line.y1(),line.y2())
        #     ymax = max(ymax,line.y1(),line.y2())
        # for line in cache.nodeShapePaths.values():
        #     ymin = min(ymin,line.y1(),line.y2())
        #     ymax = max(ymax,line.y1(),line.y2())
        return (height-self.tip_border*2-self.image_border*2)*1.0/(ymax-ymin)

    def xoffset(self,cache,width,xscale):
        xmin = 100000
        for line in cache.branchPaths.values():
            xmin = min(xmin,line.x1()*xscale,line.x2()*xscale)
        # for line in cache.nodeLabelPaths.values() :
        #     xmin = min(xmin,line.x1(),line.x2())
        # for line in cache.tipLabelPaths.values():
        #     xmin = min(xmin,line.x1(),line.x2())
        # for line in cache.nodeShapePaths.values():
        #     xmin = min(xmin,line.x1(),line.x2())
        return -xmin+self.tip_border+self.image_border

    def yoffset(self,cache,height,yscale):
        ymin = 100000
        for line in cache.branchPaths.values():
            ymin = min(ymin,line.y1()*yscale,line.y2()*yscale)
        # for line in cache.nodeLabelPaths.values() :
        #     ymin = min(ymin,line.y1(),line.y2())
        # for line in cache.tipLabelPaths.values():
        #     ymin = min(ymin,line.y1(),line.y2())
        # for line in cache.nodeShapePaths.values():
        #     ymin = min(ymin,line.y1(),line.y2())
        return -ymin+self.tip_border+self.image_border


    def total_branchlength(self, tree, node=None):
        if node == None:
            node = tree.node(tree.root)
        sumbr=node.data.branchlength
        for n in node.succ:
            sumbr+=self.total_branchlength(tree,tree.node(n))
        return sumbr

    "Returns the point in the plane, at distance d from point 1, in the direction of point 2" 
    def label_coordinates(self,p1, p2, d,scale,xoffset,yoffset):
        x1 = p1.x*scale+xoffset
        y1 = p1.y*scale+yoffset
        x2 = p2.x*scale+xoffset
        y2 = p2.y*scale+yoffset
        lengthp1p2 = math.sqrt(math.pow(x2 - x1,2) + math.pow(y2-y1,2))
        print(lengthp1p2)
        out_x = d * 1.0 * (x2 - x1) / lengthp1p2 + x1
        out_y = d * 1.0 * (y2 - y1) / lengthp1p2 + y1
        return(Point(out_x,out_y))
