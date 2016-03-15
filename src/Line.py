from Point import Point

class Line:

    point1 = Point()
    point2 = Point()

    def __init__(self, point1, point2):
        self.point1=point1
        self.point2=point2

    def x1(self):
        return(self.point1.x)

    def x2(self):
        return(self.point2.x)

    def y1(self):
        return(self.point1.y)

    def y2(self):
        return(self.point2.y)

    def __str__(self):
        return "("+str(self.x1())+","+str(self.y1())+"),("+str(self.x2())+","+str(self.y2())+")"
