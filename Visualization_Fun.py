import sys, os
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin


class MyGLViewer(GLSimulationPlugin):
    def __init__(self, world, config_init, velocity_init):
        #create a world from the given files
        self.world = world
        self.robot = world.robot(0)
        self.robot.setConfig(config_init)
        self.robot.setVelocity(velocity_init)
        GLSimulationPlugin.__init__(self,world)

    def control_loop(self):
        pass

    def display(self):
        self.sim.updateWorld()
        self.world.drawGL()

    def mousefunc(self,button,state,x,y):
        #Put your mouse handler here
        #the current example prints out the list of objects clicked whenever
        #you right click
        print "mouse",button,state,x,y
        # if button==2:
        #     if state==0:
        #         print [o.getName() for o in self.click_world(x,y)]
        #     return
        GLSimulationPlugin.mousefunc(self,button,state,x,y)

    def motionfunc(self,x,y,dx,dy):
        return GLSimulationPlugin.motionfunc(self,x,y,dx,dy)
