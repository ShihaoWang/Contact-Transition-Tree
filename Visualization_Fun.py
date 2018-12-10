import sys, os, time
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
from klampt.model.trajectory import Trajectory
from Seed_Guess_Fun import *
from Contact_Force_Fun import *
from scipy.interpolate import interp1d

# This file contains functions related to

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.starp = False

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed",c)
        if c == 'q':
            self.quit = True
            return True
        if c == 's':
            self.starp = not self.starp
            return True
        return False

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]


def Robot_Motion_Plot(world, DOF, control_force_len, contact_link_dictionary, terr_model, robot_option, grids, optimized_solution):
    # This function is used to plot the robot motion
    # The optimized solution is used to plot the robot motion and the contact forces

    # Initialize the robot motion viewer
    robot_viewer = MyGLPlugin(world)
    # robot_viewer.drawContacts = True
    # robot_viewer.drawSensors = True

    # Here it is to unpack the robot optimized solution into a certain sets of the lists
    Duration, State_List_Array, Control_List_Array, Contact_Force_List_Array = Seed_Guess_Unzip(optimized_solution, DOF, control_force_len, grids)

    interpolated_times = 5
    Inter_Time_Array, Inter_State_Array = Trajectory_Interpolation(Duration, grids, State_List_Array, interpolated_times)
    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

    robot_viewer_time = time.time()
    # robot_sim = robot_viewer.sim
    sim_robot = world.robot(0)
    # sim_robot_controller = robot_sim.controller(0)
    contact_link_list = contact_link_dictionary.keys()
    robot_mass = Robot_Total_Mass(world)
    robot_gravity = robot_mass * 9.81

    force_unit_length = 1
    while vis.shown():
        # This is the main plot program
        # ipdb.set_trace()

        for i in range(0, grids * interpolated_times):
            vis.lock()

            Robotstate_Traj_i = Inter_State_Array[i];
            RobotConfig_Traj_i = Robotstate_Traj_i[0:DOF]
            # print Robotstate_Traj_i
            sim_robot.setConfig(RobotConfig_Traj_i)
            # Now it is the plot of the contact force at the contact extremities
            # Contact_Force_Traj_i = Contact_Force_List_Array[i]
            # # Here Contact_Link_PosNVel_List is a list of dictionaries and the list is
            # Contact_Link_PosNVel_List = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)
            # for i in range(0, len(contact_link_list)):
            #     Contact_Force_Index = 0
            #     Contact_Link_i_PosNVel = Contact_Link_PosNVel_List[i]
            #     Contact_Link_i_Pos_List = Contact_Link_i_PosNVel["Pos"]
            #     for j in range(0, len(Contact_Link_i_Pos_List)):
            #         Contact_Link_i_Pos_j = Contact_Link_i_Pos_List[j]
            #         Contact_Link_i_Pos_j_Contact_Force = Contact_Force_Element_from_Index(Contact_Force_Traj_i, Contact_Force_Index)
            #         contact_start_pos, contact_terminal_pos = Contact_Force_Mag_2_Vec(Contact_Link_i_Pos_j, Contact_Link_i_Pos_j_Contact_Force, robot_gravity, force_unit_length)
            #         vis.add(" ", Trajectory([0, 1], [contact_start_pos, contact_terminal_pos]))
            #         Contact_Force_Index = Contact_Force_Index + 1
            # COMPos_start = sim_robot.getCom()
            # COMPos_end = COMPos_start
            # COMPos_end[2] = COMPos_end[2] + 100
        # vis.add("Center of Mass",  Trajectory([0, 1], [COMPos_start, COMPos_end]))
            vis.unlock()
            time.sleep(0.1)

def Contact_Force_Mag_2_Vec(contact_start_pos, contact_force, robot_gravity, force_unit_length):
    # This function is used to calculate the vector of certain contact from the link contact
    # Since this contact force is a 3 by 1 vector, here a little bit hard coding will be used
    # ipdb.set_trace()
    contact_force_length_offset_x = contact_force[0]/robot_gravity * force_unit_length
    contact_force_length_offset_y = contact_force[1]/robot_gravity * force_unit_length
    contact_force_length_offset_z = contact_force[2]/robot_gravity * force_unit_length

    contact_terminal_pos = contact_start_pos[:]

    contact_terminal_pos[0] = contact_start_pos[0] + contact_force_length_offset_x
    contact_terminal_pos[1] = contact_start_pos[1] + contact_force_length_offset_y
    contact_terminal_pos[2] = contact_start_pos[2] + contact_force_length_offset_z

    return contact_start_pos, contact_terminal_pos

def Trajectory_Interpolation(duration, grids, raw_data_list, interpolated_times):
    # This function is used to interpolate the given trajectories according to some rule
    # Here raw_data is a numpy array where each row is a DOf by 1 array
    raw_number, col_number = raw_data_list.shape          # each row is the data at a certain time
    inter_data_list = np.zeros((raw_number * interpolated_times, col_number))
    for i in range(0, col_number):
        raw_data_i = raw_data_list[:,i]
        raw_data_i_time, raw_data_i_inter = Trajectory_i_Interpolation(duration, grids, raw_data_i, interpolated_times)
        inter_data_list[:,i] = raw_data_i_inter
    return raw_data_i_time, inter_data_list

def Trajectory_i_Interpolation(duration, grids, raw_data, interpolated_times):
    # This function is used for trajectory interpolation given a 1D array
    # The first job is to retrieve the duration array
    x = np.linspace(0, duration, num=grids, endpoint=True)
    y = raw_data
    f_inter = interp1d(x, y, kind='cubic')

    # Here f_inter is a interpolation object

    inter_num = grids * interpolated_times
    xnew = np.linspace(0, duration, num=inter_num, endpoint=True)
    ynew = f_inter(xnew)

    return xnew, ynew
