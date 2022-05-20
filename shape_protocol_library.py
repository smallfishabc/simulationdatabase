# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 17:43:42 2022

@author: ShaharGroup-fyu
"""
import mdtraj as md
import numpy as np
import sys
import entropy_library as el

# Here we can set a distance between the cone tip and the first alpha carbon
def compute_forbidden_distance(distance,angle,trajectory,location_of_alpha,number_of_frames):
    # We can move the following part outside this function for better debugging
    # efficiency in the Jupyter notebook
    k=0
    current_frame=0
    prohibited_frames=[]
    # In previous function, we consider the cone's apex is the first alpha carbon.
    # Here, we defined a distance between the cone's apex and the first alpha carbon. 
    # However, the symmetry axis of the cone still goes through the first and the second alpha carbon
    # Therefor, the angle limit should be calculate according to each atom.
    while current_frame<number_of_frames:
        # The xyz coordinate of the first and the second alpha carbon
        a1=trajectory.xyz[current_frame,location_of_alpha[0],:]
        a2=trajectory.xyz[current_frame,location_of_alpha[1],:]
        # First alpha-carbon is our starting point
        # The second alpha-carbon determine our standard vector. We can change this defination
        standard_v=el.vector(a1,a2)
        standard_v_norm=el.maginitude(standard_v)
        for i in location_of_alpha:
            ak=trajectory.xyz[current_frame,i,:]
            vk=el.vector(a1,ak)
            # Skip the starting residue
            if i>location_of_alpha[1]:
                #vk_norm=el.maginitude(vk)
                d_k=el.projection(vk,standard_v)
                d_kb=el.maginitude(np.cross(standard_v,vk))/standard_v_norm
                d_limit=el.limit_distance(distance,angle,d_k)
            else:
                continue
            if d_kb>=d_limit:
                prohibited_frames.append(current_frame)
                k+=1
                break
        current_frame+=1
        ratio=(number_of_frames-k)
    return (prohibited_frames,ratio)


# Enhanced sampling
    # Compared to other scripts, here the angle_theta = (np.pi-angle_theta
def compute_forbiden_rotation_enhanced(distance,angle,trajectory,location_of_alpha,number_of_frames):
    current_frame=0
    # Rotation Interval 
    # This number determined
    rot_angle=60
    # How many rotations should be made
    rot_repeat=int(360/rot_angle)-1
    # Convert degree to rad
    rot_angle_pi=rot_angle/360*2*3.141593
    # Create a np array to save memory processing time
    #This line is for the debugging which will give out the residue that is forbiddened
    #result=np.zeros((number_of_frames,rot_repeat+1,2))
    #This line will only give out the True or False value
    result=np.zeros((number_of_frames,rot_repeat+1))
    while current_frame<number_of_frames:
        objectframe=trajectory.xyz[current_frame,:,:]
        # Shift the object to make a1 as the origin
        object_shifted=el.change_origin(objectframe,objectframe[location_of_alpha[0],:])
        a1a2=objectframe[location_of_alpha[1],:]-objectframe[location_of_alpha[0],:]
        #print(a1a2)
        a1a3=objectframe[location_of_alpha[2],:]-objectframe[location_of_alpha[0],:]
        #print(a1a3)
        a1a2a3_norm_vec=el.unit_normal_vector(a1a2,a1a3)
        #print(a1a2a3_norm_vec)
        # Rotate a1a2 by theta along the norm_vec we can get a1a2p
        # Which is parallel to the constraint surface
        a1a2p=el.vector_rotation(a1a2a3_norm_vec,a1a2,angle)
        #print(a1a2p)
        # Here we created a plane that is parallel to the constrain plane and go through the a1
        a1_r_plane_norm_vec=el.unit_normal_vector(a1a2p,a1a2a3_norm_vec)
        #print(a1_r_plane_norm_vec)
        # The distance between a2 and the constraint plane should be larger than given distance
        sign_vector=np.dot(a1a2,a1_r_plane_norm_vec)
        sign=np.sign(sign_vector)
        #now = datetime.now().time()
        #print("now =", now)
        current_rotation=0
        rot_mat=el.rotation_matrix(rot_angle_pi,a1a2)
        conformation_rotated=object_shifted
        while current_rotation <= rot_repeat:
            limit_ak=[True,0]
            #print(current_rotation)
            for index,i in enumerate(location_of_alpha[2:]):
                ak=conformation_rotated[i,:]
                #print(ak)
                if not el.limitation(sign,a1_r_plane_norm_vec,distance,ak):
                    limit_ak=[False,(index+4)]
                    break
            conformation_rotated=np.matmul(conformation_rotated,rot_mat.T)
            #This line is for the debugging which will give out the residue that is forbiddened
            #result[current_frame,current_rotation]=limit_ak
            #This line will only give out the True or False value
            result[current_frame,current_rotation]=limit_ak[0]
            #print(limit_ak)
            #print(result[current_frame,current_rotation])
            current_rotation+=1
        current_frame+=1   
        #now = datetime.now().time()
        #print("now =", now)
        #print(result.shape)
    return result

# Sphere curvature protocol
# Compute how many possible conformations is prohibitted by constraint surface (Using curvature to simulate the membrane)
# Here we define curvature using the radius. We draw a sphere as the constraint surface.
# The C center of the sphere is along the same line as vector a1a2. The distance between a1 and C is the radius of the sphere
def compute_forbidden_curvature(radius, trajectory, location_of_alpha, number_of_frames):
    # We can move the following part outside this function for better debugging
    # efficiency in the Jupyter notebook
    k = 0
    current_frame = 0
    allowed_frames = []
    prohibited_frames = []
    prohibited_test = []
    # Here, we defined the radius as the distance between the C center of the sphere and the first alpha carbon.
    while current_frame < number_of_frames:
        # The xyz coordinate of the first and the second alpha carbon
        # Make a1 as the origin of the coordinate system
        objectframe = trajectory.xyz[current_frame, :, :]
        # Shift the object to make a1 as the origin
        object_shifted = el.change_origin(objectframe, objectframe[location_of_alpha[0], :])
        a1a2 = objectframe[location_of_alpha[1], :] - objectframe[location_of_alpha[0], :]
        # First alpha-carbon is our starting point
        # The second alpha-carbon determine our standard vector and then we get a unit vector.
        standard_v_unit_vector = a1a2 / el.maginitude(a1a2)
        # Calculate the location of the center of the sphere
        c = -1 * radius * standard_v_unit_vector
        # Here we only consider that every alpha carbon should locate outside the sphere.
        # In the not simplified version, we considered that the bonds between each atom should also locate outside the sphere.
        # We can further identify whether every atom is outside the sphere in the future
        for index, i in enumerate(location_of_alpha[2:]):
            ak = object_shifted[i, :]
            dist_ak_c = np.linalg.norm(ak - c)
            if dist_ak_c <= radius:
                # Here we start from a3 (index=0). a1 in the vmd is residue 2 so a3 is residue 4. So we plus 4 here
                # For VMD debugging only
                # prohibited_frames.append([current_frame, index + 4])
                prohibited_frames.append(current_frame)
                k += 1
                break
        current_frame += 1
    return prohibited_frames, k