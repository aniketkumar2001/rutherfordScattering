
import random as r
import numpy as np
import pygame
from astropy import constants as const
from sys import exit
pygame.init()
clock = pygame.time.Clock()

SCREEN = pygame.display.set_mode((1550, 800))          
f_rate = 60

p1 = pygame.Surface((2,2))
p1.fill((255,255,255))
p2 = pygame.Surface((3,3))
p2.fill((255,0,0))


aidedfactor = 6e-15         #magnification factor
dt = 1e-22                  #time-step
y_center=(SCREEN.get_height()//2)*aidedfactor
x_center=(SCREEN.get_width()//2)*aidedfactor

# positions of gold nuclei
p2_position = np.array([[x_center-0.5e-12,y_center-0.5e-12], [x_center,y_center-0.5e-12], [x_center+0.5e-12,y_center-0.5e-12],
                        [x_center-0.5e-12,y_center], [x_center,y_center], [x_center+0.5e-12,y_center],
                        [x_center-0.5e-12,y_center+0.5e-12], [x_center,y_center+0.5e-12], [x_center+0.5e-12,y_center+0.5e-12]
                        ])


totalBodies = 20 # no. of incoming bodies 

# parameters for initial positions of incoming particles
xrange = 0.5e-11
yrange = 1e-12

# Define physical constants
epsilon_0 = const.eps0.value  # Vacuum permittivity
q_alpha = 2 * const.e.si.value  # Charge of alpha particle
m_alpha = 4 * const.m_p.si.value  # Mass of alpha particle
q_gold = 79*const.e.si.value # charge of gold nucleus
k = (q_alpha * q_gold)/(4 * np.pi * epsilon_0 * m_alpha)


# incoming particles list

obj0=[]
i=1
while i<=totalBodies: 
    randomNo = (r.randrange(-100,100)/100)*yrange
    obj0.append([np.array([x_center-xrange,y_center+randomNo]),np.array([1e7,0])])
    i+=1
obj0=np.array(obj0, float)


# function for interaction b/w particles and nuclei

def addG(obj,nucleuspos): 
    s=obj[0]-nucleuspos
    sMod=np.linalg.norm(s)
    acc=(k*s)/(sMod**3)
    obj[1]+=[acc[0]*dt,acc[1]*dt]
    obj[0]+=obj[1]*dt
    return obj   


# main:

while True :
    for e in pygame.event.get():
        if (e.type == pygame.QUIT):
            pygame.quit()
            exit()
    SCREEN.fill((0,0,0))   #comment this line to see the path of incoming particles
    for j in p2_position:
        SCREEN.blit(p2,j/aidedfactor)
    
    i=0
    while i<totalBodies:
        SCREEN.blit(p1,obj0[i][0]/aidedfactor)  
        for j in p2_position:      
            obj0[i] = addG(obj0[i],j)
        i+=1

    pygame.display.update()
    clock.tick(f_rate)