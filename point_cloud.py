from pylab import *
from scipy import spatial
from sys import stdout
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def threeDvision(tar, cur):
    
    tarX=[]
    tarY=[]
    tarZ=[]
    curX=[]
    curY=[]
    curZ=[]
    for i in range(len(tar)):
        if i%50==0:
            line = tar[i]
            tarX.append(line[0])
            tarY.append(line[1])
            tarZ.append(line[2])

            line2 = cur[i]
            curX.append(line2[0])
            curY.append(line2[1])
            curZ.append(line2[2])
            '''
    tar = tar.T
    cur = cur.T
    tarX=tar[0]
    tarY=tar[1]
    tarZ=tar[2]
    curX=cur[0]
    curY=cur[1]
    curZ=cur[2]
    '''
    ax=plt.subplot(111,projection='3d')
    ax.scatter(tarX,tarY,tarZ,c='g') 
    ax.scatter(curX,curY,curZ,c='r')
    ax.set_zlabel('Z') 
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()

def getRotation(qR):
    q0 = qR[0]
    q1 = qR[1]
    q2 = qR[2]
    q3 = qR[3]
    rotation = asarray([
            [q0**2+q1**2-q2**2-q3**2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
            [2*(q1*q2+q0*q3), q0**2+q2**2-q1**2-q3**2, 2*(q2*q3-q0*q1)],
            [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0**2+q3**2-q1**2-q2**2]
            ])
    return rotation

def getReg(tararray, closearray):
    avetar0 =0
    avetar1 =0
    avetar2 =0    
    avecur0 =0
    avecur1 =0
    avecur2 =0
    avetar0= sum(tararray.T[0])
    avetar1= sum(tararray.T[1])
    avetar2= sum(tararray.T[2])
    avecur0= sum(closearray.T[0])
    avecur1= sum(closearray.T[1])
    avecur2= sum(closearray.T[2])
    '''
    for i in range(len(tararray)):
        tt=tararray[i]
        cc=closearray[i]
        avetar0 += tt[0]
        avetar1 += tt[1]
        avetar2 += tt[2]
    
    
        avecur0 += cc[0]
        avecur1 += cc[1]
        avecur2 += cc[2]
        '''
    up = np.array([avetar0/len(tararray), avetar1/len(tararray), avetar2/len(tararray)])
    ux = np.array([avecur0/len(closearray), avecur1/len(closearray), avecur2/len(closearray)])
    
    #covxx = ((x - avetar)**2 for x in
    #crosscov = (tararray.T.dot(closearray) - up.T.dot(ux))/ len(tararray)
    
    '''
    crosscov = zeros([3,3])
    for i in range(len(tararray)):
        crosscov = crosscov + (tararray[i]-up).dot((closearray[i]-ux).T)
    print crosscov
    crosscov = crosscov/len(tararray)
    '''
    crosscov = zeros([3,3])
    for i in range(len(tararray)):
        crosscov = crosscov + (tararray[i].dot(closearray[i].T))
    crosscov = divide(crosscov,len(tararray))-up.T.dot(ux)
    print crosscov
    crosscov = crosscov/len(tararray)
    
    print crosscov
    tr = crosscov[0][0] + crosscov[1][1] + crosscov[2][2]
    #AA = np.array(crosscov - crosscov.T)
    AA = subtract(crosscov, crosscov.T)
    II = np.array(crosscov + crosscov.T - tr*identity(3))
    #II = subtract(add(corsscov, corsscov.T), multiply(identity(3), tr))

    Q = asarray([
            [tr, AA[1][2], AA[2][0], AA[0][1]],
            [AA[1][2], II[0][0], II[0][1], II[0][2]],
            [AA[2][0], II[1][0], II[1][1], II[1][2]],
            [AA[0][1], II[2][0], II[2][1], II[2][2]]
            ])
    
    '''
    Q = asarray([
            [tr, AA[1][2], AA[2][0], AA[0][1]],
            [AA[1][2], II[0][0], 0, 0],
            [AA[2][0], 0, II[1][1], 0],
            [AA[0][1], 0, 0, II[2][2]]
            ])
    '''
    eigenvalue, eigenvector = eig(Q)
    maxeigenvalue = argmax(eigenvalue)
    qR = eigenvector[:, maxeigenvalue]
    rotation = getRotation(qR)
    trans = (ux.T - (rotation.dot(up.T))).T
    trans = subtract(ux.T, rotation.dot(up.T)).T
    
    #print "*************************"
    #print rotation
    #print trans
    #print "*************************"
    #print "*************************"

    return qR, rotation, trans

def MSE(tar, cur):
    error = 0
    tol = len(tar)
    for i in range(tol):
        error += ((tar[i][0]-cur[i][0])**2 +(tar[i][1]-cur[i][1])**2 +(tar[i][2]-cur[i][2])**2)
    return error/tol

                

def ICP(tar, cur):
    kround = 1
    mse = 10
    threshold = 0.0005
    tararray= np.array(tar)
    rotation=[]
    trans=[]
    T= identity(4)
    closearray=[]
    newtar = tararray
    kdtree = spatial.KDTree(cur)
    for i in range(kround):
        print float(i)/kround,"begin"
        close=[]
        for j in range(len(tar)):
            index = kdtree.query(newtar[j])[1]
            close.append(cur[index])
        closearray= np.array(close)
        #closearray=closeorder(newtar, cur)
        #closearray=cur
        qR, rotation, trans = getReg(tararray, closearray)

        #a=tararray.dot(rotation.T)
        #T= applied(trans, T)
        #newtar = T* tararray
        newtar =  add(tararray.dot(rotation), trans) #apply Reg
        #newtar =  add(tararray, trans)
        newmse = MSE(newtar, closearray)
        if abs(newmse - mse) < threshold:
            break
        mse = newmse

        #threeDvision(newtar,cur)

    threeDvision(newtar,cur)
    return rotation, trans.reshape(3,1)


def main():
    #rfile: probe data   rfile2: link data
    rfile="/Users/zhanglida/Downloads/point_cloud_registration/point_cloud_registration/pointcloud1.fuse"
    rfile2="/Users/zhanglida/Downloads/point_cloud_registration/point_cloud_registration/pointcloud2.fuse"
	
    #read cloud1
    cloud1 = open(rfile).readlines()
    print "Load Cloud1"
    #read cloud2
    cloud2 = open(rfile2).readlines()
    print "Load Cloud2"
    
    tar=[]
    for line in cloud1[:50000]:
        nextline = []
        lines = line.split(" ")
        nextline.append(float(lines[0]))
        nextline.append(float(lines[1]))
        nextline.append(float(lines[2]))
        #nextline.append(float(1))
        tar.append(nextline)
    tar= np.array(tar)
    
    cur=[]
    for line in cloud2[:50000]:
        nextline = []
        lines = line.split(" ")
        nextline.append(float(lines[0]))
        nextline.append(float(lines[1]))
        nextline.append(float(lines[2]))
        #nextline.append(float(1))
        cur.append(nextline)
    cur= np.array(cur)
    #threeDvision(tar,cur)
    
    rotation, trans = ICP(tar, cur)
    #print len(tar),len(cur)
    #print rotation
    #print trans
    
    
    
if __name__ == "__main__":
	main()
