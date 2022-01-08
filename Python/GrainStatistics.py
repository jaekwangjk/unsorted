import numpy as np 
import matplotlib.pyplot as plt

class GrainStatistics:
     #Initialization
    def __init__(self, areaFile, ADJFile, lcount):
        self.area = np.loadtxt(areaFile)
        self.ADJ = np.loadtxt(ADJFile,delimiter=',')
        
        self.side = np.zeros(lcount)
        for i in range (0,5000):
            self.side[i] = sum(self.ADJ[i,:])
 
        extract_ID = np.nonzero(self.area)
        self.area =self.area[extract_ID]
        self.side =self.side[extract_ID]
        self.nAlive = self.area.shape[0]
        
    def sideStatistics(self):

        bins = [2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,10.5,11.5] 
        x=[3,4,5,6,7,8,9,10,11]
        count_hist, _ = np.histogram(self.side, bins)
        y = count_hist/self.nAlive
        
        return x,y 
    
    def areaStatistics(self):

        bins = [0,0.5,1.5,2.5,3.5,4.5,5.5,6.5] 
        x=[0.25,1,2,3,4,5,6]
        count_hist, _ = np.histogram(self.area/np.mean(self.area), bins)
        y= count_hist/self.nAlive
        
        return x,y

