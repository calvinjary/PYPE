# E:\Dropbox\9 PhD Research\1 PyPE
# C:\Users\Calvin Jary\Dropbox\9 PhD Research\1 PyPE
# C:\Research 

'''
Welcome to PyPE by Calvin Jary
Features include reporting progress along with timestamps
Writing out results dynamically
'''

# Importing the libraries, note numpy is faster than pandas for this application
import pandas as pd
import numpy as np
import struct
import os
from time import localtime, strftime

#################################################### User configurable variables  ###############################################################

# name of directory which contains PPI matrix outputs in the format: Q16576 to Q8N680.txt   this will also be the name of the output text file
#directory = 'groundtruth'

# name of output file
outputname = 'G1'

# Class label, such as: GroundTruth, Random, or some other class label? Recommend to keep it to 1 character
label = 'G'

###################################################### PyPE main algorithm  ######################################################################

if label == 'G':
    try:
        directory = 'groundtruth'
        pdRandomPairsList = pd.read_csv('RandomPairsGT.txt', sep = '\t', header = None, engine='c')
    except:
        print('could not read RandomPairsGT.txt!')
        
elif label == 'R':
    try:
        directory = 'randompairs'
        pdRandomPairsList = pd.read_csv('RandomPairsRP.txt', sep = '\t', header = None, engine='c')
    except:
        print('could not read RandomPairsRP.txt!')
     
#Results = ([],) * 205   sadly does not work! must instead use:
Results = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], 

LoopCounter = int(-1)
ProteinA_Name, ProteinB_Name = '', ''


while LoopCounter < 10000:   #10000
    LoopCounter += 1
    
    #for filename in os.listdir(directory):       
    ProteinA_Name = pdRandomPairsList.iloc[LoopCounter, 0]
    ProteinB_Name = pdRandomPairsList.iloc[LoopCounter, 1]
    
        
#######################################  reading in the matrix file and solvent accessiblity files  ############################################        
        
    # the current matrix file
    MatrixOutputFile = ProteinA_Name + ' to ' + ProteinB_Name + '.txt'
    
    # reading in the matrix output file
    try:
        pdMatrixDataFrame = pd.read_csv(directory + '/' + MatrixOutputFile, sep = '\s+', header = None, engine='c')
    except:
        print('could not read matrix file: ', MatrixOutputFile)
        continue
     
    # convert pandas dataframe to numpy array to improve performance
    # can also use dataframe = dataframe.as_matrix()
    npMatrix = pdMatrixDataFrame.values
    MatrixShape = npMatrix.shape
    
    # reading in solvent accessbility for protein A  
    try:
        pdProteinA_NSP = pd.read_csv('ProtDCal/' + ProteinA_Name + '.txt', sep = '\t', header = None, engine='c')
    except:
        print('could not read solvent accessibility file for protein: ', ProteinA_Name)
        continue
        
    # checking to make sure protein A solvent accessiblity table entries are complete
    if len(pdProteinA_NSP) != pdProteinA_NSP.iloc[len(pdProteinA_NSP)-1,1]:
        print('the solvent accessiblity file for Protein ', ProteinA_Name, ' is corrupted and has insufficient entries!')
        continue
                
    # reading in solvent accessbility for protein B  
    try:
        pdProteinB_NSP = pd.read_csv('ProtDCal/' + ProteinB_Name + '.txt', sep = '\t', header = None, engine='c')
    except:
        print('could not read solvent accessibility file for protein: ', ProteinB_Name)
        continue

    # checking to make sure protein B solvent accessiblity table entries are complete
    if len(pdProteinB_NSP) != pdProteinB_NSP.iloc[len(pdProteinB_NSP)-1,1]:
        print('the solvent accessiblity file for Protein ', ProteinB_Name, ' is corrupted and has insufficient entries!')
        continue
            
    
#################################################  reading GenTab file for Protein A  ##########################################################
  
    try:
        f = open('gentab/' + ProteinA_Name, 'rb')
    except:
        print('could not read GenTab file for protein: ', ProteinA_Name)
        continue        
            
    buf = f.read()
    gen = struct.iter_unpack('i', buf)  # could be 'l'
            
    array = list(gen)
    out = [array[0][0]]
    index = 1
                
    while index < len(array):
        out.append(array[index][0])
        index = index + 1 + array[index][0]
    #print(out)
                
    ProteinAGenTab = []   
                
    for i in range(1, (len(out)-19)):
        ProteinAGenTab.append(out[i])
         
################################################  reading GenTab file for Protein B  #############################################################

    try:
        f = open('gentab/' + ProteinB_Name, 'rb')
    except:
        print('could not read GenTab file for protein: ', ProteinB_Name)
        continue          
        
    buf = f.read()
    gen = struct.iter_unpack('i', buf)  # could be 'l'
            
    array = list(gen)
    out = [array[0][0]]
    index = 1
                
    while index < len(array):
        out.append(array[index][0])
        index = index + 1 + array[index][0]
    #print(out)
                
    ProteinBGenTab = []   
                
    for i in range(1, (len(out)-19)):
        ProteinBGenTab.append(out[i])    
    
###############################################  Reporting to the user that all files have been loaded in #########################################
    
    if LoopCounter == 0:
        print('PyPE has successfully initialized at time: ')
        print(strftime("%Y-%m-%d %H:%M:%S", localtime()))    
    
################################################  Solvent Accessbility Protein A ###################################################################
    
    # Step 1: Ingest both Relative Solvent Accessiblity and Absolute Solvent Accessiblity outputs to then Modulate
    npProteinA_Mod = np.zeros((len(pdProteinA_NSP),166))
    
    for column in range(166):
        npProteinA_Mod[:,column] = pdProteinA_NSP.iloc[:,column+3] # setting up relative solvent accessiblity as coloumns 0 to 4 inclusive        
     
    '''
        if column >= 3:
            npProteinA_Mod[:,column] = pdProteinA_NSP.iloc[:,5] # setting up absolute solvent accessiblity as coloumns 5 to 9 inclusive
    
    # Step 2: Choosing a Modulation Method
    for column in range(5):
        if column == 0:    # Raw values
            npProteinA_Mod[:,column] = npProteinA_Mod[:,column] 
        if column == 1:     # square rooting 
            npProteinA_Mod[:,column] = np.sqrt(npProteinA_Mod[:,column]) 
        if column == 2:    # Cube rooting 
            npProteinA_Mod[:,column] = np.cbrt(npProteinA_Mod[:,column]) 
        # absolute SA
        if column == 3:    # Raw values
            npProteinA_Mod[:,column] = npProteinA_Mod[:,column] 
        if column == 4:    # square rooting
            npProteinA_Mod[:,column] = np.sqrt(npProteinA_Mod[:,column]) 
    
    '''
    
    # Step 3: Choosing a Window filtering method
    npProteinA_WF = np.zeros((len(pdProteinA_NSP)-19,166))   # - 19 because a 20 amino acid window is summarized by 1 value
    
    # implementing mean window filtering on the first 10 SA modulations
    for column in range(0, 166, 1):   
        for row in range(len(npProteinA_WF)):  
            array = np.array(npProteinA_Mod[row:row+20,column])      # setting up a 20 number long array (without using a loop)
            npProteinA_WF[row, column] = np.mean(array)              # taking the mean value of the array and storing that  
                
            
    # Taking the absolute value
    npProteinA_WF = np.absolute(npProteinA_WF)
    # This is the normalization step
    npProteinA_WF = npProteinA_WF / npProteinA_WF.max(axis=0)
    
    #from sklearn.preprocessing import normalize
    #npProteinA_WF = normalize(npProteinA_WF, axis=0, norm='max', copy=False )    

###############################################  Solvent Accessbility Protein B ##################################################################
    
    # Step 1: Ingest both Relative Solvent Accessiblity and Absolute Solvent Accessiblity outputs to then Modulate
    npProteinB_Mod = np.zeros((len(pdProteinB_NSP),166))
    
    for column in range(166):
        npProteinB_Mod[:,column] = pdProteinB_NSP.iloc[:,column+3] # setting up relative solvent accessiblity as coloumns 0 to 4 inclusive        

    '''
    for column in range(5):
        if column <= 2:
            npProteinB_Mod[:,column] = pdProteinB_NSP.iloc[:,4] # setting up relative solvent accessiblity as coloumns 0 to 4 inclusive        
        if column >= 3:
            npProteinB_Mod[:,column] = pdProteinB_NSP.iloc[:,5] # setting up absolute solvent accessiblity as coloumns 5 to 9 inclusive
    
    # Step 2: Choosing a Modulation Method
    for column in range(5):
        if column == 0:    # Raw values
            npProteinB_Mod[:,column] = npProteinB_Mod[:,column] 
        if column == 1:    # square rooting
            npProteinB_Mod[:,column] = np.sqrt(npProteinB_Mod[:,column]) 
        if column == 2:    # Cube rooting 
            npProteinB_Mod[:,column] = np.cbrt(npProteinB_Mod[:,column]) # alternative code npProteinA_SA[:,i] * npProteinA_SA[:,i] * npProteinA_SA[:,i]
        # absolute SA
        if column == 3:    # Raw values
            npProteinB_Mod[:,column] = npProteinB_Mod[:,column] 
        if column == 4:    # square rooting
            npProteinB_Mod[:,column] = np.sqrt(npProteinB_Mod[:,column]) 
    
    '''
    
    # Step 3: Choosing a Window filtering method
    npProteinB_WF = np.zeros((len(pdProteinB_NSP)-19,166))   # - 19 because a 20 amino acid window is summarized by 1 value
    
    # implementing mean window filtering on the first 10 SA modulations
    for column in range(0, 166, 1):   
        for row in range(len(npProteinB_WF)):  
            array = np.array(npProteinB_Mod[row:row+20,column])      # setting up a 20 number long array (without using a loop)
            npProteinB_WF[row, column] = np.mean(array)              # taking the mean value of the array and storing that  
                
  
    # Taking the absolute value
    npProteinB_WF = np.absolute(npProteinB_WF)
    # This is the normalization step
    npProteinB_WF = npProteinB_WF / npProteinB_WF.max(axis=0)    
    
    
    #from sklearn.preprocessing import normalize
    #npProteinA_WF = normalize(npProteinA_WF, axis=0, norm='max', copy=False ) 

########################################  computing the conventional (outdated) PIPE score  ######################################################
    '''
    #   Conventional PIPE score
    # - there is a 3x3 matrix centered on you. if 5 or more of these 9 elements are nonzero, you get to be a 1. otherwise you are a 0.
    # - this means if your raw matrix score is 0, but 5 or more of your 8 possible neighbors are nonzero, you get to be a 1
    # - if you are on the border, missing neighbours are a 0, but you can still be 1 if you assemble 5 or more nonzero elements
    # - Final PIPE score is calculated as how many matrix elements were 1's divded by total number of matrix elements
    
    # creating a new matrix that has a border of 0's to account for iterative scoring
    npNewMatrix = np.zeros((MatrixShape[0]+2,MatrixShape[1]+2))
    for i in range(MatrixShape[0]):
        for j in range(MatrixShape[1]):
            npNewMatrix[i+1, j+1] = npMatrix[i,j]
    
    NonZeroElement, OneScoreCounter, ConventionalPipeScore = int(0), int(0), float(0)
    
    # looping through all the elements of interest
    for i in range(1, MatrixShape[0]+1):  # iterating through all the rows
        for j in range(1, MatrixShape[1]+1): # iterating through all the coloums
            NonZeroElement = 0
            # determing total number of nonzero elements in 3x3 matrix centered on element of interest
            for k in range(i-1, i+2, 1):
                for l in range(j-1, j+2, 1):
                    if npNewMatrix[k,l] > 0:
                        NonZeroElement += 1
            if NonZeroElement >= 5:   # npNewMatrix[i,j] > 0
                OneScoreCounter += 1
            
    ConventionalPipeScore = OneScoreCounter / ((MatrixShape[0]) * (MatrixShape[1]))
    '''        
#################################################  computing Sim-weighted score  ##################################################################
    '''
    IntTotalScore, FloatTotalScore, NowScore = float(0), float(0), float(0)
    IntScore = np.int16(0)
    FloatSimWeightedScore, IntSimWeightedScore = float(0), float(0)
    
    maxH = int(0)
    
    for i in range((MatrixShape[0])):  # iterating through Protein A
        for j in range((MatrixShape[1])):  # iterating through Protein B
            # Calculating sim weighted score for this particular element
            NowScore = 1000 * ((npMatrix[i,j]) / (ProteinAGenTab[i] * ProteinBGenTab[j]))  
            
            # Accurate floating point section
            FloatTotalScore += NowScore   # (ProteinA[i] * ProteinB[j])
            
            
            # Inaccurate integer rounding error section
            #IntScore = np.int16(NowScore)  # introducing integer rounding errors
            #IntTotalScore += IntScore
            
            # Calculating matrix max (highest peak of the entire matrix)
            #if (npMatrix[i,j]) > maxH:
            #    maxH = (npMatrix[i,j])
    
    # Calculating accurate Sim Weighted score by averaging over total number of matrix elements
    FloatSimWeightedScore =  FloatTotalScore / ((MatrixShape[0]) * (MatrixShape[1]))
    
    
    # Calculating inaccurate Sim Weighted integer score by averaging over total number of matrix elements
    #IntSimWeightedScore = IntTotalScore / ((MatrixShape[0]) * (MatrixShape[1]))
    '''
###############################################  200 solvent accessiblity approaches  ###########################################################
 
    RunningTotal = np.zeros(200)
    ConventionalSimWeightedScore = float(0)
       
    for i in range((MatrixShape[0])):  # iterating through Protein A
        for j in range((MatrixShape[1])): # iterating through Protein B   
            
            # conventional sim-weighted score
            
            try:
                ConventionalSimWeightedScore = 1000 * (npMatrix[i,j] / (ProteinAGenTab[i] * ProteinBGenTab[j]))
                RunningTotal[0] += ConventionalSimWeightedScore
            except:
                continue
            
            for w in range(166):  # should be 166
                try:
                    RunningTotal[w+1] += ConventionalSimWeightedScore * (4*(    max(npProteinA_WF[i,w], npProteinB_WF[j,w]) / min(npProteinA_WF[i,w], npProteinB_WF[j,w])      ))             
                except:
                    continue
  
    FinalResults = RunningTotal / (MatrixShape[0] * MatrixShape[1])
    
##################################################  tabulating and looping  #####################################################################    

    # Tabulating Results
    Results[0].append(str(ProteinA_Name))
    Results[1].append(str(ProteinB_Name))
    Results[2].append(str(label)) # class label configurable by the user
    
    for i in range(180):
        Results[i+3].append(str(FinalResults[i]))  #str
          
    ResultsString = [", ".join(map(str, x)) for x in Results]
    
    Results_PandasDataframe = pd.DataFrame(ResultsString)

    Results_PandasDataframe.to_csv(outputname + '.txt', sep = '\t', index = False, index_label = False)    
    
    '''
    Results_NumPyArray = (np.array(dude))
    Results_PandasDataframe = pd.DataFrame(Results_NumPyArray)
    Results_NumPyArray = (np.asarray(Results, dtype = str))
    '''
    
    # reporting
    if (LoopCounter % 100) == 0:
        print((LoopCounter / 100), '% of PPIs successfully calculated at time: ')
        print(strftime("%Y-%m-%d %H:%M:%S", localtime()))
    
    
######################################################  tabulating results  ####################################################################    
    

print('PyPE has successfully completed at time: ')
print(strftime("%Y-%m-%d %H:%M:%S", localtime()))



