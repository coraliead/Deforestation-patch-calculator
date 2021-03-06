import numpy as np
import xarray as xr

filepath = '/home/coralie/bash_project/Africa/'

ForestLossXr = xr.open_dataarray(filepath + 'ForestLoss2000-19 lon15-16 lat5-6.nc')

data = ForestLossXr.values
ForestLossXr = []

numOfRows = np.shape(data)[0]
numOfCols = np.shape(data)[1]

#create a matching 2d array of -1s to track where we have checked
mask = np.zeros_like(data, dtype = 'uint64')

#start iterating at 0,0
row = 0
col = 0
numOfClump = 0
cellsToLookAt = []
deforestYr = 16

while row < numOfRows and col < numOfCols:
  #this is tests if that cell has been deforested and whether it has already been accounted for
  if data[row][col] == deforestYr and mask[row][col] == 0:
    numOfClump += 1
    
    #create an empty queue of items to look at
    cellsToLookAt = [[row, col]]
    
    #these below 10 lines set the limiters for possible adjacencies - 
    #normally 9, but 4 if a corner and 6 if an edge
    minX, maxX = -1, 2
    minY, maxY = -1, 2
    
    if row == 0:
      minX = 0
    elif row == numOfRows - 1:
      maxX = 1
    
    if col == 0:
      minY = 0
    elif col == numOfCols - 1:
      maxY = 1
    
    #first pass will have the original cell, will continue until empty
    while cellsToLookAt:
      
      thisItem = cellsToLookAt.pop(0)
      cItemX = thisItem[0]
      cItemY = thisItem[1]
      if cItemY == numOfCols - 1: continue
      if cItemX == numOfRows - 1: continue
      
      for x in range(minX, maxX):
        for y in range(minY, maxY):
          #test the adjacent 8 - if we enter then one of the surrounding is deforested
          if data[cItemX+x][cItemY+y] == deforestYr and mask[cItemX+x][cItemY+y] == 0:
            cellsToLookAt.append([cItemX+x, cItemY+y])

            mask[cItemX+x][cItemY+y] = numOfClump

         
  if col == numOfCols-50: print(row, col)
  if col <= numOfCols - 2:
    col += 1

  else:
    row += 1
    col = 0
#%%
# this section is calculating the area of the clumps
PixelSize = 30*30
clumpSizeStorage = np.zeros([numOfClump], dtype = 'float64')

# clump no = 13,5535 so 1355 times it will print the clump number

for clump in range(1,numOfClump):
    clumpCells = np.size(mask[mask == clump])
    clumpSizeStorage[clump] = clumpCells * PixelSize
    
    if clump % 100 == 0:
        print(clump)
        
# ClumpArr is an array of each clump's size in m2

# setting 0's to nan to not skew the data. I think the zeros are present bc they are cells which had no adjacency
# and so should actually be 900m. however, need to check that this is correct before changing the algorithm

clumpSizeStorage[clumpSizeStorage == 0] = np.nan
# convert m2 to km2
ClumpArr_km = clumpSizeStorage / 1000000
# this step is only neccessary if i haven't redone the clumpSizeStorage bit again (the loop used to start at 0 and
# so it counted zero as a clump)
ClumpArr_km1 = ClumpArr_km[1:]

SizeMean = np.nanmean(ClumpArr_km1)
SizeMed = np.nanmedian(ClumpArr_km1)

SizeSD = np.nanstd(ClumpArr_km1)
Size = np.size(ClumpArr_km1)

SizeSE = SizeSD / np.sqrt(Size)
z_score = 1.96  # this is 95%
CI = SizeSE * z_score
CI1 = SizeMean - CI
CI2 = SizeMean + CI
