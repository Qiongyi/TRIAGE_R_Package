#pip install numpy matplotlib pandas

import numpy as np
import numpy.matlib as matlib
import pandas as pd
import matplotlib.pyplot as plt

def find_elbow(input_list):
    # takes a list of values and finds the elbow point (as index)    
    # in reverse order 
    curve = input_list
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec ** 2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(vecFromFirst * matlib.repmat(lineVecNorm, nPoints, 1), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint

filename = 'human_rts_epimap.txt'
data = pd.read_csv(filename, sep='\t', header=None)
values = data[1].values
elbow_index = find_elbow(values)

plt.figure(figsize=(10, 6))
plt.plot(range(len(values)), values, marker='o', linestyle='-', color='b', markersize=2, label='Data Points')
plt.scatter(elbow_index, values[elbow_index], color='r', s=100, label=f'Elbow Point (Index: {elbow_index})')
plt.title('Elbow Point Detection')
plt.xlabel('Index')
plt.ylabel('RTS Values')
plt.legend()
plt.grid(True)
#plt.show()

plt.annotate(f'Elbow Point\nIndex: {elbow_index}',
             xy=(elbow_index, values[elbow_index]),
             xytext=(elbow_index + 5, values[elbow_index] + 0.05),
             arrowprops=dict(facecolor='black', arrowstyle='->', connectionstyle='arc3,rad=0.2'))

plt.savefig('elbow_point_detection.png')
plt.savefig('elbow_point_detection.svg')
plt.savefig('elbow_point_detection.pdf', format='pdf')
