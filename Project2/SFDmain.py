import numpy as np
import matplotlib.pyplot as plt


def convertRepresentation(image, N, M):
    if len(image.shape) > 1: #if the image is in a matrix representation
        vectorRep = np.zeros(N*M)
        for i in range(N):
            for j in range(M):
                vectorRep[M*i + j] = image[i, j]
        return vectorRep
    elif len(image.shape) == 0:
        matrixRep = np.zeros((N, M))
        for i in range(N*M):
            matrixRep[, ] = image[i]