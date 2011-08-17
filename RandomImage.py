import numpy as np


def RandomImage(pointCnt, maxRadius, maxPower, domainShape) :
    xPos = np.random.randint(0, domainShape[0], pointCnt)
    yPos = np.random.randint(0, domainShape[1], pointCnt)
    theRadius = np.random.randint(1, maxRadius + 1, pointCnt)
    thePower = np.random.rand(pointCnt) * maxPower

    dataVals = np.zeros(domainShape, dtype='i')

    for x in xrange(domainShape[0]) :
        for y in xrange(domainShape[1]) :
            dists = np.sqrt((x - xPos)**2 + (y - yPos)**2)
            dataVals[x, y] = int(np.sum(thePower * np.exp(-(dists**2) / 
                                                    (2*(theRadius/2.0)**2))))

    return dataVals

