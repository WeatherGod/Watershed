
import numpy as np
from itertools import chain

# Also depends on scipy.ndimage for a caching mechanism

def Watershed_Transform(image, maxDepth, saliency) :
    centers = _find_centers(image)

    return _find_basins(image, centers, maxDepth, saliency)




class GreyLevel_Blob :
    def __init__(self) :
        self.polarity = True
        self.scale_level = None
        self.extremum = None
        self.saddle = None
        self.support_region = None



class Extremum_Region :
    def __init__(self, pos=(), val=None, greyblob=None) :
        self.position = pos
        self.grey_val = val
        self.grey_blob = greyblob

#   def isNeighbor(self, pos) :
#       return(min([(pos[0] - aPoint[0])**2 + (pos[1] - aPoint[1])**2
#               for aPoint in self.positions]) <= 1.0)
#
#   def addPoint(self, pos) :
#       self.positions.append(pos)

class Saddle_Region :
    def __init__(self, pos=(), val=None, greyblobs=()) :
        self.position = pos
        self.grey_val = val
        self.grey_blobs = greyblobs


class Support_Region :
    def __init__(self) :
        self.blob_area = None
        self.first_moment = None
        self.second_moment = None
        self.pixels = None
        self.atBoundary = None

GLOBBED = -3
UNMARKED = -1


def _neighbors(pos, rangeVal, shape, inclSelf=False) :
        startx = max(pos[0] - rangeVal, 0)
        starty = max(pos[1] - rangeVal, 0)
        endx = min(pos[0] + rangeVal + 1, shape[0])
        endy = min(pos[1] + rangeVal + 1, shape[1])
        if inclSelf :
                return [(x, y) for x in range(startx, endx)
                               for y in range(starty, endy)]
        else :
                return [(x, y) for x in range(startx, endx)
                               for y in range(starty, endy) if
                        (x != pos[0] or y != pos[1])]

#   neighbors = []
#   for x in range(max(pos[0] - rangeVal, 0), 
#              min(pos[0] + rangeVal + 1, shape[0])) :
#       for y in range(max(pos[1] - rangeVal, 0), 
#                  min(pos[1] + rangeVal + 1, shape[1])) :
#           if (inclSelf or x != pos[0] or y != pos[1]) :
#               neighbors.append((x, y))
#
#   return neighbors


def _find_centers(image) :
    # Assume that image is a normalized, quantized image...
    pixels = [ [] for i in xrange(image.max() + 1)]
    centers = [ [] for i in xrange(image.max() + 1)]

    # loading *pixels* with coordinates of pixels with associated values
    xs = range(image.shape[0])
    ys = range(image.shape[1])
    for x in xs :
        for y in ys :
            pixels[image[x, y]].append((x, y))

    # Exclude any pixels that are not peaks.
    # Also, need to deal with degenerate cases of two neighboring
    # pixels having the same value.
    marked = np.zeros(image.shape, dtype=bool)
    for q in range(image.max(), -1, -1) :
        for p in pixels[q] :
            if not marked[p] :
                isCenter = False
                markedSoFar = []
                for point in _neighbors(p, 2, image.shape, inclSelf=True) :
                    if not marked[point] :
                        marked[point] = True
                        markedSoFar.append(point)
                        isCenter = True
                    else :
                        # p touches an already marked point,
                        # so it can't be a center, but also don't
                        # want to bother with it again.
                        #marked[p] = 1
                        isCenter = False
                        break


                if isCenter : 
                    centers[q].append(Extremum_Region(p, q))
                    marked[p] = True
                else :
                    # time to undo the markings
                    #for aPos in markedSoFar :
                    #    marked[aPos] = False
                    marked[zip(*markedSoFar)] = False

#       print q, len(centers[q])
    return centers



def _find_basins(image, centers, maxDepth, saliency) :
    # Reset caching
    _capture._markedSoFar = {}

    # Zero for background, -1 for unchecked, positive values for blob number
    basins = UNMARKED * np.ones(image.shape, dtype='i')

    # Initializing the basin number
    basinNumber = 1
    
    globs = []

    deferredToNext = []
    for level in range(len(centers) - 1, maxDepth - 1, -1) :
        # Hysterisis level.  Don't let it get below 0.
        hlevel = level - maxDepth
        print hlevel
        #print level, len(centers[level]), len(deferredToNext)

        centersTmp = centers[level] + deferredToNext
        deferredToNext = []

        foothills = []

        for centIndex, aCenter in enumerate(centersTmp) :
            if basins[aCenter.position] == UNMARKED :
                (basin, captured) = _capture(image, basins, aCenter,
                                             basinNumber, hlevel, saliency,
                                             foothills)
                if not captured :
                    # Defer to next iteration to see if it will get big enough
                    centersTmp[centIndex].grey_val -= 1
                    deferredToNext.append(centersTmp[centIndex])
                elif basin is not None :
                    globs.append(basin)
                    basinNumber += 1

            
        #print "%3d  Centers: %4d  Deferred: %3d  Globs: %4d  Foothills: %4d  " %  (level, len(centers[level]), len(deferredToNext), len(globs), len(foothills))
        _remove_foothills(image, basins, hlevel, centers, foothills)
    return globs, basins


def _remove_foothills(image, basins, hlevel, centers, foothills) :
    """
    Set points in the foothills as globbed (as opposed to any particular
    basin.
    This effectively removes the foothill points from consideration.
    """
    for (foothill, center) in foothills :
        binthresh = int(image[center.position] // 2)
        # Examine centers with peaks that have grey values of
        # half or less of the current center.
        centers_slice = centers[binthresh:]
        xs = np.array([otherCent.position[0] for otherCent in
                       chain(*centers_slice)])
        ys = np.array([otherCent.position[1] for otherCent in
                       chain(*centers_slice)])


        while len(foothill) > 0 :
            pixel = foothill.pop()
            if basins[pixel] != UNMARKED : continue
            basins[pixel] = GLOBBED
            is_closest = _is_closest(pixel, center, xs, ys)

            # Checking the neighbors
#            for point in _neighbors(pixel, 1, image.shape) :
#                if basins[point] == UNMARKED :
#                    if ( image[point] >= 0 and image[point] < hlevel
#                         and (image[point] <= image[pixel] or _is_closest(point, center, centers))) :
#                        foothill.append(point)
            pts = [point for point in _neighbors(pixel, 1, image.shape) if
                    (basins[point] == UNMARKED and
                     0 <= image[point] < hlevel and
                     (image[point] <= image[pixel] or is_closest))]
            foothill.extend(pts)

def _is_closest(pixel, center, xs, ys) :
    mydist = ((pixel[0] - center.position[0])**2 +
              (pixel[1] - center.position[1])**2)
#    mydist = np.hypot(pixel[0] - center.position[0],
#                      pixel[1] - center.position[1])



#    xs, ys = np.array(zip(*[(otherCent.position[0], otherCent.position[1]) for
#                            otherCent in chain(*centers[binthresh:])]))
#    dists = np.hypot(pixel[0] - xs, pixel[1] - ys)
    dists = ((pixel[0] - xs)**2 + (pixel[1] - ys)**2)
#    dists = np.hypot([pixel[0] - otherCent.position[0] for 
#                      otherCent in chain(*centers[binthresh:])],
#                     [pixel[1] - otherCent.position[1] for
#                      otherCent in chain(*centers[binthresh:])])
    return np.all(mydist <= dists)

#    for otherCenter in chain(*centers) :
#        if mydist > np.hypot(pixel[0] - otherCenter.position[0],
#                             pixel[1] - otherCenter.position[1]) :
#        if mydist > ((pixel[0] - otherCenter.position[0])**2 +
#                     (pixel[1] - otherCenter.position[1])**2) :
#            return False
#
#    return True


def _find_boundary(image) :
    from scipy.ndimage import binary_dilation, binary_erosion
    return image - binary_erosion(image)

def _capture(image, basins, center, basinNumber, hlevel, saliency, foothills) :


    foothill = []
    markedSoFar = _capture._markedSoFar.get(center.position, [])


    # Caching mechanism to reduce the amount of redundant
    # processing when a center is reconsidered
    if len(markedSoFar) == 0 :
        neighbors = [center.position]
    else :
        # Just to make sure I am not stepping on any other
        # basin's toes.
        tempy = np.zeros(image.shape, dtype=bool)
        tempy[zip(*markedSoFar)] = True
        edge = _find_boundary(tempy)

        # Rebuild markedSoFar to be set of points before the edge contour
        # this is a conservative caching mechanism in order to avoid assuming
        # that edge points has not been declared for someone else.
        markedSoFar = [tuple(pt) for pt in np.argwhere(tempy - edge) if
                       basins[tuple(pt)] == UNMARKED]
        del tempy
        neighbors = [tuple(pt) for pt in np.argwhere(edge) if
                     basins[tuple(pt)] == UNMARKED]

        basins[zip(*markedSoFar)] = basinNumber

    # Is this basin eligable to be considered again?
    # In other words, if the basin-growing did not
    # yeild a basin down to the hysteresis level that
    # was large enough, then it might be eligible to
    # be reconsidered at a lower hlevel.
    willBeConsideredAgain = False
    
    while len(neighbors) > 0 :
        pixel = neighbors.pop()
        if basins[pixel] != UNMARKED : continue # already processed

        basins[pixel] = basinNumber
        markedSoFar.append(pixel)

        # Checking the neighbors
        for point in _neighbors(pixel, 1, image.shape) :
            if (basins[point] == UNMARKED) :
                # Because the center's grey_val gets decremented when
                # reconsidered, this check sees if we have found any
                # possible change in the basin growing by finding a
                # neighbor point further down. If there is room to grow,
                # then the basin is eligible for reconsideration.
                if (not willBeConsideredAgain and
                    center.grey_val > image[point]) :
                    willBeConsideredAgain = True

                # Capturing all points that are above hlevel
                if image[point] >= hlevel :
                    neighbors.append(point)
                else :
                    # Too deep to be part of the basin, so it is a foothill
                    foothill.append(point)

    # If we are already at bottom, then it won't be re-considered
    if center.grey_val == 0 :
        willBeConsideredAgain = False

    # So, did this region grow large enough to be finished?
    bigEnough = (len(markedSoFar) >= saliency)

    basin = None
    
    if bigEnough :
        foothills.append((foothill, center))
        basin = GreyLevel_Blob()
        basin.extremum = center
        basin.support_region = Support_Region()
        basin.support_region.pixels = markedSoFar
        _capture._markedSoFar.pop(center.position, None)
    elif willBeConsideredAgain :
        # Basin has not been captured
        # Now I need to undo what I have done...
#        for p in markedSoFar :
#            basins[p] = UNMARKED
        basins[zip(*markedSoFar)] = UNMARKED
        _capture._markedSoFar[center.position] = markedSoFar
    else :
        #print "Not being deferred!"
        # So, it is not big enough, and it won't be considered again,
        # Then mark them as globbed, so they won't cause confusion.
#        for p in markedSoFar :
#            basins[p] = GLOBBED
        basins[zip(*markedSoFar)] = GLOBBED
        _capture._markedSoFar.pop(center.position, None)

    return (basin, (bigEnough or not willBeConsideredAgain))

# State variable for _capture() for caching.
_capture._markedSoFar = {}

if __name__ == '__main__' :
    #from RandomImage import RandomImage
    #np.random.seed(32)
    #i = RandomImage(200, 50, 95, (1000, 1000))

    from scipy.ndimage import imread, gaussian_filter
    i = imread("/home/bvr/SatData/2011.07.01.12.00.png")
    i = gaussian_filter(i, 1)[300:800, 300:i.shape[1]/2]
    print "Shape:", i.shape, "  DType:", i.dtype

    globs, basins = Watershed_Transform(i, 5, 250)
    print "Glob Cnt:", len(globs)
    #import cProfile
    #cProfile.run("Watershed_Transform(i, 40, 1000)", "watershed_lg5_profile")


    basins = np.ma.masked_array(basins, mask=(basins <= 0))
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    ax.imshow(i, cmap=plt.gray(), interpolation='none')
    ax = fig.add_subplot(1, 2, 2)
    ax.imshow(basins, vmin=0, cmap=plt.prism(), interpolation='none')
    plt.show()
