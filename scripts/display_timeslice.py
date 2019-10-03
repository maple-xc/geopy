
import sys 
sys.path.append(".") 

import numpy as np
import os

import matplotlib.pyplot as plt

from segpy.reader import create_reader, read_binary_reel_header, read_trace_header
from segpy.trace_header import TraceHeaderRev1
from segpy.packer import HeaderPacker
from seismic.inputoutput import readSeis2DMatFromSegyNoInfo

from seismic.analysis import analysis as seis_ays
# from segpy.packer import make_header_packer


NUMPY_DTYPES = {'ibm':     np.dtype('f4'),

                'int32':   np.dtype('i4'),

                'int16':   np.dtype('i2'),

                'float32': np.dtype('f4'),

                'int8':    np.dtype('i1')}





def make_dtype(data_sample_format): # TODO: What is the correct name for this arg?

    """Convert a SEG Y data sample format to a compatible numpy dtype.



    Note :

        IBM float data sample formats ('ibm') will correspond to IEEE float data types.



    Args:

        data_sample_format: A data sample format string.



    Returns:

        A numpy.dtype instance.



    Raises:

        ValueError: For unrecognised data sample format strings.

    """

    try:

        return NUMPY_DTYPES[data_sample_format]

    except KeyError:

        raise ValueError("Unknown data sample format string {!r}".format(data_sample_format))



def extract_timeslice(segy_filename, slice_index, dtype= None, null=0):

    """Extract a timeslice from a 3D SEG Y file and return it as a Numpy array.



    Args:

        segy_filename: Filename of a SEG Y file.




        slice_index: The zero-based index (increasing with depth) of the slice to be extracted.



        dtype: Optional Numpy dtype for the result array. If not provided a dtype compatible with

            the SEG Y data will be used.



        null: Optional sample value to use for missing or short traces. Defaults to zero.

    """

    with open(segy_filename, 'rb') as segy_file:



        segy_reader = create_reader(segy_file)

        print(segy_reader.xline_range())
        
        

        if dtype is None:
    
            dtype = make_dtype(segy_reader.data_sample_format)

        if segy_reader.dimensionality != 3:

            raise DimensionalityError("Cannot slice {n} dimensional seismic.".format(segy_reader.dimensionality))



        i_size = segy_reader.num_inlines()

        x_size = segy_reader.num_xlines()

        t_size = segy_reader.max_num_trace_samples()



        if not (0 <= slice_index < t_size):

            raise ValueError("Time slice index {0} out of range {} to {}".format(slice_index, 0, t_size))



        timeslice = np.full((i_size, x_size), null, dtype)



        for inline_num, xline_num in segy_reader.inline_xline_numbers():

            trace_index = segy_reader.trace_index((inline_num, xline_num))

            trace = segy_reader.trace_samples(trace_index)



            try:

                sample = trace[slice_index]

            except IndexError:

                sample = null



            i_index = segy_reader.inline_range().index(inline_num)

            x_index = segy_reader.xline_range().index(xline_num)



            timeslice[i_index, x_index] = sample



        return timeslice





def nullable_dtype(s):

    return None if s == "" else np.dtype(s)


def normalized(self, data):
    
    arr = data.copy()
    arr = arr - arr.mean(axis=0)
    arr = arr / np.abs(arr).max(axis=0)

    return arr  








def plotSeisZSliceFrom2DMat(seis2dmat, zsls=None, datacol=3,

                            inlcol=0, xlcol=1, zcol=2,

                            colormap= 'seismic',                            

                            titlesurf='', colorbaron=True,

                            verbose=False):

    """

    Plot seismic z slices from 2D matrix
    Codes modified from GeoPy: https://github.com/geopyteam/geopy/blob/master/src/seismic/visualization.py

    Args:

        seis2dmat:  2D matrix representing seismic data

                    It contains at least four columns, [IL, XL, Z, Value, ...]

        zsls:       depth/time No. for plotting

                    Plot all z slices if not specified

        datacol:    index of data column for plotting in 2D matrix (indexing from 0)

                    Plot the fourth column if not specified (3)

        inlcol:     index of inline column. Default is the first column (0)

        xlcol:      index of crossline column. Default is the second column (1)

        zcol:       index of z column. Default is the third column (2)

        colormap:   colormap name for seismic data visualization, such as 'seismic'

                    Use the default colormap by vis_cmap.makeColorMap if not specified

        

        titlesurf:  surfix for the title. Default is blank

        colorbaron: colorbar display. Default is false

        verbose:    flag for message display. Default is True

    Return:

        None

    Note:

        Negative z is used in the vertical direction

    """



    # Check input matrix

    if np.ndim(seis2dmat) != 2:

        print('ERROR in plotSeisZSliceFrom2DMat: 2D seismic matrix expected')

        return

    if datacol < 0 or len(seis2dmat[0, :]) <= datacol:

        print('ERROR in plotSeisZSliceFrom2DMat: Not data column found in 2D seismic matrix')

        return

    if inlcol < 0 or len(seis2dmat[0, :]) <= inlcol:

        print('ERROR in plotSeisZSliceFrom2DMat: Not inline column found in 2D seismic matrix')

        return

    if xlcol < 0 or len(seis2dmat[0, :]) <= xlcol:

        print('ERROR in plotSeisZSliceFrom2DMat: Not crossline column found in 2D seismic matrix')

        return

    if zcol < 0 or len(seis2dmat[0, :]) <= zcol:

        print('ERROR in plotSeisZSliceFrom2DMat: Not z column found in 2D seismic matrix')

        return



    seisinfo = seis_ays.getSeisInfoFrom2DMat(seis2dmat,

                                            inlcol=inlcol, xlcol=xlcol, zcol=zcol)

    seis3dmat = seis_ays.convertSeis2DMatTo3DMat(seis2dmat,

                                                datacol=datacol,

                                                inlcol=inlcol, xlcol=xlcol, zcol=zcol)



    inlrange = seisinfo['ILRange']

    xlrange = seisinfo['XLRange']

    zrange = seisinfo['ZRange']

    inlstart = seisinfo['ILStart']

    inlend = seisinfo['ILEnd']

    xlstart = seisinfo['XLStart']

    xlend = seisinfo['XLEnd']

    zstart = seisinfo['ZStart']

    zstep = seisinfo['ZStep']

    znum = seisinfo['ZNum']

    if znum == 1:

        zstep = -1



    if zsls is None:

        print('WARNING in plotSeisZSliceFrom2DMat: to plot all a slices in 2D seismic matrix')

        zsls = zrange



    if np.ndim(zsls) != 1:

        print('ERROR in plotSeisZSliceFrom2DMat: 1D array of z slices expected')

        return



    x, y = np.meshgrid(xlrange, inlrange)



    nzsls = len(zsls)

    if verbose:

        print('Plot ' + str(nzsls) + ' z slices')

    for idx in zsls:        

        # idx = np.round((z - zstart) / zstep).astype(np.int32)

        if idx >= 0 and idx < znum:

            seisdata = seis3dmat[idx, :, :]

            seisdata = seisdata.transpose()

            valuemin = seisdata.min()
            valuemax = seisdata.max()

            plt.figure(facecolor='white')

            plt.pcolormesh(x, y, seisdata,

                           cmap=colormap,

                           shading='gouraud',

                           vmin=valuemin, vmax=valuemax)

            plt.xlim([xlstart, xlend])

            plt.ylim([inlstart, inlend])

            plt.title('Depth/Time at ' + str(zrange[idx]) + titlesurf)

            plt.xlabel('Crossline No.')

            plt.ylabel('Inline No.')

            if colorbaron:

                plt.colorbar()

    plt.show()



    return

def main():

    
    segy_file = 'C:/Users/ZHAO/Desktop/Programming-Geophysics-in-Python/Datasets/F3/original_seismic_small_t1500-1512.sgy'
    
    
    slice_index = 0

    if not os.path.exists('timeslice.npy'):
        timeslice = extract_timeslice(segy_file, slice_index)
        np.save('timeslice.npy', timeslice)

    else:

        timeslice = np.load('timeslice.npy')



    # print(timeslice.shape)

    plt.imshow(np.flipud(timeslice), cmap= 'seismic')
    plt.colorbar()
    plt.show()

    if not os.path.exists('seis2dmat.npy'):
        seis2dmat = readSeis2DMatFromSegyNoInfo(segy_file)
        np.save('seis2dmat.npy', seis2dmat)

    else:

        seis2dmat = np.load('seis2dmat.npy')

    
    plotSeisZSliceFrom2DMat(seis2dmat, zsls = np.array([slice_index]))






if __name__ == '__main__':
    
    main()
    
