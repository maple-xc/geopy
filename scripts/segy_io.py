
import sys 
sys.path.append(".") 

from segpy.reader import create_reader

def main():

    

    segy_in_file = 'C:/Users/ZHAO/Desktop/Programming-Geophysics-in-Python/Datasets/F3/original_seismic_small_t1500-1512.sgy'

    segy_file = open(segy_in_file, 'rb')

    segy_reader = create_reader(segy_file, endian='>')

    
    tid = 0
    
    print("=== BEGIN TRACE #1 HEADER ===")

    header = segy_reader.trace_header(tid)

    print('Inline Number： ', header.inline_number)
    print('Crossline Number： ', header.crossline_number)

    print("=== END TRACE #1 HEADER ===")

    segy_file.close()

if __name__ == '__main__':
    
    main()
    
